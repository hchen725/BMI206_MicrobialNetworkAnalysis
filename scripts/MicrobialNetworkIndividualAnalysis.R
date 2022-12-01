# Script used in BMI 206 Individual Project
# Analysis of a metabolic reaction network graph created from the adjacency matrices
# in Kim et al. (BMC Bioinformatics, 2019)

# ==== Setup ====
# Load packages used
library(here) # Paths
library(tidyverse) 
library(pbapply) # Apply function with progress bars

library(pheatmap) # For creating heatmaps
library(viridis) # Heatmap fill colors
library(igraph) # Node metric calculations


# ---- * Functions ----
# Used to calculate different node metrics from the adjacency matrix
source(here("scripts/MicrobialNetworkFunctions.R"))

# Other functions:

#' get_rxn_overlap
#'
#' @param mat List of adjacency matricies
#' @description Calculate the number of reactions shared between species
#' @return Tibble with all reactions and the number of species a reaction exists in
get_rxn_overlap <- function(mat){
  
  rxn_list <- lapply(names(mat), function(m){
    rxn_list <- data.frame(rxn = rownames(mat[[m]]),
                           exists = 1)
    colnames(rxn_list) <- c("rxn", m)
    return(rxn_list)
  })
  
  rxn_list <- rxn_list %>%
    reduce(full_join, by = "rxn")
  
  rxn_list <- rxn_list %>%
    replace(is.na(.), 0) %>%
    mutate(rxn_sum = rowSums(across(where(is.numeric)))) %>%
    arrange(desc(rxn_sum)) %>%
    tibble()
  
  return(rxn_list)
}

#' get_jaccard_index
#'
#' @param a object 1 to compare
#' @param b object 2 to compare
#' @desciption computes the jaccard index for two items
#' @return jaccard index value
get_jaccard_index <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  jaccard <- intersection/union
  return (jaccard)
}

#' sort_hclust
#'
#' @param ... Input passed to function as.dendrogram
sort_hclust <- function(...){
  as.hclust(dendsort::dendsort(as.dendrogram(...)))
}


# ---- * Aesthetics ----
# Because things that look nice look nice
theme_set(theme_bw())
# Color map used for metrics
cmap_metrics <- c("btw_centrality" = "#04a3bd",
                  "clustering_coefficient" = "#f0be3d",
                  "in-degree" = "#931e18",
                  "CascadeNum" = "#da7901",
                  "eigenvector_centrality" = "#247d3f",
                  "out-degree" = "#20235b")


# ---- * Paths ----
file_dir <- here("data/") # Directory containing adjacency matrix xlsx
output_dir <- here("output/") # Output folder for plots 


# ---- * Load adjacency matrix ----
adj_mat <- read_excel_allsheets(filename = here(file_dir, 
                                                "AdjacencyMatrix.xlsx"),
                                tibble = FALSE)

# Reformat adjacency matrix and format as matrices
adj_mat_names <- names(adj_mat)
adj_mat <- lapply(names(adj_mat), function(x){
  m <- adj_mat[[x]] %>%
    column_to_rownames("...1")
})
names(adj_mat) <- adj_mat_names





# ==== Reaction overlap ====
# Get reactions that overlap between species
rxn_list <- get_rxn_overlap(adj_mat)
table(rxn_list$rxn_sum)

# Get combinations of species and calculate the jaccard index
species_comp <- combn(names(adj_mat), 2) %>%
  t() %>%
  data.frame() 

for (i in 1:nrow(species_comp)){
  s1 <- species_comp$X1[i]
  s2 <- species_comp$X2[i]
  
  # Get number of shared reactions
  sub <- rxn_list %>%
    filter(rxn_sum != 1) %>%
    select(rxn, all_of(c(s1, s2))) %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    filter(sum == 2)
  num_rxn <- nrow(sub)
  
  # Get jaccard index between
  ji <- get_jaccard_index(a = names(adj_mat[[s1]]),
                          b = names(adj_mat[[s2]]))
  
  species_comp$shared[i] <- num_rxn
  species_comp$jaccard[i] <- ji
  
  # Species comparison models
  species_comp$model[i] <- paste0(s2, "~", s1)
}


# ---- * Plot jaccard index ----
# Visualize similarities between species purely by the jaccard index

# Reformat jaccard index to a symmetrical matrix for heatmap
jaccard_mat <- bind_rows(species_comp %>%
                          rename(X2_new = X1,
                                 X1 = X2) %>%
                          rename(X2 = X2_new) %>%
                          select(X1, X2, shared, jaccard),
                        species_comp %>%
                          select(-model)) %>%
  select(-shared) %>%
  pivot_wider(names_from = "X2",
              values_from = "jaccard",
              values_fill = 1) %>%
  column_to_rownames("X1") %>%
  as.matrix()

# Organize dendrogram clusters 
mat_cluster_cols <- hclust(dist(t(jaccard_mat)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
mat_cluster_rows <- hclust(dist(jaccard_mat))
mat_cluster_rows <- sort_hclust(mat_cluster_rows)

p <- pheatmap(jaccard_mat,
              cluster_cols = mat_cluster_cols,
              cluster_rows = mat_cluster_rows,
              color = viridis(n = 256))
ggsave(plot = p,
       filename = here(output_dir,
                       "JaccardIndex_Heatmap.pdf"),
       height = 6, width = 6,
       dpi = 600, useDingbats = FALSE)





# ==== Node Metric Calculations ====
# Calculate different node metrics for each species' metabolic reactions
# Node metric calculations functions in "MicrobialNetworkFunctions.R"
# Metrics calculated: Cascade number, degrees, clustering coefficient, eigenvector
#   centrality, and betweenness centrality


# ---- * Cascade number ----
# Read computed cascade numbers from script GetCascadeNumber
cascade_num <- read_csv(here(file_dir, "AllCascadeNum.csv"),
                        show_col_types = FALSE) %>%
  mutate(metric = "CascadeNum") %>%
  select(reaction, metric, Bsubtilis:Scerevisiae) %>%
  rename(rxn = reaction)


# ---- * Node degrees ----
node_degrees <- lapply(names(adj_mat), function(m){
  node_degrees <- get_all_degrees(adj_mat[[m]]) %>%
    mutate(species = m)
})
node_degrees <- combine_species_metrics(metric_list = node_degrees,
                                        metric_val_col = "k")


# ---- * Clustering coefficient ----
clustering_coefficients <- lapply(names(adj_mat), function(m){
  ccs <- get_clustering_coefficient(adj_mat[[m]]) %>%
    mutate(species = m)
})
clustering_coefficients <- combine_species_metrics(metric_list = clustering_coefficients,
                                                   metric_val_col = "clustering_coefficient")


# ---- * Eigenvector centrality ----
ev_centrality <- lapply(names(adj_mat), function(m){
  evs <- get_eigenvector_centrality(adj_mat[[m]]) %>%
    mutate(species = m)
})
ev_centrality <- combine_species_metrics(metric_list = ev_centrality,
                                         metric_val_col = "eigenvector_centrality")


# ---- * Betweenness centrality ----
btw_centrality <- lapply(names(adj_mat), function(m){
  btw <- get_betweenness_centrality(adj_mat[[m]]) %>%
    mutate(species = m)
})
btw_centrality <- combine_species_metrics(metric_list = btw_centrality,
                                          metric_val_col = "btw")


# ---- * Combine metrics ----
# Combine metrics together into a single table
node_metrics <- list()
node_metrics <- bind_rows(cascade_num,
                          node_degrees,
                          clustering_coefficients,
                          ev_centrality,
                          btw_centrality) %>%
  filter(metric !="total-degree") # total degree a bit repetitive with in and out degrees






# ==== Pearson Correlations ====
# Calculate the correlations of metric values between species
cor_res <- pblapply(unique(node_metrics$metric), function(met){
  
  # Subset for current metric
  df_met <- node_metrics %>%
    filter(metric == met)
  
  cor_res <- list() # store correlation output
  for (i in 1:nrow(species_comp)){
    s1 <- species_comp$X1[i]
    s2 <- species_comp$X2[i]
    df_met2 <- df_met %>%
      select(rxn, metric, all_of(c(s1, s2)))
    
    res <- cor.test(df_met2[[s1]],
                    df_met2[[s2]],
                    use = "pairwise.complete.obs",
                    method = "pearson")
    
    # Extract correlation value and confidence interval
    output <- data.frame(species1 = s1,
                         species2 = s2,
                         model = species_comp$model[i],
                         metric = met,
                         r = res$estimate,
                         r_low = res$conf.int[1],
                         r_high = res$conf.int[2],
                         pval = res$p.value)
    
    cor_res <- bind_rows(cor_res,
                         output)
    rownames(cor_res) <- NULL
  }
  return(cor_res)
})
cor_res <- do.call(bind_rows,
                   cor_res) %>%
  tibble()


# ---- * Variance vs jaccard index ----
# Get order of species comparison based on jaccard index
ord <- (species_comp %>%
  arrange(desc(jaccard)))$model

# Plot a heatmap of jaccard index
p1 <- ggplot(data = species_comp %>%
               mutate(dummy = 1) %>%
               mutate(model = factor(model, levels = ord)),
             mapping = aes(x = as.factor(dummy), y = model, fill = jaccard)) +
  geom_tile() +
  geom_text(mapping = aes(label = sprintf("%.2f", jaccard)),
            color = "white") +
  scale_fill_viridis(limits = c(0,1)) +
  xlab("Jaccard Index") +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

# Plot correlation results for each metric
p2 <- cor_res %>%
  left_join(species_comp, by = c("model")) %>%
  mutate(model = factor(model, levels = ord)) %>%
  ggplot(mapping = aes(x = r, y = model, 
                       fill = metric)) +
  geom_point(size = 7, alpha = 0.8, shape = 21, stroke = NA) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  xlab("Pearson Correlation") +
  scale_fill_manual(values = cmap_metrics)

# Combine and save plot
p <- egg::ggarrange(p1, p2, 
                    ncol = 2, widths = c(0.25, 2))
ggsave(filename = here(output_dir,
                       "MetricCorrelationVariance.pdf"),
       plot = p,
       height = 7, width = 9,
       dpi = 600, useDingbats = FALSE)


# ---- * 95% CI ----
# Plot the 95% CI of correlation values
ggdf <- cor_res %>%
  mutate(model = factor(model, levels = ord)) %>%
  pivot_longer(cols = all_of(c("r", "r_low", "r_high")),
               names_to = "r_vals",
               values_to = "values")

ggplot() +
  geom_point(data = ggdf %>%
               filter(r_vals == "r"),
             mapping = aes(x = values, y = model, fill = metric),
             shape = 21, stroke = NA, size = 3) +
  geom_line(data = ggdf,
            mapping = aes(x = values, y = model, color = metric)) +
  facet_wrap(~metric, ncol = 3, scales = "free_x") +
  xlab("Pearson Correlation with 95% intervals") +
  theme(axis.title.y = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = cmap_metrics) +
  scale_color_manual(values = cmap_metrics)
ggsave(filename = here(output_dir,
                       "MetricCorrelation95CI.pdf"),
       height = 6, width = 12, dpi = 600, useDingbats = FALSE)






# ==== Linear Model ====
# Look at the actual linear relationship of metrics between species
# For each metric, calculate the correlation of metric values between species

dir.create(here(output_dir, "MetricLinearRegression"), 
           showWarnings = FALSE)

lm_res <- invisible(pblapply(unique(node_metrics$metric), function(met){
  dir <- paste0(here(output_dir, "MetricLinearRegression/", met))
  dir.create(dir, showWarnings = FALSE)
  
  # Subset for current metric
  df_met <- node_metrics %>%
    filter(metric == met) 
  
  lm_res <- list() # store linear model output
  for (i in 1:nrow(species_comp)){
    
    s1 <- species_comp$X1[i]
    s2 <- species_comp$X2[i]
    design <- species_comp$model[i] # model design
    df_met2 <- df_met %>%
      select(rxn, metric, all_of(c(s1, s2)))
    
    # Linear regression
    lm_mod <- lm(as.formula(design), data = df_met)
    lm_sum <- summary(lm_mod)
    
    # Coefficient and 95% CI
    coeff <- coefficients(lm_sum)[2, 1]
    coeff_ci <- confint(lm_mod, level = 0.95)
    coeff_low <- coeff_ci[2,1]
    coeff_high <- coeff_ci[2,2]
    
    
    # calculate pvalue
    f <- lm_sum$fstatistic
    pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    
    # store results
    res <- data.frame(species1 = s1,
                      species2 = s2,
                      metric = met,
                      model = design,
                      r2 = lm_sum$adj.r.squared,
                      pval = pval,
                      coeff = coeff,
                      coeff_low = coeff_low,
                      coeff_high = coeff_high)
    lm_res <- bind_rows(lm_res,
                        res)
    rownames(lm_res) <- NULL
    
    # Plot correlations between species
    x_pos <- max(df_met[[s1]], na.rm = TRUE) - (max(df_met[[s1]], na.rm = TRUE) * 0.8)
    y_pos <- max(df_met[[s2]], na.rm = TRUE) * 0.8 # loc of R^2 label
    r_lab <- paste0("R^2 = ", sprintf("%.2f", lm_sum$adj.r.squared))
    
    ggplot(data = df_met,
           mapping = aes(x = !!as.symbol(s1),
                         y = !!as.symbol(s2))) +
      geom_smooth(method = "lm",
                  se = FALSE,
                  formula = "y~x") +
      geom_point(alpha = 0.8) +
      geom_text(x = x_pos, y = y_pos,
                label = r_lab) +
      ggtitle(met)
    
    ggsave(filename = paste0(dir, "/",
                             met, "_", s1, "-", s2, ".pdf"),
           height = 3, width = 3, dpi = 600, useDingbats = FALSE)
  }
  return(lm_res)
}))
# Combine linear model results into a table
lm_res <- do.call(bind_rows,
                  lm_res) %>%
  tibble()






# ---- * 95% CI ----
ggdf <- lm_res %>%
  mutate(model = factor(model, levels = ord)) %>%
  pivot_longer(cols = all_of(c("coeff", "coeff_low", "coeff_high")),
               names_to = "coeff_vals",
               values_to = "values")

ggplot() +
  geom_vline(xintercept = 1,
             linetype = "dashed") +
  geom_point(data = ggdf %>%
               filter(coeff_vals == "coeff"),
             mapping = aes(x = values, y = model, 
                           color = metric)) +
  geom_line(data = ggdf,
            mapping = aes(x = values, y = model, 
                          color = metric)) +
  facet_wrap(~metric, ncol = 2, scales = "free")
ggsave(filename = here(output_dir,
                       "Metric95CI.pdf"),
       height = 9, width = 16)
