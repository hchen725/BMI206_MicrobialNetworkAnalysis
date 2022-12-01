# Recreate Figure 2

# ==== Setup ====
library(here)
library(tidyverse)
library(scales)

theme_set(theme_bw())

# Get the node degree calculation functions
source(here("scripts/MicrobialNetworkFunctions.R"))

# ---- * Load data ----
file_dir <- here("data/")
adj_mat <- read_excel_allsheets(filename = here(file_dir, 
                                                "AdjacencyMatrix.xlsx"),
                                tibble = FALSE)

adj_mat_names <- names(adj_mat)
# Make all adj_mat into matrices
adj_mat <- lapply(names(adj_mat), function(x){
  m <- adj_mat[[x]] %>%
    column_to_rownames("...1")
})
names(adj_mat) <- adj_mat_names

# ==== Calculate node degrees ====
k_all <- lapply(names(adj_mat), function(m){
  k_all <- get_all_degrees(adj_mat[[m]])
})
names(k_all) <- names(adj_mat)

# Degree distributions 
output_dir <- here("output/k_distributions")
dir.create(output_dir)
invisible(lapply(names(k_all), function(species){
  
  k_sub <- k_all[[species]]
  
  ggplot(data = k_sub,
         mapping = aes(x = k, color = metric, fill = metric)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = c("in-degree" = "red",
                                  "out-degree" = "blue",
                                  "total-degree" = "black")) +
    scale_fill_manual(values = c("in-degree" = "red",
                                 "out-degree" = "blue",
                                 "total-degree" = "black")) +
    ggtitle(species) +
    xlab("Degree (k)")
  ggsave(filename = paste0(output_dir, "/Raw_", species, ".png"),
         height = 4, width = 4.6)
  
  ggplot(data = k_sub,
         mapping = aes(x = k, color = metric, fill = metric)) +
    geom_density(alpha = 0.5) +
    scale_x_log10(limits = c(1, 200),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  name = "Degree (k)") +
    scale_color_manual(values = c("in-degree" = "red",
                                  "out-degree" = "blue",
                                  "total-degree" = "black")) +
    scale_fill_manual(values = c("in-degree" = "red",
                                 "out-degree" = "blue",
                                 "total-degree" = "black")) +
    ggtitle(species)
  
  ggsave(filename = paste0(output_dir, "/LogTrans_", species, ".png"),
         height = 4, width = 4.6)
}))


# ==== Slopes ====
# Calculating slope for each degree metric

glm_coef_all <- lapply(names(k_all), function(species){
  glm_coef <- lapply(c("in-degree", "out-degree", "total-degree"), function(deg){
    k_sub <- k_all[[species]] %>%
      filter(metric == deg,
             k != 0)
    # Assume a gaussian glm on log10 values of prk and k
    glm_mod <- glm(log10(prk)~log10(k),
                   data = k_sub,
                   family = gaussian)
    glm_sum <- summary(glm_mod)
    #b0 <- coef(glm_sum)[1,1]
    b1 <- coef(glm_sum)[2,1]
    glm_coef <- data.frame(species = species,
                           metric = deg,
                           coeff = b1)
    return(glm_coef)
  })
  glm_coef <- do.call(bind_rows,
                      glm_coef)
})
glm_coef_all <- do.call(bind_rows,
                        glm_coef_all)



# ==== Figure 2 ====

# Degree distribution in the reaction-centric metabolic network
# X-axis is the degree (k) and the y-axis is the probability that a metabolite
#   has degree k (Pr(K))
# Degree of a node (k) - sum of all links connected to node
# Degree distribution (Pr_k) - probability that a randomly picked node in a network
#   has a degree of k
output_dir <- here("output/Figure2")
dir.create(output_dir)



invisible(lapply(names(adj_mat), function(m){
  coeffs <- glm_coef_all %>%
    filter(species == m)
  
  lab <- format_degree_coeffs(coeffs)
  k_all[[m]] %>%
    filter(k != 0) %>%
    ggplot(mapping = aes(x = k, y = prk, 
                         color = metric,
                         shape = metric,
                         group = metric)) +
    stat_smooth(method = "lm",
                se = FALSE,
                fullrange = TRUE) +
    geom_point() +
    # Add gamma coefficients
    annotate("text",
             x = 10^1.25, y = 10^-0.3,
             label = lab$total_deg,
             color = "black",
             hjust = 0,
             size = 3) +
    annotate("text",
             x = 10^1.25, y = 10^-0.5,
             label = lab$in_deg,
             color = "red",
             hjust = 0,
             size = 3) +
    annotate("text",
             x = 10^1.25, y = 10^-0.7,
             label = lab$out_deg,
             color = "blue",
             hjust = 0,
             size = 3) +
    ggtitle(m) +
    scale_x_log10(limits = c(1, 200),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  name = "Degree (k)") +
    scale_y_log10(limits = c(0.0001, 1),
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  name = "Pr(k)") +
    scale_color_manual(values = c("in-degree" = "red",
                                  "out-degree" = "blue",
                                  "total-degree" = "black")) +
    scale_shape_manual(values = c("in-degree" = 15, # square
                                  "out-degree" = 17, # triangle
                                  "total-degree" = 16)) # circle 
  
  ggsave(filename = paste0(output_dir, "/", m, ".png"),
         height = 4, width = 4.6)
}))


# Ecoli example without log transformation
k_all[["Ecoli"]] %>%
  filter(k != 0) %>%
  ggplot(mapping = aes(x = k, y = prk, 
                       color = metric,
                       shape = metric,
                       group = metric)) +
  geom_point() +
  ggtitle("Ecoli") +
  scale_color_manual(values = c("in-degree" = "red",
                                "out-degree" = "blue",
                                "total-degree" = "black")) +
  scale_shape_manual(values = c("in-degree" = 15, # square
                                "out-degree" = 17, # triangle
                                "total-degree" = 16)) + # circle 
  xlab("Degree (k)") + ylab("Pr(k)")
ggsave(filename = paste0(output_dir, "/Raw_Ecoli.png"),
       height = 4, width = 4.6)
