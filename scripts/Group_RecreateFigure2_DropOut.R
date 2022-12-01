# Recreate Figure 2 with drop outs

# ==== Setup ====
library(here)
library(tidyverse)
library(pbapply)
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

# ==== Observed ====
# Node degrees
k_all <- lapply(names(adj_mat), function(m){
  k_all <- get_all_degrees(adj_mat[[m]])
})
names(k_all) <- names(adj_mat)

# Coefficients
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


# ==== Iterations ====

# Save output directory
output_dir <- here("output/han/Fig2_DropOut")
dir.create(output_dir)

# Filename to save the iterations as
filename <- "iterations.rds"
if(file.exists(paste0(output_dir, "/", filename))){
  
  cat("Reading previously computed iterations.",
      sep = "\n")
  iteration_output <- read_rds(file = paste0(output_dir, "/", filename))
  
} else {
  
  cat("Performing drop out iterations.",
      sep = "\n")
  
  # Sample adjacency matrix
  sample_fraction <- c(0.1, 0.25, 0.5, 0.8, 0.9, 0.95)
  num_iterations <- 1000
  
  iteration_output <- lapply(names(adj_mat), function(species){
    
    # Get adjacency matrix for the species
    mat <- adj_mat[[species]]
    rxns <- colnames(mat)
    
    # Sample at different fractions
    sampled_mat <- lapply(sample_fraction, function(frac){
      
      cat(paste0("Current species: ", species),
          paste0("Performing ", num_iterations, " iterations with sample fraction: ", frac),
          sep = "\n")
      
      # Perform iterations
      coeff_its <- pblapply(1:num_iterations, function(i){
        
        rxns_sample <- sample(rxns,
                              size = length(rxns) * frac,
                              replace = FALSE)
        mat_sample <- mat[rxns_sample, rxns_sample]
        
        # Calculate node degrees
        k <- get_all_degrees(mat_sample)
        
        # Calculate coefficients
        glm_coeff <- lapply(c("in-degree", "out-degree", "total-degree"), function(deg){
          k_sub <- k %>%
            filter(metric == deg,
                   k != 0)
          # Assume a gaussian glm on log10 values of prk and k
          glm_mod <- glm(log10(prk)~log10(k),
                         data = k_sub,
                         family = gaussian)
          glm_sum <- summary(glm_mod)
          glm_coef <- data.frame(species = species,
                                 metric = deg,
                                 coeff = coef(glm_sum)[2,1])
          return(glm_coef)
        })
        
        glm_coeff <- do.call(bind_rows,
                             glm_coeff) %>%
          mutate(iteration = i,
                 sample_fraction = frac)
        
      })
      coeff_its <- do.call(bind_rows,
                           coeff_its)
      return(coeff_its)
      
    })
    names(sampled_mat) <- sample_fraction
    return(sampled_mat)
  })
  names(iteration_output) <- names(adj_mat)
  # Save output
  write_rds(iteration_output,
            file = paste0(output_dir, "/", filename))
}


# ==== Plot ====
invisible(lapply(names(adj_mat), function(s){
  ord <- rev(sample_fraction)
  species_iteration <- do.call(bind_rows,
                               iteration_output[[s]]) %>%
    mutate(sample_fraction = factor(sample_fraction,
                                    levels = ord))
  coef_lines <- glm_coef_all %>%
    dplyr::filter(species == s)
    
  ggplot() +
    # Plot iterations
    geom_density(data = species_iteration,
                 mapping = aes(x = coeff, 
                               color = metric, 
                               fill = metric),
                 alpha = 0.5) +
    # Plot observed value
    geom_vline(data = coef_lines,
               mapping = aes(xintercept = coeff,
                             color = metric),
               linetype = "dashed") +
    facet_wrap(metric~sample_fraction, nrow = 3, scales = "free_y") +
    ggtitle(s) +
    scale_color_manual(values = c("in-degree" = "red",
                                  "out-degree" = "blue",
                                  "total-degree" = "black")) +
    scale_fill_manual(values = c("in-degree" = "red",
                                 "out-degree" = "blue",
                                 "total-degree" = "black")) 
  ggsave(filename = paste0(output_dir, "/", s, ".png"),
         height = 7, width = 15)
}))



