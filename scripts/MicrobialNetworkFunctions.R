# Microbial Network Functions
# ==== General ====

#' read_excel_allsheets
#'
#' @param filename .XLSX file to read
#' @param tibble T/F, store list elements as a tibble or dataframe
#'
#' @return List containing all sheets present in excel files
#' @description Reads all the sheets present in an excel file
read_excel_allsheets <- function(filename,
                                 tibble = TRUE){
  # Get all the sheets present in filename
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets,
              function(x){
                readxl::read_excel(filename, sheet = x)
              })
  if(!tibble){
    x <- lapply(x, as.data.frame)
  }
  names(x) <- sheets
  return(x)
}


# ==== Node Metrics ====

# Format all metric calculation the same where
# rows = reactions
# Minimum columns to have
# columns = reaction name, metric column, metric name, species


# https://rpubs.com/pjmurphy/313180
# local metrics in igraph:
# degree centrality
# eigenvector centrality
# hubs and authorities
# closeness centrality
# reach centrality
# betweenness centrality
combine_species_metrics <- function(metric_list,
                                    reaction_col = "rxn",
                                    metric_col = "metric",
                                    species_col = "species",
                                    metric_val_col = "value"){
  
  metric_list <- do.call(bind_rows,
                         metric_list) %>%
    select(all_of(c(reaction_col, metric_col, species_col, metric_val_col))) %>%
    pivot_wider(names_from = all_of(species_col),
                values_from = all_of(metric_val_col))
  
  return(metric_list)
}



# ---- * Degrees (k) and distributions (Pr(k)) ----
# Calculate out-degrees, in-degrees, and total-degrees from an adjacency matrix
get_degree_out <- function(mat){
  rxn_names <- rownames(mat)
  total_nodes <- nrow(mat)
  k_out <- mat %>%
    colSums()
  prk_out <- table(k_out) %>%
    data.frame() %>%
    mutate(prk_out = Freq/total_nodes,
           k_out = as.numeric(as.character(k_out)))
  k_out <- data.frame(k_out = k_out) %>%
    left_join(prk_out,
              by = c("k_out"))
  rownames(k_out) <- rxn_names
  return(k_out)
}

get_degree_in <- function(mat){
  rxn_names <- rownames(mat)
  total_nodes <- nrow(mat)
  k_in <- mat %>%
    rowSums()
  prk_in <- table(k_in) %>%
    data.frame() %>%
    mutate(prk_in = Freq/total_nodes,
           k_in = as.numeric(as.character(k_in)))
  k_in <- data.frame(k_in = k_in) %>%
    left_join(prk_in,
              by = c("k_in"))
  rownames(k_in) <- rxn_names
  return(k_in)
}

get_degree_total <- function(k_in, k_out){
  rxn_names <- k_in$rxn
  total_nodes <- nrow(k_in)
  k_in <- k_in %>%
    select(k)
  k_out <- k_out %>%
    select(k)
  k_total <- k_in + k_out
  colnames(k_total) <- c("k_total")
  
  prk_total <- table(k_total) %>%
    data.frame() %>%
    mutate(prk_total = Freq/total_nodes,
           k_total = as.numeric(as.character(k_total)))
  k_total <- data.frame(k_total = k_total) %>%
    left_join(prk_total,
              by = c("k_total"))
  rownames(k_total) <- rxn_names
  return(k_total)
}

# Calculate all the degree metrics for an adjacency matrix in a single function
get_all_degrees <- function(mat){
  k_in <- get_degree_in(mat) %>%
    rownames_to_column("rxn") %>%
    rename(k = k_in,
           prk = prk_in) %>%
    mutate(metric = "in-degree")
  k_out <- get_degree_out(mat) %>%
    rownames_to_column("rxn") %>%
    rename(k = k_out,
           prk = prk_out) %>%
    mutate(metric = "out-degree")
  k_total <- get_degree_total(k_in, k_out) %>%
    rownames_to_column("rxn") %>%
    rename(k = k_total,
           prk = prk_total) %>%
    mutate(metric = "total-degree")
  k_all <- bind_rows(k_in,
                     k_out,
                     k_total)
  return(k_all)
}

# Format the degree coefficients for labeling of plots (used in recreation of Figure 2)
format_degree_coeffs <- function(coeffs){
  
  in_deg <- sprintf("%.2f", 
                    coeffs[which(coeffs$metric == "in-degree"), 3])
  out_deg <- sprintf("%.2f", 
                     coeffs[which(coeffs$metric == "out-degree"), 3])
  total_deg <- sprintf("%.2f", 
                       coeffs[which(coeffs$metric == "total-degree"), 3])
  
  lab <- list(total_deg = paste0("\u03b3 (total) = ", total_deg),
              in_deg = paste0( "\u03b3 (in) = ", in_deg),
              out_deg = paste0("\u03b3 (out) = ", out_deg))
  
  return(lab)
}


# ---- * Clustering Coefficient ----
get_clustering_coefficient <- function(mat){
  # Get graph from adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(mat), 
                                   mode = "directed")
  
  # Calculate the local clustering coefficient 
  cc <- transitivity(graph = g,
                     type = "local") 
  
  cc <- tibble(rxn = colnames(mat),
               clustering_coefficient = cc,
               metric = "clustering_coefficient")
  
  return(cc)
}



# ==== Centralization ====
get_eigenvector_centrality <- function(mat){
  # Get graph from adjacency matrix
  g <- graph_from_adjacency_matrix(as.matrix(mat),
                                   mode = "directed") 
  
  # Calculating eigenvector centrality
  eig <- evcent(g)$vector
  
  eig <- tibble(rxn = colnames(mat),
                eigenvector_centrality = eig,
                metric = "eigenvector_centrality")
  
  return(eig)
}

get_betweenness_centrality <- function(mat){
  g <- graph_from_adjacency_matrix(as.matrix(mat),
                                   mode = "directed")
  
  btw <- betweenness(g)
  btw <- tibble(rxn = colnames(mat),
                btw = btw,
                metric = "btw_centrality")
  return(btw)
}
