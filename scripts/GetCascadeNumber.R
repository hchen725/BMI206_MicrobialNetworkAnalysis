# ==== Setup ====
# ---- * Load packages ----
library(here)
library(matlabr)
library(tidyverse)




# ==== Functions ====
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





# ==== Split Adjacency Matrix ====
# Load adjacency matrix downloaded from supplemental table 2 of Kim et al.
# Split excel sheets into csv files that can be used for cascade number
# calculations performed by CascadeNum.M

adj_mat_path <- here("data/AdjacencyMatrix.xlsx")
dir.create(here("data/AdjacencyMatrixCSV"))
dir.create(here("data/RxnListCSV"))

adj_mat <- read_excel_allsheets(filename = adj_mat_path)

for (m in names(adj_mat)){

  mat <- adj_mat[[m]] %>%
    column_to_rownames(var = "...1") 
  rxn_names <- data.frame(rxn = colnames(mat))

  # Write CSV
  write_csv(mat,
            file = here("data/AdjacencyMatrixCSV",
                        paste0("AdjMat_", m, ".csv")),
            col_names = FALSE)
  write_csv(rxn_names,
            file = here("data/RxnListCSV",
                        paste0("RxnList_", m, ".csv")))
}





# ==== RunCascade.M ====
# Run the MATLAB script to calculate the cascade number for each species
cascade_script <- here("scripts/RunCascade.m")
run_matlab_script(fname = cascade_script)

# ==== Combine with reaction list ====
# Combine calculated cascade number with the reaction list
cascade_path <- here("data/CascadeNumbers")
rxn_path <- here("data/RxnListCSV")
ls <- data.frame(species = c("Bsubtilis", 
                             "Ecoli", 
                             "Gmetallireducens", 
                             "Kpneumoniae", 
                             "Scerevisiae"),
                 cascade = list.files(cascade_path),
                 rxn = list.files(rxn_path))

cascade_num <- lapply(1:nrow(ls), function(i){
  s <- ls$species[i]
  rxn <- read_csv(here(rxn_path, ls$rxn[i]),
                  show_col_types = FALSE)
  rxn <- rxn$rxn
  cs <- read_csv(here(cascade_path, ls$cascade[i]),
                 col_names = FALSE,
                 show_col_types = FALSE)
  
  colnames(cs) <- rxn
  cs <- data.frame(t(cs)) %>%
    rownames_to_column("reaction")
  colnames(cs) <- c("reaction", s)
  return(cs)
})
names(cascade_num) <- ls$species

# Combine all cascade numbers together
all_cascade_num <- cascade_num[[1]]
for (i in 2:length(cascade_num)){
  all_cascade_num <- full_join(all_cascade_num, cascade_num[[i]],
                                by = c("reaction"))
}

write_csv(all_cascade_num,
          file = here("data/AllCascadeNum.csv"))
