# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

# location that the all_peptides.csv is in
setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# check that the all_peptides.csv is in the working directory
list.files()

proteome("WT vs KO all pep.csv", group_names = c("WT", "KO"), 
         organism = "mouse", group = TRUE, norm = TRUE, 
         norm_method = "protein", norm_pro = "B2RTM0", csv = FALSE)








