# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

# location that the all_peptides.csv is in
setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")

# check that the all_peptides.csv is in the working directory
list.files()

z <- proteome("WT vs KO all pep.csv", organism = "mouse", 
              wt_samps = "all", group = TRUE, group_names = c("WT", "KO"), 
              norm = TRUE, norm_method = "protein", 
              norm_pro = "B2RQQ1; Q91Z83", csv = FALSE)










