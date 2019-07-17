# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
library(proteomixr)

# location that the all_peptides.csv is in
setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")

# check that the all_peptides.csv is in the working directory
list.files()

z <- proteome("WT vs KO all pep.csv", organism = "mouse", 
              wt_samps = "all", group = TRUE, group_names = c("WT", "KO"), 
              norm = TRUE, norm_method = "protein", 
              norm_pro = "B2RQQ1; Q91Z83", csv = FALSE)


# location that the all_peptides.csv is in
setwd("C:/Users/PrevBeast/Documents/R/HCM")

# check that the all_peptides.csv is in the working directory
list.files()

z <- proteome("HCM_duplicate_2_all_peptides.csv", organism = "human", 
              wt_samps = "all", group = TRUE, group_names = c("GOL", "DIS_1", 
              "DIS_2", "DIS_3"), norm = TRUE, norm_method = "sum_total")







