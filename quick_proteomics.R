# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

# location that the all peptides csv is in
setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# file name of the all peptides csv
data_raw <- read.csv("Kowalski_F_w1_w8_all_peptides.csv")


proteome <- function(dat, samp_id = "F", max_pep = 3){
  
  # data frame check
  if (colnames(dat)[1:3] != c("Annotated.Sequence", "Modifications", 
                              "Master.Protein.Accessions"))
    stop("invalid data frame: must have Annotated.Sequence, Modifications, 
         Master.Protein.Accessions")
  
  
  wt_grp <- grep(samp_id, colnames(dat))
  data_raw$wt_grp_med <- by_group(dat, wt_grp)
  
  data <- tbl_df(data_raw) %>%
    group_by(Master.Protein.Accessions) %>%
    top_n(n = max_pep, wt = wt_grp_med)
  
  
  
  
}


# testing out conditional error throwing
test <- function(norm) {
  
  if (norm != FALSE & norm != TRUE) 
    stop ("invalid input for norm argument: must contain TRUE or FALSE")
  
print("hello")
  
}

test(FALSE)

test("df")
