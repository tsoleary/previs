# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

# location that the all_peptides.csv is in
setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# file name of the all_peptides.csv
data_raw <- read.csv("Kowalski_F_w1_w8_all_peptides.csv")


proteome <- function(dat, group_names, file_gene, samp_id = "F", top_pep = 3, min_pep = 3,
                     norm = FALSE, 
                     norm_method = "sum_total", norm_pro = NULL){
  
  # check the arguments to make sure they are all in the correct format
  
  # data frame check
  if (colnames(dat)[1] != "Annotated.Sequence" &
      colnames(dat)[2] != "Modifications" &
      colnames(dat)[3] != "Master.Protein.Accessions")
    stop("invalid data frame: must have Annotated.Sequence, Modifications, 
         Master.Protein.Accessions")
  
  #samp_id check
  
  # top_pep
  
  
  
  # determine the top ionizers and sort by them taking only the top_pep
  wt_grp <- grep(samp_id, colnames(dat))
  dat$wt_grp_med <- by_group(dat, wt_grp)
  
  df <- tbl_df(dat) %>%
    group_by(Master.Protein.Accessions) %>%
    top_n(n = top_pep, wt = wt_grp_med)
  
  
  # do normalization if indicated
  if (norm == FALSE){
    
    if (norm_method == "sum_total") {
    norm_value <- colSums(df[, wt_grp], na.rm = TRUE)
    
    raw_abun_mat <- as.matrix(df[, wt_grp])
    
    norm_abun <- t(t(raw_abun_mat)/norm_value)
    colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
    norm_test <- as.data.frame(norm_abun)
    df <- as_tibble(df)
    df <- cbind(df, norm_test)
    
    }
    
    if (norm_method == "protein"){
      
    norm_pep <- subset(df, df$Master.Protein.Accessions == norm_pro)
    numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
    raw_abun <- numeric_cols[-length(numeric_cols)]
    norm_value <- sapply(norm_pep[, raw_abun], mean, na.rm = TRUE)
    raw_abun_mat <- as.matrix(df[, raw_abun])
    norm_abun <- t(t(raw_abun_mat)/norm_value)
    colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
    norm_test <- as.data.frame(norm_abun)
    df <- as_tibble(df)
    df <- cbind(df, norm_test)
    
    
    }
  }

  
  pro_df <- by_protein(df, colnames(df)[grep(group_names, colnames(df))]) %>%
    as.data.frame %>%
    rownames_to_column("Master.Protein.Accessions")

  pro_df$gene <- mpa_to_gene(pro_df, read.csv(file_gene))
  
  pro_df$peptides <- table(dat$Master.Protein.Accessions)
  pro_df <- filter(pro_df, pro_df$peptides >= min_pep)
  
  write.csv(pro_df, 
            paste0("r_output_", format(Sys.time(), "%d %b %Y"), ".csv"))
  
  return(pro_df)
}



x <- proteome(data_raw, "norm", "Kowalski_F_w1_w8_gene_names.csv")

