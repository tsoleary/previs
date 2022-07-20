# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Beach/Bulk Abundances")
data_raw <- read.csv("D3 Cells mouse DB All Peptides.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("0Hr", colnames(data_raw))

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

data_raw <- filter(data_raw, (is.na(ctrl_raw_med) == FALSE))

# set the max number of peptides used in analysis
max_pep <- 3

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)


# proteins used for normalization
GAPDH <- "A0A0A0MQF6"

norm_pro <- GAPDH

norm_pep <- subset(data, data$Master.Protein.Accessions == norm_pro)
numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
raw_abun <- numeric_cols[-length(numeric_cols)]
norm_value <- sapply(norm_pep[, raw_abun], mean, na.rm = TRUE)

raw_abun_mat <- as.matrix(data[, raw_abun])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data <- as_tibble(data)
data <- cbind(data, norm_test)


# # sum of every column
# norm_value <- colSums(data_raw[, ctrl_raw], na.rm = TRUE)
# 
# raw_abun_mat <- as.matrix(data[, ctrl_raw])
# 
# norm_abun <- t(t(raw_abun_mat)/norm_value)
# colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
# norm_test <- as.data.frame(norm_abun)
# data <- as_tibble(data)
# data <- cbind(data, norm_test)

# # TA muscle mass (mg) and sum-total (one operation)
# 
# sample_id_mass <- read.csv("PSD proteomic sample info M w1_8.csv")
# ta_mass <- sample_id_mass$TA.mass
# names(ta_mass) <- sample_id_mass$ID
# 
# 
# norm_value <- colSums(data_raw[, ctrl_raw], na.rm = TRUE)
# 
# raw_abun_mat <- as.matrix(data[, ctrl_raw])
# 
# 
# norm_abun <- t(t(raw_abun_mat)*ta_mass/norm_value)
# colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
# norm_test <- as.data.frame(norm_abun)
# data <- as_tibble(data)
# data <- cbind(data, norm_test)

# Data frame with only top few ionizing peptides -------------------------------

# pep_top <- 3
# data_top <-
#   tbl_df(data) %>%
#   group_by(Master.Protein.Accessions) %>%
#   top_n(n = pep_top, wt = ctrl_raw_med) %>%
#   as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- colnames(data)[grep("norm", colnames(data))]

protein <- by_protein(data, group_names) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv('D3 Cells mouse DB genelist.csv')

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 2 
protein$peptides <- table(data_raw$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "Top 3 GAPDH Norm.csv")
write.csv(data, "Top 3 GAPDH Norm Peptides.csv")
