# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------
# Markus Meyer HFpEF samples

library(tidyverse)

# source the functions in the proteomics_functions.R script!

setwd("C:/Users/PrevBeast/Documents/R/Meyer")
data_raw <- read.csv("meyer_all_peptides.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("Control", colnames(data_raw))

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

# set the max number of peptides used in analysis
max_pep <- 15

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# proteins used for normalization
myosin <- "P12883; P13533"
histones <- "B4DR52; P06899"

norm_pro <- histones

norm_pep <- subset(data, data$Master.Protein.Accessions == norm_pro)
numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
raw_abun <- numeric_cols[-length(numeric_cols)]
norm_value <- sapply(norm_pep[, raw_abun], mean)

raw_abun_mat <- as.matrix(data[, raw_abun])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data <- as_tibble(data)
data <- cbind(data, norm_test)

# Median, sd, & ratio of peptides ----------------------------------------------

group1 <- grep("HFpEF_norm", colnames(data))
ctrl <- grep("Control_norm", colnames(data))

data$group1_med <- by_group(data, group1)
data$ctrl_med <- by_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- by_group(data, group1, FUN = sd)
data$ctrl_sd <- by_group(data, ctrl, FUN = sd)

data$ratio <- abun_ratio(data, "group1_med")

data$group1_df <- count_df(data, group1)
data$ctrl_df <- count_df(data, ctrl)

# Removing rows with NA values for the median
data <- data[!(is.na(data$ratio)), ]

# Remove outliers --------------------------------------------------------------

pro_out <- by_protein(data, "ratio") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

data <- rm_outliers(data, pro_out, "ratio")

# Oxidized peptides Only -------------------------------------------------------

# data_ox <- filter(data, grepl("Oxid", Modifications))
# data <- data_ox

# Collagen peptides without oxidation ------------------------------------------
# converting protein accession to gene symbol
# gene_df <- read.csv('meyer_protein_accession_gene_list.csv')
# 
# data$gene <- mpa_to_gene(data, gene_df)
# 
# data_col <- filter(data, grepl("COL", gene))
# data <- filter(data_col, !grepl("Oxid", Modifications))

# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- c("group1_med", "ctrl_med")

protein <- by_protein(data_top, group_names) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

protein$Master.Protein.Accessions <-
  protein$Master.Protein.Accessions %>%
  as.character

data_top$group1_sd_df <- square_x_df(data_top, "group1_sd", "group1_df")
data_top$ctrl_sd_df <- square_x_df(data_top, "ctrl_sd", "ctrl_df")

# Creating a temp data frame to calculate the ratio sd
temp_df <- NULL 
sd_sum <- c("group1_sd_df", "ctrl_sd_df", "group1_df", "ctrl_df")
temp_df <- by_protein(data_top, sd_sum, FUN = sum) %>% as.data.frame

group1_pooled_sd <- temp_df$group1_sd_df / temp_df$group1_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
protein <- cbind(protein, group1_pooled_sd, ctrl_pooled_sd)

# Grouped relative protein abundance ratio
protein$ratio <- by_protein(data, "ratio") %>% as.numeric

# Standard deviation relative protein abundance ratio
data$ratio_sd <- ratio_sd(data, "ratio", "group1_sd", "group1_med", 
                             "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio
data$ratio_df <- data$group1_df + data$ctrl_df
data$ratio_sd_df <- square_x_df(data, "ratio_sd", "ratio_df")

sd_calc <- c("ratio_df", "ratio_sd_df")
temp_df <- by_protein(data, sd_calc, FUN = sum) %>% as.data.frame

ratio_sd <- temp_df$ratio_sd_df / temp_df$ratio_df
protein <- cbind(protein, ratio_sd)

# Statistics -------------------------------------------------------------------

log_norm <- log(data[, c(group1, ctrl)])
colnames(log_norm) <- paste(colnames(log_norm), sep = "_", "log")
data <- cbind(data, log_norm)

log_cols <- grep("log", colnames(data))
group1_log_cols <- grep("HFpEF_norm_log", colnames(data))
ctrl_log_cols <- grep("Control_norm_log", colnames(data))

# Problem if different number of samples be in each group ----------------------

length(group1_log_cols)
length(ctrl_log_cols)

#add two cols to HFpEF samples
data$dummy1_HFpEF_norm_log <- NA
data$dummy2_HFpEF_norm_log <- NA

group1_log_cols <- grep("HFpEF_norm_log", colnames(data))
ctrl_log_cols <- grep("Control_norm_log", colnames(data))

stacked <- data.frame(data[, "Master.Protein.Accessions"], 
                           stack(data[, group1_log_cols]),
                           stack(data[, ctrl_log_cols]))

colnames(stacked) <- c("Master.Protein.Accessions", "group1_log", "samp", 
                        "ctrl_log", "ctrl")

p_vals <- pval_ttest(stacked, "group1_log", "ctrl_log")
colnames(p_vals)[1] <- "Master.Protein.Accessions"

p_vals$Master.Protein.Accessions <-
  p_vals$Master.Protein.Accessions %>% 
    as.character

protein <- full_join(protein, p_vals, by = "Master.Protein.Accessions")

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv('meyer_protein_accession_gene_list.csv')

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 5 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "meyer_HFpEF_hist_norm_col_non_ox_pep.csv")