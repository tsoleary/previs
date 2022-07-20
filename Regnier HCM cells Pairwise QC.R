# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------
# SNORD KO vs WT at 6wk and 12wk

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Spees/Titin Splicing")

data_raw <- read.csv("6wk and 12wk Ventricles Male only all peps.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("WT_6wk", colnames(data_raw))
samps <- grep("F", colnames(data_raw))

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

colnames(data_raw)[1:3] <-  c("Annotated.Sequence", "Modifications", 
                              "Master.Protein.Accessions")

# set the max number of peptides used in analysis
max_pep <- 15

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# proteins used for normalization
# titin <- "D3DPG0"
# alphabeta_shared <- "P12883; P13533"
# 
# norm_pro <- alphabeta_shared
# 
# norm_pep <- subset(data, data$Master.Protein.Accessions == norm_pro)
# numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
# raw_abun <- numeric_cols[-length(numeric_cols)]
# norm_value <- sapply(norm_pep[, raw_abun], mean, na.rm = TRUE)
# 
# raw_abun_mat <- as.matrix(data[, raw_abun])
# 
# norm_abun <- t(t(raw_abun_mat)/norm_value)
# colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
# norm_test <- as.data.frame(norm_abun)
# data <- as_tibble(data)
# data <- cbind(data, norm_test)

# sum of every column
norm_value <- colSums(data_raw[, samps], na.rm = TRUE)

raw_abun_mat <- as.matrix(data[, samps])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data <- as_tibble(data)
data <- cbind(data, norm_test)

# Median, sd, & ratio of peptides ----------------------------------------------

group1 <- grep("WT_12wk_norm", colnames(data))
group2 <- grep("KO_6wk_norm", colnames(data))
group3 <- grep("KO_12wk_norm", colnames(data))
ctrl <- grep("WT_6wk_norm", colnames(data))

data$group1_med <- by_group(data, group1)
data$group2_med <- by_group(data, group2)
data$group3_med <- by_group(data, group3)
data$ctrl_med <- by_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- by_group(data, group1, FUN = sd)
data$group2_sd <- by_group(data, group2, FUN = sd)
data$group3_sd <- by_group(data, group3, FUN = sd)
data$ctrl_sd <- by_group(data, ctrl, FUN = sd)

# ratios
data$ratio1 <- abun_ratio(data, "group1_med")
data$ratio2 <- abun_ratio(data, "group2_med")
data$ratio3 <- abun_ratio(data, "group3_med")

# Degrees of freedom
data$group1_df <- count_df(data, group1)
data$group2_df <- count_df(data, group2)
data$group3_df <- count_df(data, group3)
data$ctrl_df <- count_df(data, ctrl)

# Removing rows with NA values for the median
data <- data[!(is.na(data$ratio1)), ]
data <- data[!(is.na(data$ratio2)), ]
data <- data[!(is.na(data$ratio3)), ]

# Remove outliers --------------------------------------------------------------

pro_out1 <- by_protein(data, "ratio1") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

data <- rm_outliers(data, pro_out1, "ratio1")

pro_out2 <- by_protein(data, "ratio2") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

data <- rm_outliers(data, pro_out2, "ratio2")

pro_out3 <- by_protein(data, "ratio3") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

data <- rm_outliers(data, pro_out3, "ratio3")

# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- c("group1_med", "group2_med", "group3_med", "ctrl_med")

protein <- by_protein(data_top, group_names) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

protein$Master.Protein.Accessions <-
  protein$Master.Protein.Accessions %>%
  as.character

data_top$group1_sd_df <- square_x_df(data_top, "group1_sd", "group1_df")
data_top$group2_sd_df <- square_x_df(data_top, "group2_sd", "group2_df")
data_top$group3_sd_df <- square_x_df(data_top, "group3_sd", "group3_df")
data_top$ctrl_sd_df <- square_x_df(data_top, "ctrl_sd", "ctrl_df")

# Creating a temp data frame to calculate the ratio1 sd
temp_df <- NULL 
sd_sum <- c("group1_sd_df", "ctrl_sd_df", "group1_df", "ctrl_df")
temp_df <- by_protein(data_top, sd_sum, FUN = sum) %>% as.data.frame

group1_pooled_sd <- temp_df$group1_sd_df / temp_df$group1_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
protein <- cbind(protein, group1_pooled_sd)

# Creating a temp data frame to calculate the ratio2 sd
temp_df <- NULL 
sd_sum <- c("group2_sd_df", "ctrl_sd_df", "group2_df", "ctrl_df")
temp_df <- by_protein(data_top, sd_sum, FUN = sum) %>% as.data.frame

group2_pooled_sd <- temp_df$group2_sd_df / temp_df$group2_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
protein <- cbind(protein, group2_pooled_sd)

# Creating a temp data frame to calculate the ratio3 sd
temp_df <- NULL 
sd_sum <- c("group3_sd_df", "ctrl_sd_df", "group3_df", "ctrl_df")
temp_df <- by_protein(data_top, sd_sum, FUN = sum) %>% as.data.frame

group3_pooled_sd <- temp_df$group3_sd_df / temp_df$group3_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
protein <- cbind(protein, group3_pooled_sd, ctrl_pooled_sd)

# Grouped relative protein abundance ratio
protein$ratio1 <- by_protein(data, "ratio1") %>% as.numeric
protein$ratio2 <- by_protein(data, "ratio2") %>% as.numeric
protein$ratio3 <- by_protein(data, "ratio3") %>% as.numeric

# Standard deviation relative protein abundance ratio
data$ratio1_sd <- ratio_sd(data, "ratio1", "group1_sd", "group1_med", 
                          "ctrl_sd", "ctrl_med")
data$ratio2_sd <- ratio_sd(data, "ratio2", "group2_sd", "group2_med", 
                          "ctrl_sd", "ctrl_med")
data$ratio3_sd <- ratio_sd(data, "ratio3", "group3_sd", "group3_med", 
                          "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio1
data$ratio1_df <- data$group1_df + data$ctrl_df
data$ratio1_sd_df <- square_x_df(data, "ratio1_sd", "ratio1_df")

# Standard deviation grouped relative protein abundance ratio2
data$ratio2_df <- data$group2_df + data$ctrl_df
data$ratio2_sd_df <- square_x_df(data, "ratio2_sd", "ratio2_df")

# Standard deviation grouped relative protein abundance ratio3
data$ratio3_df <- data$group2_df + data$ctrl_df
data$ratio3_sd_df <- square_x_df(data, "ratio3_sd", "ratio3_df")

# ratio1 sd_calc
sd_calc <- c("ratio1_df", "ratio1_sd_df")
temp_df <- by_protein(data, sd_calc, FUN = sum) %>% as.data.frame

ratio1_sd <- temp_df$ratio1_sd_df / temp_df$ratio1_df
protein <- cbind(protein, ratio1_sd)

# ratio2 sd_calc
sd_calc <- c("ratio2_df", "ratio2_sd_df")
temp_df <- by_protein(data, sd_calc, FUN = sum) %>% as.data.frame

ratio2_sd <- temp_df$ratio2_sd_df / temp_df$ratio2_df
protein <- cbind(protein, ratio2_sd)

# ratio3 sd_calc
sd_calc <- c("ratio3_df", "ratio3_sd_df")
temp_df <- by_protein(data, sd_calc, FUN = sum) %>% as.data.frame

ratio3_sd <- temp_df$ratio3_sd_df / temp_df$ratio3_df
protein <- cbind(protein, ratio3_sd)

# Statistics -------------------------------------------------------------------

log_norm <- log(data[, c(group1, group2, group3, ctrl)])
colnames(log_norm) <- paste(colnames(log_norm), sep = "_", "log")
data <- cbind(data, log_norm)

log_cols <- grep("log", colnames(data))
group1_log_cols <- grep("WT_12wk_norm_log", colnames(data))
group2_log_cols <- grep("KO_6wk_norm_log", colnames(data))
group3_log_cols <- grep("KO_12wk_norm_log", colnames(data))
ctrl_log_cols <- grep("WT_6wk_norm_log", colnames(data))


# Problem if different number of samples be in each group
length(group1_log_cols)
length(group2_log_cols)
length(group3_log_cols)
length(ctrl_log_cols)

# add 0 dummy columns to WT_12wk (1 for M, 0 for F)
# data$dummy1_WT_12wk_norm_log <- NA
# data$dummy2_WT_12wk_norm_log <- NA
# data$dummy3_WT_12wk_norm_log <- NA
# data$dummy4_WT_12wk_norm_log <- NA
# data$dummy5_WT_12wk_norm_log <- NA
# data$dummy6_WT_12wk_norm_log <- NA
# data$dummy7_WT_12wk_norm_log <- NA
# data$dummy8_WT_12wk_norm_log <- NA
# data$dummy9_WT_12wk_norm_log <- NA
# data$dummy10_WT_12wk_norm_log <- NA
# data$dummy11_WT_12wk_norm_log <- NA
# data$dummy12_WT_12wk_norm_log <- NA

# add 10 dummy columns to KO_6wk (6 for M, 5 for F)
data$dummy1_KO_6wk_norm_log <- NA
data$dummy2_KO_6wk_norm_log <- NA
data$dummy3_KO_6wk_norm_log <- NA
data$dummy4_KO_6wk_norm_log <- NA
data$dummy5_KO_6wk_norm_log <- NA
# data$dummy6_KO_6wk_norm_log <- NA
# data$dummy7_KO_6wk_norm_log <- NA
# data$dummy8_KO_6wk_norm_log <- NA
# data$dummy9_KO_6wk_norm_log <- NA
# data$dummy10_KO_6wk_norm_log <- NA
# data$dummy11_KO_6wk_norm_log <- NA
# data$dummy12_KO_6wk_norm_log <- NA

# add 2 dummy columns to KO_12wk (0 for M, 3 for F)
data$dummy1_KO_12wk_norm_log <- NA
data$dummy2_KO_12wk_norm_log <- NA
data$dummy3_KO_12wk_norm_log <- NA
# data$dummy4_KO_12wk_norm_log <- NA
# data$dummy5_KO_12wk_norm_log <- NA
# data$dummy6_KO_12wk_norm_log <- NA
# data$dummy7_KO_12wk_norm_log <- NA
# data$dummy8_KO_12wk_norm_log <- NA
# data$dummy9_KO_12wk_norm_log <- NA
# data$dummy10_KO_12wk_norm_log <- NA
# data$dummy11_KO_12wk_norm_log <- NA
# data$dummy12_KO_12wk_norm_log <- NA
# data$dummy13_KO_12wk_norm_log <- NA
# data$dummy14_KO_12wk_norm_log <- NA
# data$dummy15_KO_12wk_norm_log <- NA
# data$dummy16_KO_12wk_norm_log <- NA
# data$dummy17_KO_12wk_norm_log <- NA
# data$dummy18_KO_12wk_norm_log <- NA

# add 8 dummy columns to WT_6wk (6 for M, 3 for F)
data$dummy1_WT_6wk_norm_log <- NA
data$dummy2_WT_6wk_norm_log <- NA
data$dummy3_WT_6wk_norm_log <- NA
# data$dummy4_WT_6wk_norm_log <- NA
# data$dummy5_WT_6wk_norm_log <- NA
# data$dummy6_WT_6wk_norm_log <- NA
# data$dummy7_WT_6wk_norm_log <- NA
# data$dummy8_WT_6wk_norm_log <- NA
# data$dummy9_WT_6wk_norm_log <- NA
# data$dummy10_WT_6wk_norm_log <- NA
# data$dummy11_WT_6wk_norm_log <- NA
# data$dummy12_WT_6wk_norm_log <- NA
# data$dummy13_WT_6wk_norm_log <- NA
# data$dummy14_WT_6wk_norm_log <- NA
# data$dummy15_WT_6wk_norm_log <- NA
# data$dummy16_WT_6wk_norm_log <- NA
# data$dummy17_WT_6wk_norm_log <- NA
# data$dummy18_WT_6wk_norm_log <- NA
# data$dummy19_WT_6wk_norm_log <- NA
# data$dummy20_WT_6wk_norm_log <- NA

group1_log_cols <- grep("WT_12wk_norm_log", colnames(data))
group2_log_cols <- grep("KO_6wk_norm_log", colnames(data))
group3_log_cols <- grep("KO_12wk_norm_log", colnames(data))
ctrl_log_cols <- grep("WT_6wk_norm_log", colnames(data))

stacked <- data.frame(data[, "Master.Protein.Accessions"], 
                      stack(data[, group1_log_cols]),
                      stack(data[, group2_log_cols]),
                      stack(data[, group3_log_cols]),
                      stack(data[, ctrl_log_cols]))


colnames(stacked) <- c("Master.Protein.Accessions", "group1_log", "samp1", 
                       "group2_log", "samp2","group3_log", "samp3",
                       "ctrl_log", "ctrl")

# t-test for group1
p_vals1 <- pval_ttest(stacked, "group1_log", "ctrl_log")
colnames(p_vals1) <- c("Master.Protein.Accessions", "pvalue1")

p_vals1$Master.Protein.Accessions <-
  p_vals1$Master.Protein.Accessions %>% 
  as.character

protein <- full_join(protein, p_vals1, by = "Master.Protein.Accessions")

# t-test for group2
p_vals2 <- pval_ttest(stacked, "group2_log", "ctrl_log")
colnames(p_vals2) <- c("Master.Protein.Accessions", "pvalue2")

p_vals2$Master.Protein.Accessions <-
  p_vals2$Master.Protein.Accessions %>% 
  as.character

protein <- full_join(protein, p_vals2, by = "Master.Protein.Accessions")

# t-test for group3
p_vals3 <- pval_ttest(stacked, "group3_log", "ctrl_log")
colnames(p_vals3) <- c("Master.Protein.Accessions", "pvalue3")

p_vals3$Master.Protein.Accessions <-
  p_vals3$Master.Protein.Accessions %>% 
  as.character

protein <- full_join(protein, p_vals3, by = "Master.Protein.Accessions")

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv("6wk and 12wk Ventricles with August Reps Gene list.csv")

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 3 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "Female only SNORD KO WRT WT 6wk.csv")



sixwk <- read.csv("Female only SNORD KO WRT WT 6wk reduced.csv")
twelvewk <- read.csv("Female only SNORD KO WRT WT 12wk ratio only.csv")

fullsum <- full_join(sixwk, twelvewk, by = c("gene", "Master.Protein.Accessions"))
write.csv(fullsum, "Female only SNORD KO summary.csv")
