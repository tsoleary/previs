# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)

setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")
data_raw <- read.csv("WT vs KO_pep.csv")

# Normalization ----------------------------------------------------------------

# Function applied to each peptide for a group of data with median as default
as_group <- function (dat, col, FUN = median){
  list <- NULL
  for (i in 1:nrow(dat)){
    temp <- FUN(as.numeric(dat[i, col], na.rm = TRUE))
    list <- c(list, temp)
  }
  return(list)
}

ctrl_raw <- grep("Control", colnames(data_raw))
data_raw$ctrl_raw_med <- as_group(data_raw, ctrl_raw)

# set the max number of peptides used in analysis
max_pep <- 15

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# proteins used for normalization
norm_pro <- "B2RQQ1; Q91Z83"

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

group1 <- grep("Sample_norm", colnames(data))
ctrl <- grep("Control_norm", colnames(data))

data$group1_med <- as_group(data, group1)
data$ctrl_med <- as_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- as_group(data, group1, FUN = sd)
data$ctrl_sd <- as_group(data, ctrl, FUN = sd)

# Relative abundance ratio for each peptide
ratio <- data$group1_med / data$ctrl_med
data$ratio <- ratio

# Function to count degrees of freedom in a group, defaults witout duplicates
count_df <- function (dat, col, dup = 1){
  list <- NULL
  for (i in 1:nrow(dat)){
    temp <- length(which(!is.na(dat[i, col])))
    temp <- temp / dup
    if(temp >= 1){
      temp <- temp - 1
    }
    list <- c(list, temp)
  }
  return(list)
}

# Degrees of freedom
data$group1_df <- count_df(data, group1, dup = 3)
data$ctrl_df <- count_df(data, ctrl, dup = 3)

# Removing rows with NA values for the median
data <- data[!(is.na(data$ratio)), ]

# Data frame with only top few ionizing peptides
pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

# Relative abundance of each protein using top ionizers
protein_group <- function (dat, groups, FUN = mean){
  tab <- NULL
  for (col in groups){
    temp <- tapply(dat[, col],
                   dat$Master.Protein.Accessions,
                   FUN,
                   na.rm = TRUE)
    tab <- cbind(tab, temp)
  }
  colnames(tab) <- groups
  return(tab)
}

group_names <- c("group1_med", "ctrl_med")

protein <- protein_group(data_top, group_names) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

protein$Master.Protein.Accessions <-
  protein$Master.Protein.Accessions %>%
  as.character

# Standard deviation of grouped relative abundance using top ionizers
square_x_df <- function (dat, group_sd, group_df){
  (dat[, group_sd])^2 * (dat[, group_df])
}

data_top$group1_sd_df <- square_x_df(data_top, "group1_sd", "group1_df")
data_top$ctrl_sd_df <- square_x_df(data_top, "ctrl_sd", "ctrl_df")

# Creating a temp data frame to calculate the ratio sd

temp_df <- NULL # initialize data frame
sd_sum <- c("group1_sd_df", "ctrl_sd_df", "group1_df", "ctrl_df")
temp_df <- protein_group(data_top, sd_sum, FUN = sum)
temp_df <- as.data.frame(temp_df)

group1_pooled_sd <- temp_df$group1_sd_df / temp_df$group1_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
protein <- cbind(protein, group1_pooled_sd, ctrl_pooled_sd)

# Grouped relative protein abundance ratio
protein$ratio <- protein_group(data, "ratio")

# Taylor Expansion to get the sd of the ratio of two means
ratio_sd <- function (dat, ratio, group1_sd, group1_med, ctrl_sd, ctrl_med){
  ratio_sd_pep <- NULL
  for (i in 1:nrow(dat)){
    temp1 <- (dat[i, group1_sd] / dat[i, group1_med])^2
    temp2 <- (dat[i, ctrl_sd] / dat[i, ctrl_med])^2
    result <- dat[i, ratio] * sqrt(temp1 + temp2)
    ratio_sd_pep <- c(ratio_sd_pep, result)
  }
  return(ratio_sd_pep)
}

# Standard deviation relative protein abundance ratio
data$ratio_sd <- ratio_sd(data, "ratio", "group1_sd", "group1_med", 
                             "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio
data$ratio_df <- data$group1_df + data$ctrl_df
data$ratio_sd_df <- (data$ratio_sd_pep)^2 * data$ratio_df

temp_df <- NULL

ratio_df_sum <- tapply(data$ratio_df,
                       data$Master.Protein.Accessions,
                       sum,
                       na.rm = TRUE)

ratio_sd_df_sum <- tapply(data$ratio_sd_df,
                       data$Master.Protein.Accessions,
                       sum,
                       na.rm = TRUE)

temp_df <- cbind(ratio_df_sum, ratio_sd_df_sum)
temp_df <- as.data.frame(temp_df)
ratio_sd <- temp_df$ratio_sd_df_sum / temp_df$ratio_df_sum

protein_table <- cbind(protein_table, ratio_sd)
colnames(protein_table) <- c(pro_med, sd_calc, pro_data, "ratio_sd")

# Statistics -------------------------------------------------------------------
log_norm <- log(data[, c(group1, ctrl)])
colnames(log_norm) <- paste(colnames(log_norm), sep = "_", "log")
data <- cbind(data, log_norm)

log_cols <- grep("log", colnames(data))
group1_log_cols <- grep("Sample_norm_log", colnames(data))
ctrl_log_cols <- grep("Control_norm_log", colnames(data))

group1_stack <- data.frame(data[, 3], stack(data[, group1_log_cols]))
colnames(group1_stack) <- c("Master.Protein.Accessions", "group1_log", "samp")

ctrl_stack <- data.frame(data[, 3], stack(data[, ctrl_log_cols]))
colnames(ctrl_stack) <- c("Master.Protein.Accessions", "ctrl_log", "ctrl")


stacked <- cbind(group1_stack, ctrl_stack[, 2:3])

protein <- NULL
p_vals <- NULL
for (pro in unique(stacked$Master.Protein.Accessions)){
  temp <- filter(stacked, stacked$Master.Protein.Accessions == pro)
  p_val_temp <- t.test(temp$group1_log, temp$ctrl_log)$p.value
  protein <- c(protein, pro)
  p_vals <- c(p_vals, p_val_temp)
}

pro_pvals <- as.data.frame(cbind(protein, p_vals))
colnames(pro_pvals) <- c("Master.Protein.Accessions", "p-value")

protein_df <- as.data.frame(protein_table)

protein_df <- full_join(protein_df, pro_pvals, by = "Master.Protein.Accessions")
# Warning - joining character vector and factor
protein_df$Master.Protein.Accessions <-
  as.factor(protein_df$Master.Protein.Accessions)

# Minimum number of peptides for each protein group

min_pep <- 5 
protein_df$peptides <- table(data$Master.Protein.Accessions)
protein_min <- filter(protein_df, protein_df$peptides >= min_pep)

# Converting Protein Accession to Gene Symbol ----------------------------------
gene_df <- read.csv('mouse_PD_accession_gene.csv')

mpa_to_gene <- function (dat, gene_dat){
  dat$gene <- dat$Master.Protein.Accessions
  for (i in 1:nrow(dat)){
    temp <- which(dat$gene[i] == gene_dat$Accession, TRUE)
    if (length(temp) == 1){
      dat$gene <- gsub(dat$gene[i], gene_dat$Gen[temp], dat$gene)
    }
  }
  return(dat$gene)
}

data$gene <- mpa_to_gene(data, gene_df)
protein_df$gene <- mpa_to_gene(protein_df, gene_df)

# Remove outliers --------------------------------------------------------------

# Removing rows with NA values for the median
data <- data[!(...), ]

#home
