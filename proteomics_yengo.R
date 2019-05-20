# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------
# Yengo

library(tidyverse)

setwd("C:/Users/PrevBeast/Documents/R/Yengo")
data_raw <- read.csv("Yengo_2_6_12_hr_infected_051719_all_peptides.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("Control", colnames(data_raw))

by_group <- function (dat, col, FUN = median) {
  list <- NULL
  for (i in 1:nrow(dat)) {
    temp <- FUN(as.numeric(dat[i, col]), na.rm = TRUE)
    list <- c(list, temp)
  }
  return(list)
}

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

# set the max number of peptides used in analysis
max_pep <- 15

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# proteins used for normalization
histones <- "B2RTM0"

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

group1 <- grep("infected_norm", colnames(data))
ctrl <- grep("Control_norm", colnames(data))

data$group1_med <- by_group(data, group1)
data$ctrl_med <- by_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- by_group(data, group1, FUN = sd)
data$ctrl_sd <- by_group(data, ctrl, FUN = sd)

# Relative abundance ratio for each peptide
abun_ratio <- function (dat, group, ctrl = "ctrl_med"){
  dat[, group] / dat[, ctrl]
}

data$ratio <- abun_ratio(data, "group1_med")

# Degrees of freedom
count_df <- function (dat, col, rep = 1){
  list <- NULL
  for (i in 1:nrow(dat)){
    temp <- length(which(!is.na(dat[i, col])))
    temp <- temp / rep
    if(temp >= 1){
      temp <- temp - 1
    }
    list <- c(list, temp)
  }
  return(list)
}

data$group1_df <- count_df(data, group1)
data$ctrl_df <- count_df(data, ctrl)

# Removing rows with NA values for the median
data <- data[!(is.na(data$ratio)), ]

# Remove outliers --------------------------------------------------------------
by_protein <- function (dat, groups, FUN = mean){
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

pro_out <- by_protein(data, "ratio") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

rm_outliers <- function (dat, pro_df, ratio, mult = 2){
  
  sd_ratio_temp_df <- by_protein(dat, ratio, FUN = sd) %>%
    as.data.frame %>%
    rownames_to_column("Master.Protein.Accessions") %>%
    'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))
  
  pro_temp <- dplyr::full_join(pro_df, sd_ratio_temp_df,
                               by = "Master.Protein.Accessions")
  
  pro_temp$max_ratio <- pro_temp$ratio + mult * pro_temp$sd_ratio
  pro_temp$min_ratio <- pro_temp$ratio - mult * pro_temp$sd_ratio
  
  data_rm_out <- NULL
  
  for (pro in unique(dat$Master.Protein.Accessions)){
    temp <- dplyr::filter(dat, dat$Master.Protein.Accessions == pro)
    
    rm_high <- which(temp[, ratio] > pro_temp$max_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)])
    
    rm_low <- which(temp[, ratio] < pro_temp$min_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)])
    
    rm <- c(rm_high, rm_low)
    if (length(rm) > 0){
      temp_rm <- temp[-rm, ]
      data_rm_out <- dplyr::bind_rows(data_rm_out, temp_rm)
    } else {
      data_rm_out <- dplyr::bind_rows(data_rm_out, temp)
    }
  }
  return(data_rm_out)
}


data <- rm_outliers(data, pro_out, "ratio")

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

square_x_df <- function (dat, group_sd, group_df){
  (dat[, group_sd])^2 * (dat[, group_df])
}

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

data$ratio_sd <- ratio_sd(data, "ratio", "group1_sd", "group1_med", 
                          "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio
data$ratio_df <- data$group1_df + data$ctrl_df
data$ratio_sd_df <- square_x_df(data, "ratio_sd", "ratio_df")

sd_calc <- c("ratio_df", "ratio_sd_df")
temp_df <- by_protein(data, sd_calc, FUN = sum) %>% as.data.frame

ratio_sd <- temp_df$ratio_sd_df / temp_df$ratio_df
protein <- cbind(protein, ratio_sd)

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv('mouse_uniprot_gene_yengo.csv')

mpa_to_gene <- function (dat, gene_dat){
  dat$gene <- dat$Master.Protein.Accessions
  for (i in 1:nrow(dat)){
    temp <- which(dat$gene[i] == gene_dat$Accession, TRUE)
    if (length(temp) == 1){
      dat$gene <- gsub(dat$gene[i], gene_dat$Gene[temp], dat$gene)
    }
  }
  return(dat$gene)
}

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 5 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "yengo_hist_norm.csv")