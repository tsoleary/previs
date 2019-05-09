# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------
# HCM DIS GOL samples

library(tidyverse)

setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")
data_raw <- read.csv("DIS_GOL_proteomics_all_pep.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("GOL", colnames(data_raw))

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
myosin <- "A5YM51; P12883; P13533"
histones <- "B4DR52; P06899"

norm_pro <- myosin

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

group1 <- grep("DIS_1_norm", colnames(data))
group2 <- grep("DIS_2_norm", colnames(data))
group3 <- grep("DIS_3_norm", colnames(data))
ctrl <- grep("GOL_norm", colnames(data))

data$group1_med <- by_group(data, group1)
data$group2_med <- by_group(data, group2)
data$group3_med <- by_group(data, group3)
data$ctrl_med <- by_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- by_group(data, group1, FUN = sd)
data$group2_sd <- by_group(data, group2, FUN = sd)
data$group3_sd <- by_group(data, group3, FUN = sd)
data$ctrl_sd <- by_group(data, ctrl, FUN = sd)

# Relative abundance ratio for each peptide
abun_ratio <- function (dat, group, ctrl = "ctrl_med"){
  dat[, group] / dat[, ctrl]
}

data$ratio1 <- abun_ratio(data, "group1_med")
data$ratio2 <- abun_ratio(data, "group2_med")
data$ratio3 <- abun_ratio(data, "group3_med")

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
data$group2_df <- count_df(data, group2)
data$group3_df <- count_df(data, group3)
data$ctrl_df <- count_df(data, ctrl)

# Removing rows with NA values for the median
data <- data[!(is.na(data$ratio1)), ]
data <- data[!(is.na(data$ratio2)), ]
data <- data[!(is.na(data$ratio3)), ]

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

# Remove outlier without function #######ratio 1################################

pro_out1 <- by_protein(data, "ratio1") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

sd_ratio_temp_df <- by_protein(data, "ratio1", FUN = sd) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions") %>%
  'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))

pro_temp <- dplyr::full_join(pro_out1, sd_ratio_temp_df,
                             by = "Master.Protein.Accessions")

pro_temp$max_ratio <- pro_temp$ratio1 + 2 * pro_temp$sd_ratios
pro_temp$min_ratio <- pro_temp$ratio1 - 2 * pro_temp$sd_ratios

data_rm_out <- NULL

for (pro in unique(data$Master.Protein.Accessions)){
  temp <- dplyr::filter(data, data$Master.Protein.Accessions == pro)
  
  rm_high <- which(temp$ratio1 > pro_temp$max_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm_low <- which(temp$ratio1 < pro_temp$min_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm <- c(rm_high, rm_low)
  if (length(rm) > 0){
    temp_rm <- temp[-rm, ]
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp_rm)
  } else {
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp)
  }
}

data <- data_rm_out

################################################################################

# Remove outlier without function #######ratio 2################################

pro_out2 <- by_protein(data, "ratio2") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

sd_ratio_temp_df <- by_protein(data, "ratio2", FUN = sd) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions") %>%
  'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))

pro_temp <- dplyr::full_join(pro_out2, sd_ratio_temp_df,
                             by = "Master.Protein.Accessions")

pro_temp$max_ratio2 <- pro_temp$ratio2 + 2 * pro_temp$sd_ratios
pro_temp$min_ratio2 <- pro_temp$ratio2 - 2 * pro_temp$sd_ratios

data_rm_out <- NULL

for (pro in unique(data$Master.Protein.Accessions)){
  temp <- dplyr::filter(data, data$Master.Protein.Accessions == pro)
  
  rm_high <- which(temp$ratio2 > pro_temp$max_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm_low <- which(temp$ratio2 < pro_temp$min_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm <- c(rm_high, rm_low)
  if (length(rm) > 0){
    temp_rm <- temp[-rm, ]
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp_rm)
  } else {
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp)
  }
}

data <- data_rm_out

################################################################################

# Remove outlier without function #######ratio 3################################

pro_out3 <- by_protein(data, "ratio3") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

sd_ratio_temp_df <- by_protein(data, "ratio3", FUN = sd) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions") %>%
  'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))

pro_temp <- dplyr::full_join(pro_out3, sd_ratio_temp_df,
                             by = "Master.Protein.Accessions")

pro_temp$max_ratio3 <- pro_temp$ratio3 + 2 * pro_temp$sd_ratios
pro_temp$min_ratio3 <- pro_temp$ratio3 - 2 * pro_temp$sd_ratios

data_rm_out <- NULL

for (pro in unique(data$Master.Protein.Accessions)){
  temp <- dplyr::filter(data, data$Master.Protein.Accessions == pro)
  
  rm_high <- which(temp$ratio3 > pro_temp$max_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm_low <- which(temp$ratio3 < pro_temp$min_ratio[which(
    pro_temp$Master.Protein.Accessions == pro)])
  
  rm <- c(rm_high, rm_low)
  if (length(rm) > 0){
    temp_rm <- temp[-rm, ]
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp_rm)
  } else {
    data_rm_out <- dplyr::bind_rows(data_rm_out, temp)
  }
}

data <- data_rm_out

################################################################################


# with the rm_outliers function ################################################

pro_out1 <- by_protein(data, "ratio1") %>% as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

data <- rm_outliers(data, pro_out1, "ratio1")

# still might not work #########################################################


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

square_x_df <- function (dat, group_sd, group_df){
  (dat[, group_sd])^2 * (dat[, group_df])
}

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
group1_log_cols <- grep("DIS_1_norm_log", colnames(data))
group2_log_cols <- grep("DIS_2_norm_log", colnames(data))
group3_log_cols <- grep("DIS_3_norm_log", colnames(data))
ctrl_log_cols <- grep("GOL_norm_log", colnames(data))


# Problem if different number of samples be in each group
length(group1_log_cols)
length(group2_log_cols)
length(group3_log_cols)
length(ctrl_log_cols)

# add 6 dummy columns to DIS_2
data$dummy1_DIS_2_norm_log <- NA
data$dummy2_DIS_2_norm_log <- NA
data$dummy3_DIS_2_norm_log <- NA
data$dummy4_DIS_2_norm_log <- NA
data$dummy5_DIS_2_norm_log <- NA
data$dummy6_DIS_2_norm_log <- NA

# add 9 dummy columns to DIS_3
data$dummy1_DIS_3_norm_log <- NA
data$dummy2_DIS_3_norm_log <- NA
data$dummy3_DIS_3_norm_log <- NA
data$dummy4_DIS_3_norm_log <- NA
data$dummy5_DIS_3_norm_log <- NA
data$dummy6_DIS_3_norm_log <- NA
data$dummy7_DIS_3_norm_log <- NA
data$dummy8_DIS_3_norm_log <- NA
data$dummy9_DIS_3_norm_log <- NA

# add 12 dummy columns to DIS_3
data$dummy1_GOL_norm_log <- NA
data$dummy2_GOL_norm_log <- NA
data$dummy3_GOL_norm_log <- NA
data$dummy4_GOL_norm_log <- NA
data$dummy5_GOL_norm_log <- NA
data$dummy6_GOL_norm_log <- NA
data$dummy7_GOL_norm_log <- NA
data$dummy8_GOL_norm_log <- NA
data$dummy9_GOL_norm_log <- NA
data$dummy10_GOL_norm_log <- NA
data$dummy11_GOL_norm_log <- NA
data$dummy12_GOL_norm_log <- NA


group1_log_cols <- grep("DIS_1_norm_log", colnames(data))
group2_log_cols <- grep("DIS_2_norm_log", colnames(data))
group3_log_cols <- grep("DIS_3_norm_log", colnames(data))
ctrl_log_cols <- grep("GOL_norm_log", colnames(data))

stacked <- data.frame(data[, "Master.Protein.Accessions"], 
                      stack(data[, group1_log_cols]),
                      stack(data[, group2_log_cols]),
                      stack(data[, group3_log_cols]),
                      stack(data[, ctrl_log_cols]))


colnames(stacked) <- c("Master.Protein.Accessions", "group1_log", "samp1", 
                       "group2_log", "samp2","group3_log", "samp3",
                       "ctrl_log", "ctrl")

pval_ttest <- function (dat, group, ctrl, col = "Master.Protein.Accessions"){
  name <- NULL
  p_value <- NULL
  for (pro in unique(dat[, col])) {
    temp <- dplyr::filter(dat, dat[, col] == pro)
    pval_temp <- tryCatch(t.test(temp[, group], temp[, ctrl])$p.value, 
                          error=function(err) NA)
    name <- c(name, pro)
    p_value <- c(p_value, pval_temp)
  }
  return(as.data.frame(cbind(name, p_value)))
}

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
gene_df <- read.csv("human_PD_accession_gene.csv")

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

write.csv(protein, "HCM_DIS_GOL_myo_norm.csv")

