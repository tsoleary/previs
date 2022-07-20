# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Previs/D3L mouse study/Sumtotals")
data_raw <- 
  read.csv("All LV peptides day 56.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("F", colnames(data_raw))

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

data_raw <- filter(data_raw, (is.na(ctrl_raw_med) == FALSE))

# Specify possible dynamic modifications

dynmod <- c("Oxidation", "Phospho", "Dioxidation", "Acetyl", "Label:2H\\(3\\)")

# Merge strings in Annotated Sequence column with Master Protein Accessions
# (Accounts for special case when summing modifications, wherein peptides with
# the same exact sequence are matched to different accessions by PD)

for (i in 1:nrow(data_raw)) {
  if (i == 1) {
    temp <- str_c(data_raw$Annotated.Sequence[i], data_raw$Master.Protein.Accessions[i], sep = " ")
  }
  if (i > 1) {
    temp <- c(temp, str_c(data_raw$Annotated.Sequence[i], data_raw$Master.Protein.Accessions[i], sep = " "))
  }
}

data_raw$seqmpa <- temp

# Create logical columns to simplify information in "Modifications" column of data_raw
# Question: Does this row have modification X? (where X is any specified dynamic modification) ---------------------------------------------

# temp <- NULL
# j <- 14
# i <- 1
# 
# rm(temp, j, i)

for (i in 1:length(dynmod)) {
  temp <- NULL
  hashap <- FALSE
  for (j in 1:nrow(data_raw)) {
    if (length(temp) == 0 & !str_detect(data_raw$Modifications[[j]], dynmod[i])) {
      # print(paste("First row is a dud! No", dynmod[i], "at row", j))
      temp <- FALSE
      # print(j)
      # print(length(temp))
      next
    }
    if (length(temp) == 0 & str_detect(data_raw$Modifications[[j]], dynmod[i])) {
      # print(paste("Wow, first try!", dynmod[i], "at row", j))
      temp <- TRUE
      # print(j)
      # print(length(temp))
      next
    }
    if (length(temp) > 0 & !str_detect(data_raw$Modifications[[j]], dynmod[i])) {
      # print(paste("Blimey! No", dynmod[i], "at row", j))
      temp <- c(temp, FALSE)
      if (j != length(temp) & hashap == FALSE){
        print(j)
        print(length(temp))
        hashap <- TRUE
      }
      next
    }
    if (length(temp) > 0 & str_detect(data_raw$Modifications[[j]], dynmod[i])) {
      # print(paste("Blimey! We have", dynmod[i], "at row", j))
      temp <- c(temp, TRUE)
      if (j != length(temp) & hashap == FALSE){
        print(j)
        print(length(temp))
        hashap <- TRUE
      }
      next
    }
  }
  print(dynmod[i])
  testcol <- dynmod[i]
  data_raw[, testcol] <- temp
  # if (i == 1) {
  #   
  # }
  # name(temp) <- testcol
  # assign(testcol, temp)
  # print(temp)
  print("next mod please")
}

# Filter out rows containing undesired dynamic modifications -------------------------------------------------------------------------------

wantmod <- "Label:2H\\(3\\)"
dynmod[which(wantmod == dynmod)]

for (i in 1:length(dynmod)) {
  modname <- dynmod[i]
  tar <- grep(modname, colnames(data_raw))
  if (i == 1) {
    data_filt <- data_raw
    if (dynmod[i] != wantmod) {
      data_filt <- data_filt[data_filt[,tar] == FALSE,]
    }
    next
  }
  if (dynmod[i] != wantmod) {
    data_filt <- data_filt[data_filt[,tar] == FALSE,]
  } else {
    next
  }
}

# Sum all peptides with same primary sequence, such that all non-filtered,
# modified states are counted in a single row -----------------------------------------------------------------------------------------------

data_sum <- by_variable(data_filt, colnames(data_filt)[ctrl_raw], FUN = sum, VAR = "seqmpa") %>%
  as.data.frame %>%
  rownames_to_column("seqmpa")

# Get median of summed values for top_n operation to follow

ctrl_sum <- grep("F", colnames(data_sum))

data_sum$ctrl_sum_med <- by_group(data_sum, ctrl_sum)

# Re-separate Annotated Sequence and Master Protein Accessions from merged string
# column

str_split_fixed(c("cool beans", "cool beats", "cool fleans"), " ", 2)

data_sum$Annotated.Sequence <- str_split_fixed(data_sum$seqmpa, " ", 2)[,1]
data_sum$Master.Protein.Accessions <- str_split_fixed(data_sum$seqmpa, " ", 2)[,2]

# set the max number of peptides used in analysis
# max_pep <- 3
# 
# data <-
#   tbl_df(data_raw) %>%
#   group_by(Master.Protein.Accessions) %>%
#   top_n(n = max_pep, wt = ctrl_raw_med)

# Alternative for sum of filtered variants
max_pep <- 3

data <-
  tbl_df(data_sum) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_sum_med)


# proteins used for normalization
# titin <- "D3DPG0"
# alphabeta_shared <- "P12883; P13533"


# norm_pro <- titin
# norm_pro <- alphabeta_shared

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
# norm_value <- colSums(data_raw[, ctrl_raw], na.rm = TRUE)
# 
# raw_abun_mat <- as.matrix(data[, ctrl_raw])
# 
# norm_abun <- t(t(raw_abun_mat)/norm_value)
# colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
# norm_test <- as.data.frame(norm_abun)
# data <- as_tibble(data)
# data <- cbind(data, norm_test)

# Sum of columns, formatted for sum of modified states

norm_value <- colSums(data_raw[, ctrl_raw], na.rm = TRUE)

raw_abun_mat <- as.matrix(data[, ctrl_sum])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data <- as_tibble(data)
data <- cbind(data, norm_test)

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
gene_df <- read.csv('All LV gene list.csv')

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 1 
# protein$peptides <- table(data_raw$Master.Protein.Accessions)
protein$peptides <- table(data_sum$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "WT and KO d56 LV Top3 Sumtotal norm proteins.csv")
write.csv(data, "WT and KO d56 LV Top3 Sumtotal norm peptides.csv")
write.csv(norm_value, "sums.csv")
