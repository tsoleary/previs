# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------
# WT v MyBP-C KO Mouse samples

library(tidyverse)

# source the functions in the proteomics_functions.R script!

setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")
data_raw <- read.csv("WT vs KO all pep.csv")

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
myosin <- "B2RQQ1; Q91Z83"
histones <- "B2RTM0"

norm_pro <- histones

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

# Median, sd, & ratio of peptides ----------------------------------------------

group1 <- grep("Sample_norm", colnames(data))
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
group1_log_cols <- grep("Sample_norm_log", colnames(data))
ctrl_log_cols <- grep("Control_norm_log", colnames(data))

# Problem if different number of samples be in each group ----------------------

length(group1_log_cols) - length(ctrl_log_cols) == 0

group1_log_cols <- grep("Sample_norm_log", colnames(data))
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
gene_df <- read.csv("mouse_PD_accession_gene.csv")

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 5 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "WT_v_KO_hist_norm.csv")

################################################################################
## STOP ------------------------------------------------------------------------

# Group by sub-cellular compartment --------------------------------------------

sub_cell_comp <- read.csv("sub_cell_comp_mouse.csv")
sub_cell_comp <- subset(sub_cell_comp, !duplicated(sub_cell_comp$Gene))

# Group genes into subcellular compartments on data
data$compartment <- gene_to_comp2(data, sub_cell_comp,
                                  level = "compartment")

# data$sub_compartment <- gene_to_comp(data, sub_cell_comp, 
#                                      level = "sub_compartment")

# There is some sort of bug in the gene_to_comp function that results in this:
# > data$compartment[1]
# [1] "MitochondriaMMitochondriaiMitochondriatMitochondriaoMitochondriac
# MitochondriahMitochondriaoMitochondrianMitochondriadMitochondriarMitochondriai
# MitochondriaaMitochondria"

gene_to_comp2 <- function (dat, comp_dat, level = "compartment") {
  dat$comp <- dat$gene
  for (i in 1:nrow(dat)) {
    temp <- which(dat$gene[i] == comp_dat$Gene, TRUE)
    if (length(temp) == 1) {
      dat$comp[i] <- gsub(dat$comp[i], comp_dat[temp, level], 
                          dat$comp, ignore.case = TRUE)
    }
  }
  return(dat$comp)
}

# test code below

data$comp <- data$gene
which(data$gene[1] == sub_cell_comp$Gene, TRUE)


data$compartment <- gene_to_comp(data, sub_cell_comp,
                                 level = "compartment")




################################################################################

# Condense unique compartments into only one
compart_list <- as.character(unique(sub_cell_comp$compartment))

comp_simp <- data$compartment
list <- NULL
temp <- NULL
for (pro in compart_list) {
  for (i in 1:nrow(data)){
    if (str_detect(comp_simp[i], pro) == TRUE){
      temp <- str_extract(comp_simp[i], pro)
    } else {
      comp_simp[i] <- comp_simp[i]
    }
  }
}

data$compartment <- comp_simp

# Compartments ratio average
weighted_ratio <- cbind(sapply(split(data, data$compartment), 
                               function (x) {weighted.mean(x$ratio, x$group1_med)}))

# Statistics for compartments and sub-compartments------------------------------

stack2 <- data.frame(data[, "compartment"], stack(data[, group1_log_cols]),
                     stack(data[, ctrl_log_cols]))
colnames(stack2) <- c("compartment", "group1_log", "samp", "ctrl_log", "ctrl")

p_val2 <- pval_ttest(stack2, "group1_log", "ctrl_log", col = "compartment")
colnames(p_val2)[1] <- "compartment"