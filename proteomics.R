# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
library(devtools)

# if proteomixr has been updated on github
# remove.packages("proteomixr") # remove old proteomixr
# install_github("tsoleary/proteomixr") # to get latest version
library("proteomixr")

#getwd()
setwd("C:/Users/PrevBeast/Documents/R/WT v KO mouse")
#setwd("C:/Users/PrevBeast/Documents/GitHub/previs")
data_raw <- read.csv("WT vs KO_pep.csv")

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

data$group1_med <- by_group(data, group1)
data$ctrl_med <- by_group(data, ctrl)

# Standard deviation between samples for each peptide
data$group1_sd <- by_group(data, group1, FUN = sd)
data$ctrl_sd <- by_group(data, ctrl, FUN = sd)

# Relative abundance ratio for each peptide
data$ratio <- abun_ratio(data, "group1_med")

# Degrees of freedom
data$group1_df <- count_df(data, group1, rep = 3)
data$ctrl_df <- count_df(data, ctrl, rep = 3)

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
group1_log_cols <- grep("Sample_norm_log", colnames(data))
ctrl_log_cols <- grep("Control_norm_log", colnames(data))

group1_stack <- data.frame(data[, 3], stack(data[, group1_log_cols]))
colnames(group1_stack) <- c("Master.Protein.Accessions", "group1_log", "samp")

ctrl_stack <- data.frame(data[, 3], stack(data[, ctrl_log_cols]))
colnames(ctrl_stack) <- c("Master.Protein.Accessions", "ctrl_log", "ctrl")

stacked <- cbind(group1_stack, ctrl_stack[, 2:3])

p_vals <- pval_ttest(stacked, "group1_log", "ctrl_log")
p_vals$Master.Protein.Accessions <-
  p_vals$Master.Protein.Accessions %>% 
    as.character

protein <- full_join(protein, p_vals, by = "Master.Protein.Accessions")

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv('mouse_PD_accession_gene.csv')

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 5 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

# Group by sub-cellular compartment --------------------------------------------

human_sub_cell <- read.csv("compartments_human.csv")

# Function from Character Translation and Casefolding
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

human_sub_cell$mouse_gene <- capwords(as.character(human_sub_cell$Gene), 
                                      strict = TRUE)

# Group genes into subcellular compartments
gene_to_group <- function (dat, group_dat){
  dat$Compartment <- dat$gene
  for (i in 1:nrow(dat)){
    temp <- which(dat$gene[i] == group_dat$mouse_gene, TRUE)
    if (length(temp) == 1){
      dat$Compartment <- gsub(dat$Compartment[i], group_dat[temp, "Compartment"],
                              dat$Compartment, ignore.case = TRUE)
    }
  }
  return(dat$compartment)
}

data$Compartment <- gene_to_group(data, human_sub_cell)

# Group and average on protein level

weighted.mean(c(2, 8, 1000), c(.99, .009, .001), na.rm = TRUE)


# Group and average on peptide level


