# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Previs")
data_raw <- 
  read.csv("Extraction WT 150 200 250 All peptides.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("F", colnames(data_raw))

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

data_raw <- filter(data_raw, (is.na(ctrl_raw_med) == FALSE))

# set the max number of peptides used in analysis
max_pep <- 5

data <-
  as_tibble(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)


# proteins used for normalization
BSA <- "P02769"


# norm_pro <- titin
norm_pro <- BSA

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
gene_df <- read.csv('Extraction WT 150 200 250 Allreps Gene List.csv')

data$gene <- mpa_to_gene(data, gene_df)
protein$gene <- mpa_to_gene(protein, gene_df)

# Minimum number of peptides for each protein group ----------------------------
min_pep <- 2
protein$peptides <- table(data_raw$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

# Summary Sheet ----------------------------------------------------------------

# Gather separate file abundances to single column, and separate filename into
# several variables
pbyfile <- gather(protein, contains("F"), key = "file", value = "abundance")
septest <- separate(pbyfile, file, into = c("file", "frac","conc","salt","type"))

# Assign replicate-number by PD file-identifier
indlist <- read.csv("WT 150 200 250 Allreps 1209 identifiers.csv")
septest$rep <- file_to_ind(septest, indlist)

# Replace NA, NaN values with zeroes
septest$abundance[is.na(septest$abundance)] = 0

# Spread fraction variable into separate abundance columns "A", "B", "P", and 
# join together while dropping nonessential variables
A <- filter(septest, frac == "A") %>%
  select(., c(1:3, 6, 10, 9)) %>%
  rename(., A = abundance)
B <- filter(septest, frac == "B") %>%
  select(., c(1:3, 6, 10, 9)) %>%
  rename(., B = abundance)
P <- filter(septest, frac == "P") %>%
  select(., c(1:3, 6, 10, 9)) %>%
  rename(., P = abundance)

jointest <- full_join(A, B) %>%
  full_join(., P)

# Calculate soluble and insoluble fraction within each conc., rep., and gene
jointest <- mutate(jointest, total = A+B+P)
jointest <- mutate(jointest, insoluble = P/total)
jointest <- mutate(jointest, soluble = (A+B)/total)

# Export values in current state (no averaging)
write.csv(jointest, "extraction individual reps test.csv")

# Compute summary stats within conc. and gene

solsummary <- group_by(jointest, conc, gene) %>%
  summarise(., solavg = mean(soluble, na.rm = TRUE), solsd = sd(soluble, na.rm = TRUE),
            peptides = mean(peptides)) %>%
  ungroup(.)


# mean <- group_by(jointest, conc, gene) %>%
#   summarise(., solavg = mean(soluble, na.rm = TRUE))
# sd <- group_by(jointest, conc, gene) %>%
#   summarise(., solsd = sd(soluble, na.rm = TRUE))
# solsummary <- full_join(mean, sd) %>%
#   ungroup(.)

# Spread conc. variable into multiple mean columns

Low <- filter(solsummary, conc == 150) %>%
  select(-conc) %>%
  rename(., avg150 = solavg, sd150 = solsd)
Mid <- filter(solsummary, conc == 200) %>%
  select(-conc) %>%
  rename(., avg200 = solavg, sd200 = solsd)
High <- filter(solsummary, conc == 250) %>%
  select(-conc) %>%
  rename(., avg250 = solavg, sd250 = solsd)

Concsum <- full_join(Low, Mid) %>%
  full_join(., High) %>%
  mutate(., "200/150" = avg200/avg150, "250/150" = avg250/avg150)

# Export -----------------------------------------------------------------------

write.csv(protein, "Proteins.csv")
# write.csv(specialprotein, "H4 norm protein watchlist by Top3.csv")
write.csv(data, "Peptides.csv")
# write.csv(specialpeptide, "H4 Norm watchlist peptides Top5.csv")
write.csv(Concsum, "Full extraction summary 150 200 250 WT with peptide counts.csv")

