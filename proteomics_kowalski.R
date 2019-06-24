# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")
data_raw <- 
  read.csv("Kowalski_F_w1_w8_all_peptides.csv")

# Normalization ----------------------------------------------------------------

ctrl_raw <- grep("F", colnames(data_raw))

by_group <- function (dat, col, FUN = median) {
  list <- NULL
  for (i in 1:nrow(dat)) {
    temp <- FUN(as.numeric(dat[i, col]), na.rm = TRUE)
    list <- c(list, temp)
  }
  return(list)
}

data_raw$ctrl_raw_med <- by_group(data_raw, ctrl_raw)

# # set the max number of peptides used in analysis
# max_pep <- 100000
# 
# data <-
#   tbl_df(data_raw) %>%
#   group_by(Master.Protein.Accessions) %>%
#   top_n(n = max_pep, wt = ctrl_raw_med)

# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_raw_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- colnames(data_top)[grep("F", colnames(data))]

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

protein <- by_protein(data_top, group_names) %>%
  as.data.frame %>%
  rownames_to_column("Master.Protein.Accessions")

protein$Master.Protein.Accessions <-
  protein$Master.Protein.Accessions %>%
  as.character

# Converting protein accession to gene symbol ----------------------------------
gene_df <- read.csv('Kowalski_F_w1_w8_gene_names.csv')

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
min_pep <- 3 
protein$peptides <- table(data$Master.Protein.Accessions)
protein <- filter(protein, protein$peptides >= min_pep)

write.csv(protein, "kowalski_F_w1_w8_no_norm_r.csv")




################################################################################
# proteins changing over time --------------------------------------------------

protein <- read.csv("kowalski_F_w1_w8_no_norm_r.csv")
protein$X <- NULL

# convert protein data.frame to a tidy data.frame
samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "leg", "sex", "week"), sep = "_")

# average together the duplicates in the same week
df <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

df$week <- as.numeric(df$week)

df$group <- paste(df$Master.Protein.Accessions, df$sex, df$leg, sep = "_")

#remove the 
df <- df[!is.na(df$abundance), ]

# need to first make a plot of the whole thing
ggplot(df, aes(x = week, y = abundance)) +
  geom_jitter(mapping = aes(x = week, y = abundance, fill = leg), 
              alpha = 0.5, size = 3, pch = 21,  color = "black")


# get the slope and intercept
lin_fit <- function(dat){
  the_fit <- lm(dat$abundance ~ dat$week, dat)
  p_val <- anova(the_fit)$'Pr(>F)'[1]
  slo_int <- data.frame(t(coef(the_fit)))
  r_sq <- summary(the_fit)$r.squared
  result <- cbind(slo_int, r_sq, p_val)
  colnames(result) <- c("intercept", "slope", "r_squared", "p_value")
  return(result)
}

# do that function to each subgroup
df_fit <- df %>%
          group_by(Master.Protein.Accessions, sex, leg) %>%
          do(lin_fit(.))

# 






