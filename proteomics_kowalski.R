# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library(tidyverse)

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")
data_raw <- 
  read.csv("Kowalski_reprex_all_peptides.csv")

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

# set the max number of peptides used in analysis
max_pep <- 100000

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_raw_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- colnames(data)[grep("F", colnames(data))]

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
gene_df <- read.csv('Kowalski_F_3_LandR_w1_trial_gene_names.csv')

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

# proteins changing over time --------------------------------------------------

# convert protein data.frame to a tidy data.frame
samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "sex", "leg", "week"), sep = "_")

# average together the duplicates in the same week
df <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

df$week <- as.numeric(df$week)

df$group <- paste(df$Master.Protein.Accessions, df$sex, df$leg, sep = "_")


# need to first make a plot of the whole thing
ggplot(df, aes(x = week, y = abundance, group = group, color = leg)) +
  geom_point(color = data$Master.Protein.Accessions) +
  geom_smooth(lwd = 1, se = FALSE, method = "lm")


# get the slope and intercept
lin_fit <- function(dat) {
  the_fit <- lm(dat$abundance ~ dat$week, dat)
  setNames(data.frame(t(coef(the_fit))), c("intercept", "slope"))
}

# get the slope and intercept and pvalue
lin_fit2 <- function(dat) {
  fit <- lm(dat$abundance ~ dat$week, dat)
  result <- c(fit$coefficients[1], fit$coefficients[2], anova(fit)$'Pr(>F)'[1])
  setNames(result, c("intercept", "slope", "p_value"))
}

# apply that function to each subgroup
df_s_i2 <- df %>%
          group_by(Master.Protein.Accessions, sex, leg) %>%
          do(lin_fit2(.))


# rsquared value
df_r2 <- df %>%
          group_by(Master.Protein.Accessions, sex, leg) %>%
          summarize(correlation = cor(week, abundance) ^ 2)

df_linear <- full_join(df_s_i, df_r2, by = c("Master.Protein.Accessions" = 
                       "Master.Protein.Accessions", 
                       "sex" = "sex", "leg" = "leg"))

df_linear$group <- paste(df_linear$Master.Protein.Accessions, 
                         df_linear$sex, df_linear$leg, sep = "_")

ggplot()





# great site to help http://stat545.com/block023_dplyr-do.html





