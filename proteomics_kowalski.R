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
max_pep <- 10000

data <-
  tbl_df(data_raw) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = max_pep, wt = ctrl_raw_med)

# # proteins used for normalization
# histones <- "B2RTM0"
# 
# norm_pro <- histones
# 
# norm_pep <- subset(data, data$Master.Protein.Accessions == norm_pro)
# numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
# raw_abun <- numeric_cols[-length(numeric_cols)]
# norm_value <- sapply(norm_pep[, raw_abun], mean)
# 
# raw_abun_mat <- as.matrix(data[, raw_abun])
# 
# norm_abun <- t(t(raw_abun_mat)/norm_value)
# colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
# norm_test <- as.data.frame(norm_abun)
# data <- as_tibble(data)
# data <- cbind(data, norm_test)

# Median, sd, & ratio of peptides ----------------------------------------------

# group1 <- grep("infected_norm", colnames(data))
# group2 <- grep("HMM_norm", colnames(data))
# ctrl <- grep("Control_norm", colnames(data))
# 
# data$group1_med <- by_group(data, group1)
# data$group2_med <- by_group(data, group2)
# data$ctrl_med <- by_group(data, ctrl)

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

# convert protein data.frame to tidy data so we can do the plotting  

samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "sex", "leg", "week"), sep = "_")

# average together the duplicates in the same week
df_avg <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

df_avg$week <- as.numeric(df_avg$week)


cor(df_avg$week, df_avg$abundance)

df_avg %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(corr = cor(df_avg$week, df_avg$abundance))


# need to first make a plot of the whole thing

ggplot(df_avg, aes(x = week, y = abundance)) +
  geom_jitter() +
  geom_smooth(lwd = 3, se = FALSE, method = "lm")

ov_cor <- df_avg %>%
    cor(df_avg$week, df_avg$abundance)

gcor <- df_avg %>%
    group_by(Master.Protein.Accessions, leg) %>%
    summarize(correlation = cor(week, abundance))

lin_fit <- function(dat) {
  the_fit <- lm(dat$abundance ~ I(dat$week), dat)
  setNames(data.frame(t(coef(the_fit))), c("intercept", "slope"))
}



gfits_me <- df_avg %>%
  group_by(Master.Protein.Accessions, leg) %>% 
  do(lin_fit(.))

gfits_me <- df_avg %>%
  group_by(Master.Protein.Accessions, leg) %>%
  do(cor(df_avg$week, df_avg$abundance))

# great site to help http://stat545.com/block023_dplyr-do.html


# spreads out to columns by week
df <- as.data.frame(spread(df_avg, "week", "abundance"))








# plot the 
pro <- as.character(unique(df$Master.Protein.Accessions))

for (i in 1:length(pro)){
  
  temp <- filter(df, Master.Protein.Accessions == pro[i])
  
  pep <- as.character(unique(temp$peptide))
  
  for (pep_x in pep){
    
    temp_pep <- filter(temp, peptide == pep_x)
    
    cols <- c(grep(c("01"), colnames(temp_pep)),
              grep(c("02"), colnames(temp_pep)))
    
    for (j in cols){
      
      g1 <- ggplot(data = temp_pep) +
        geom_point(mapping = aes(x = hrs, 
                                 y = as.numeric(temp_pep[, colnames(temp_pep)[j]]), 
                                 color = Group), 
                   size = 3, alpha = 0.8) + 
        labs(title = paste(pro[i], pep_x, sep = " -- "), 
             subtitle = colnames(temp_pep)[j],
             y = "Abundance", x = "Time (Weeks)") +
        theme_classic()
      plot(g1)
    }
  }
}




