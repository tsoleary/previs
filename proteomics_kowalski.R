# Proteomic analysis: Data imported from Proteome Discoverer 2.2 ---------------

library("tidyverse")

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
norm_value <- sapply(norm_pep[, raw_abun], mean, na.rm = TRUE)

raw_abun_mat <- as.matrix(data[, raw_abun])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data <- as_tibble(data)
data <- cbind(data, norm_test)


# sum of every column 

norm_value <- colSums(data_raw[, ctrl_raw], na.rm = TRUE)

raw_abun_mat <- as.matrix(data_raw[, ctrl_raw])

norm_abun <- t(t(raw_abun_mat)/norm_value)
colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
norm_test <- as.data.frame(norm_abun)
data_raw <- as_tibble(data_raw)
data <- cbind(data_raw, norm_test)


# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(data) %>%
  group_by(Master.Protein.Accessions) %>%
  top_n(n = pep_top, wt = ctrl_raw_med) %>%
  as.data.frame

# Protein Averages -------------------------------------------------------------

group_names <- colnames(data_top)[grep("norm", colnames(data))]

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

write.csv(protein, "kowalski_F_w1_w8_norm_to_sum_total_r.csv")




################################################################################
# proteins changing over time --------------------------------------------------

# protein <- read.csv("kowalski_F_w1_w8_norm_to_sum_total_r.csv")
# protein$X <- NULL

colnames(protein) <- gsub("_norm", "", colnames(protein))
colnames(protein) <- gsub("_Sample", "", colnames(protein))

# convert protein data.frame to a tidy data.frame
samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "leg", "sex", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week)

#remove the NA's
df_tidy <- df_tidy[!is.na(df$abundance), ]

# mean together the duplicates in the same week
df <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

# median together the duplicates in the same week
df_med <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = median(abundance, na.rm = TRUE))

df$gene <- mpa_to_gene(df, gene_df)

#remove the NA's
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

df_fit$gene <- mpa_to_gene(df_fit, gene_df)

# filter for learning
df_g <- filter(df, Master.Protein.Accessions == "G0YZM8")

# PLOT
ggplot(df_g, aes(x = week, y = abundance)) +
  geom_point(mapping = aes(x = week, y = abundance, fill = leg), 
             alpha = 0.5, size = 3, pch = 21,  color = "black") +
  labs(title = "Neil", x = "Week", y = "Raw Abundance", fill = "Leg") +
  theme_classic() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(mapping = aes(color = leg), method = 'lm', se = FALSE, 
              size = 1.1, show.legend = FALSE, linetype = "dotted")

# function to get the regression eqn
lm_eqn <- function(lm_object) {
  eq <-
    substitute(
      italic(y) == a + b %.% italic(x) * "," ~ ~ italic(r) ^ 2 ~ "=" ~ r2,
      list(
        a = as.character(signif(coef(lm_object)[1], digits = 2)),
        b = as.character(signif(coef(lm_object)[2], digits = 2)),
        r2 = as.character(signif(summary(lm_object)$r.squared, digits = 3)),
        p = as.character(signif(p_value, digits = 2))
      )
    )
  as.character(as.expression(eq))
}

# plot_pro function
plot_pro <- function(dat, g_title, FUN = geom_point){
  
  df_l <- filter(dat, leg == "L")
  df_r <- filter(dat, leg == "R")
  
  lm_l <- lm(abundance ~ week, df_l)
  lm_r <- lm(abundance ~ week, df_r)
  
  eqn_l <- lm_eqn(lm_l)
  eqn_r <- lm_eqn(lm_r)
  
  g <- ggplot(dat, aes(x = week, y = abundance)) +
    FUN(mapping = aes(x = week, y = abundance, fill = leg), 
                alpha = 0.5, size = 3, pch = 21,  color = "black", width = 0.05) +
    labs(title = g_title, x = "Week", y = "Raw Abundance", fill = "Leg") +
    expand_limits(x = 0, y = 0) +
    theme_classic() +  
    expand_limits(x = 0, y = 0) +
    geom_smooth(mapping = aes(color = leg), method = 'lm', se = FALSE, 
                size = 1.1, show.legend = FALSE, linetype = "dotted") +
    annotate("text", x = 10, y = 0.75*(max(dat$abundance)), 
             label = eqn_l, parse = TRUE, color = "#F98B86") +
    annotate("text", x = 10, y = 0.25*(max(dat$abundance)), 
             label = eqn_r, parse = TRUE, color = "#53D3D7") +
    coord_cartesian(xlim = c(0, 8), clip = 'off') +
    theme(plot.margin = unit(c(1, 8, 1, 1), "lines"))
  return(g)
}

plot_pro(df_g, "???", FUN = geom_jitter)


# make a loop to make multiple plots
pros <-  c("P07310", "Q5SX40; Q5SX39; G3UW82", "P05977", "Q5SX39", "P68134", 
           "P97457", "A6ZI44", "A0A0A0MQF6", "Q5SX40; Q5SX39", "Q5SX39; G3UW82", 
           "Q8R429", "P21550", "E9Q8K5; A2ASS6", "Q7TPR4", "Q9JI91", "O88990", 
           "O88990; Q9JI91", "O88990; Q9JI91; Q7TPR4")


# convert accession to gene for each graph -------------------------------------
indiv_mpa_to_gene <- function (Acc_pro, gene_dat){
  temp <- which(Acc_pro == gene_dat$Master.Protein.Accessions, TRUE)
  gene <- gsub(Acc_pro, gene_dat$gene[temp], Acc_pro)
  return(gene)
}

for (pro in pros){
  
  gene <- indiv_mpa_to_gene(pro, protein)
  
  temp_df <- filter(df, Master.Protein.Accessions == pro)
  
  g <- plot_pro(temp_df, g_title = gene, FUN = geom_jitter)
  
  plot(g)
  
}

