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

# throw gene names into df_fit
df_fit$gene <- mpa_to_gene(df_fit, gene_df)

# filter 
df_mybpc1 <- filter(df, gene == "Mybpc1")

df_mybpc1_g <- filter(df, Master.Protein.Accessions == "G0YZM8")

df_mybpc1_q <- filter(df, Master.Protein.Accessions == "Q6P6L5")

test <- filter(df, Master.Protein.Accessions == "Q3UIJ3; P68134")

# need to first make a plot of the whole thing
ggplot(test, aes(x = week, y = abundance)) +
  geom_point(mapping = aes(x = week, y = abundance, fill = leg), 
              alpha = 0.5, size = 3, pch = 21,  color = "black") +
  labs(title = "Q3UIJ3; P68134", x = "Week", y = "Raw Abundance", fill = "Leg") +
  theme_classic()
  
# make this a function

plot_pro <- function(dat, tle = pro){
  g <- ggplot(dat, aes(x = week, y = abundance)) +
         geom_point(mapping = aes(x = week, y = abundance, fill = leg), 
                    alpha = 0.5, size = 3, pch = 21,  color = "black") +
         labs(title = title, x = "Week", y = "Raw Abundance", fill = "Leg") +
         theme_classic()
  return(g)
}

plot_pro(Myh2, title = pro)



# make a loop to make multiple plots

pros <- c("Q3UIJ3; P60710; P68134", "Q3UIJ3; P68134", "P07310", 
          "Q5SX40; Q5SX39; G3UW82", "Q5SX40; P13541; Q5SX39; B1AR69", 
          "P05977", "Q5SX39", "Q5SX40; P13541; Q5SX39; G3UW82", "B1AR69", 
          "P68134", "Q5SX40; P13541; Q5SX39; G3UW82; B1AR69", "P97457", 
          "P58774; Q545Y3", "A6ZI44", "A0A0A0MQF6", "Q5SX40; Q5SX39", 
          "Q5SX39; G3UW82", "Q5SX40; P13541; Q5SX39", 
          "Q5SX40; Q5SX39; G3UW82; B1AR69", "Q8R429", "P21550", 
          "E9Q8K5; A2ASS6", "Q7TPR4", "Q9JI91", "Q9JI91; Q7TPR4", 
          "O88990", "O88990; Q7TPR4", "O88990; Q9JI91", 
          "O88990; Q9JI91; Q7TPR4")




x <- df$Master.Protein.Accessions == "neil"

sum(x)

list_nots <- NULL

for (pro in pros){
  
  temp <- df$Master.Protein.Accessions == pro
  sum_temp <- sum(temp)
  if(sum_temp == 0){
    list_nots <- c(list_nots, pro)
    
  }
}


list_ares <- NULL

for (pro in pros){
  
  temp <- df$Master.Protein.Accessions == pro
  sum_temp <- sum(temp)
  if(sum_temp != 0){
    list_ares <- c(list_ares, pro)
    
  }
}


pros <-  c("P07310", "Q5SX40; Q5SX39; G3UW82", "P05977", "Q5SX39", "P68134", 
           "P97457", "A6ZI44", "A0A0A0MQF6", "Q5SX40; Q5SX39", "Q5SX39; G3UW82", 
           "Q8R429", "P21550", "E9Q8K5; A2ASS6", "Q7TPR4", "Q9JI91", "O88990", 
           "O88990; Q9JI91", "O88990; Q9JI91; Q7TPR4")






for (pro in pros){
  temp_df <- filter(df, Master.Protein.Accessions == pro)
  
  g <- plot_pro(temp_df)
  
  plot(g)
}


length(unique(pros))

