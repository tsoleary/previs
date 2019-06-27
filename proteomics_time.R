# proteins changing over time --------------------------------------------------

library(tidyverse)

# source the functions in the proteomics_functions.R script

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# read in the data.frame with the top 3 avg for each protein & each sample
protein <- read.csv("kowalski_F_w1_w8_no_norm_r.csv")
protein$X <- NULL 

# colnames should only have File_Leg_Sex_Week
colnames(protein) <- gsub("_norm", "", colnames(protein))
colnames(protein) <- gsub("_Sample", "", colnames(protein))

# convert protein data.frame to a tidy data.frame
samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "leg", "sex", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week)

# median together the duplicates in the same week
df <- df_tidy %>%
  group_by(Master.Protein.Accessions, sex, leg, week) %>%
  summarize(abundance = median(abundance, na.rm = TRUE))

# remove the NA's & proteins that have data in only one week or one leg
df <- df[!is.na(df$abundance), ]

df <- df %>%
  group_by(Master.Protein.Accessions, sex, leg) %>%
  do(filter(., length(unique(week)) > 1)) 

df <- df %>%
  group_by(Master.Protein.Accessions, sex) %>%
  do(filter(., length(unique(leg)) > 1))

# for linear fit data frame ----------------------------------------------------
# # need to first make a plot of the whole thing
# ggplot(df, aes(x = week, y = abundance)) +
#   geom_jitter(mapping = aes(x = week, y = abundance, fill = leg), 
#               alpha = 0.5, size = 3, pch = 21,  color = "black")
# 
# # get get the slope, intercept r_sq, and pval of each group
# df_fit <- df %>%
#   group_by(Master.Protein.Accessions, sex, leg) %>%
#   do(lin_fit(.))
# 
# df_fit$gene <- mpa_to_gene(df_fit, gene_df)

# plotting check on one protein ------------------------------------------------
# filter for learning
df_g <- filter(df, Master.Protein.Accessions == "A0A068BFR3")

# for test the plots and learning :)
plot_pro(df_g, "PLOTTING IS FUN", FUN = geom_jitter)

# plot all proteins ------------------------------------------------------------
pros <- as.character(unique(df$Master.Protein.Accessions))

plot_list <- list()

for (pro in pros){
  gene <- indiv_mpa_to_gene(pro, protein)
  temp_df <- filter(df, Master.Protein.Accessions == pro)
  g <- plot_pro(temp_df, g_title = gene)
  plot_list[[pro]] <- g
}

pdf("plot_F_w1_w8_norm_sum_total_top_3_all.pdf", width = 10.75, height = 6)

for(pro in pros){
  print(plot_list[[pro]])
}

dev.off() # pdf file will appear in working directory
