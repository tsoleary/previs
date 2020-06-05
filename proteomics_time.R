# proteins changing over time --------------------------------------------------

library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# read in the data.frame with the top 3 avg for each protein & each sample
# protein <- read.csv("kowalski_M_week_1_8_top_3_abund_norm_sum_total_r.csv")
# protein <- inner_join(read.csv("kowalski_F_w1_w8_norm_sum_total_top_3_r.csv"),
#                       read.csv("kowalski_M_week_1_8_top_3_abund_norm_sum_total_r.csv"),
#                       by = c("Master.Protein.Accessions", "gene"))

protein <- full_join(read.csv("kowalski_F_w1_w8_norm_sum_total_top_3_r.csv"),
                      read.csv("kowalski_M_week_1_8_top_3_abund_norm_sum_total_r.csv"),
                      by = "Master.Protein.Accessions")

protein$X.x <- NULL 
protein$X.y <- NULL 
protein <- rename(protein, "peptides_fem" = "peptides.x", "peptides_male" = "peptides.y")
protein$gene.y <- NULL

colnames(protein) <- gsub(".x", "", colnames(protein))

# colnames should only have File_Leg_Sex_Week
colnames(protein) <- gsub("_norm", "", colnames(protein))
colnames(protein) <- gsub("_Sample", "", colnames(protein))

# convert protein data.frame to a tidy data.frame
samples <- colnames(protein)[grep("F", colnames(protein))]

df_tidy <- protein %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("file", "leg", "sex", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week)

# New column for mouse number

mouse_dat <- read.csv("Mouse file list.csv")
mouse_dat$file <- as.character(mouse_dat$file) %>%
  str_pad(., 4, "right", pad = "-")

df_tidy$file <- as.character(df_tidy$file) %>%
  str_pad(., 4, "right", pad = "-")

df_tidy$individual <- file_to_ind(df_tidy, mouse_dat)

# median together the duplicates in the same week (STILL RUN FIRST TWO LINES
# EVEN IF MEDIAN IS DESIRED - following code requires data be
# object named 'df')

df <- df_tidy %>%
  group_by(Master.Protein.Accessions, leg, week) %>%
  summarize(abundance = median(abundance, na.rm = TRUE))

# remove the NA's & proteins that have data in only one week or one leg
df <- df[!is.na(df$abundance), ]
df <- df[df$Master.Protein.Accessions != "", ]

df <- df %>%
group_by(Master.Protein.Accessions, leg) %>%
do(filter(.,length(unique(week)) > 1))

# If not grouping by individual mouse number
# df <- df %>%
# group_by(Master.Protein.Accessions, sex) %>%
# do(filter(.,length(unique(leg)) > 1))

# If data from both legs is desired per-protein in each individual mouse 
# (Removes if data from both legs not found)
df <- df %>%
  group_by(Master.Protein.Accessions, individual) %>%
  do(filter(., length(unique(leg)) > 1))

# Sort by Mouse identifier

df$individual <- as.numeric(df$individual)
df <- arrange(df, individual)

# Export merged 'raw'

write.csv(df, file = "Protein values norm sum total from top 3 pruned.csv", row.names = FALSE)


## Stats! ----------------------------------------------------------------------

# Paired ttest per-protein, per-week, comparing legs

pval_df <- paired_leg_ttest(df)
gene_df <- read.csv('Kowalski_F_w1_w8_gene_names.csv')
pval_df$gene <- mpa_to_gene(pval_df, gene_df)
write.csv(pval_df, "Kowalski Paired Ttest Pvalues.csv", row.names = FALSE)

# # Two-way, repeated-measures ANOVA of protein abundance values by individual------
# # and leg
# 
# # prot_aov <- with(df[which(df$Master.Protein.Accessions == "A0A075DC90"), ],
# #                    aov(abundance ~ leg * week +
# #                          Error(individual / (leg * week)))
# # )
# # 
# # summary(prot_aov)
# 
# # week_as_ind <- with(df[which(df$Master.Protein.Accessions == "A0A075DC90"), ],
# #                  aov(abundance ~ leg * individual +
# #                        Error(week / (leg * individual)))
# # )
# # 
# # summary(week_as_ind)

## Per-protein average right-leg to left-leg ratio of abundance-----------------

# df$ratio <- protein_leg_ratio(df)
df$ratio <- protein_leg_ratio_ind(df)

# spread out df to make an excel file of the numbers with the data that we graph

# df_L_R <- df %>%
#   spread(., "leg", "abundance")
# 
# library(data.table)
# 
# df_csv <- dcast(setDT(df_L_R), Master.Protein.Accessions ~ week, 
#                  value.var = c("L", "R")) 
# 
# gene_df <- read.csv("Kowalski_M_week_1_8_gene_names.csv")
# 
# df_csv$gene <- mpa_to_gene(df_csv, gene_df)
# 
# write.csv(df_csv, "kowalski_M_w1_w8_med_protein_time_r.csv")

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
df_g <- filter(df, Master.Protein.Accessions == "A0A075DC90")

# for test the plots and learning :)
plot_pro(df_g, "PLOTTING IS FUN")
plot_rat(df_g, "PLOTTING IS FUN")

# plot all proteins OR R/L ratio ------------------------------------------------
pros <- as.character(unique(df$Master.Protein.Accessions))

plot_list <- list()

for (pro in pros){
  gene <- indiv_mpa_to_gene(pro, protein)
  temp_df <- filter(df, Master.Protein.Accessions == pro)
  # g <- plot_pro(temp_df, g_title = gene)
  g <- plot_rat(temp_df, g_title = gene)
  plot_list[[pro]] <- g
}

pdf("Pool Week 1-8 Protein Abundance Ratios Individual per Individual.pdf", width = 10.75, height = 6)

for(pro in pros){
  print(plot_list[[pro]])
}

dev.off() # pdf file will appear in working directory
