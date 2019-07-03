# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(broom)
library(RColorBrewer)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# Tidy the data ----------------------------------------------------------------

df_raw <- read.csv("kowalski_F_week_8_11_Myh_Hist_062819_isotopes.csv")

samples <- colnames(df_raw)[4:length(colnames(df_raw))]

df_tidy <- df_raw %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sex", "mouse", "leg", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week) - 8

# average together the mice legs in the same week lines for each peptide
df_tidy <- df_tidy %>%
  group_by(protein, peptide, isotope, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

# spread out by isotope and add normalization to sum columns
df <- df_tidy %>%
  spread(., "isotope", "abundance") %>%
  as.data.frame(.) %>%
  norm_iso(.)

# ratio calculation ------------------------------------------------------------

# corrected ratio
#     (M_1 - (M_1/M_0(T0)) * M_0))
# -------------------------------------
# (M_0 + (M_1 - ((M_1/M_0(T0)) * M_0)))

# more corrected ratio
#     (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)/M_0(T0)) * M_0))
# -------------------------------------------------------------
# (M_0 + (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)/M_0(T0)) * M_0)))

# t0_r calculation
df <- df %>%
  group_by(protein, peptide, sex, leg, week) %>%
  mutate(., t0_r = (M_1_norm + M_2_norm + M_3_norm) / M_0_norm)

# corrected ratio calculation
df$ratio <- iso_ratio_calc(df)
df <- filter(df, is.na(ratio) != TRUE)

# plot and run a non-linear regression for the corrected isotope ratio ---------

# filter for learning!
df_g <- filter(df, protein == "Shared Myosin" & peptide == "DTQLHLDDALR")
plot_iso_r(df_g, g_title = "Protein", g_subtitle = "PEPTIDE")

# plot all proteins and peptides
pros <- as.character(unique(df$protein))

plot_list <- list()

for (pro in pros){
  pro_df <- filter(df, protein == pro)
  peps <- as.character(unique(pro_df$peptide))
    for (pep in peps){
      pep_df <- filter(pro_df, peptide == pep)
      tryCatch({
        g <- plot_iso_r(pep_df, g_title = pro, g_subtitle = pep)
      }, error = function(error){})
      plot_list[[pep]] <- g
    }
}

pdf("plot_F_w8_w_11_syn_deg_curve_all.pdf", width = 10.75, height = 6)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory

# plot the changing isotopic distribution over time ----------------------------

# tidy the data frame back by isotope
df_iso <- df %>%
  gather(., isotope, abundance, colnames(.)[grep("_norm", colnames(.))]) 

df_iso$isotope <- gsub("_norm", "", df_iso$isotope)
df_iso <- df_iso[!is.nan(df_iso$abundance), ]
df_iso$week <- as.character(df_iso$week)

#filter for learning
df_g <- filter(df_iso, protein == "Hist2h4", peptide == "DAVTYTEHAK")
plot_all_isos(df_g, "PRO", "PEP")

# plot the isotope distribution for each peptide changing over time
pros <- as.character(unique(df$protein))
plot_list <- list()

for (pro in pros){
  pro_df <- filter(df_iso, protein == pro)
  peps <- as.character(unique(pro_df$peptide))
  for (pep in peps){
    pep_df <- filter(pro_df, peptide == pep)
    g <- plot_all_isos(pep_df, g_title = pro, g_subtitle = pep)
    plot_list[[pep]] <- g
  }
}

pdf("plot_F_w8_w_11_isotopomer.pdf", width = 10.75, height = 5)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory

