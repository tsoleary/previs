# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)

# source the functions in the proteomics_functions.R script!

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# Tidy the data ----------------------------------------------------------------

df_messy <- read.csv("kowalski_F_week_8_11_Myh_Hist_062819_isotopes.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df_tidy <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sex", "mouse", "leg", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week)

# # average together the mice legs in the same week lines for each peptide
# df <- df_tidy %>%
#   group_by(protein, peptide, isotope, sex, leg, week) %>%
#   summarize(abundance = mean(abundance, na.rm = TRUE))

#     (M_1 - (M_1/M_0(T0)) * M_0))
# -------------------------------------
# (M_0 + (M_1 - ((M_1/M_0(T0)) * M_0)))

# m1_m0_r ratio calculation using messy data

df_iso <- df_tidy %>%
  filter(., isotope == "M_0" | isotope == "M_1") %>%
  group_by(protein, peptide, isotope, sex, leg)

df_messy <- df_iso %>%
  spread(., "isotope", "abundance")

df_messy <- as.data.frame(df_messy)

df_avg <- df_messy %>%
  group_by(protein, peptide, sex, leg, week) %>%
  summarize(., M_0 = mean(M_0, na.rm = TRUE), M_1 = mean(M_1, na.rm = TRUE)) %>%
  mutate(., m1_m0_r = M_1 / M_0)


ratios <- NULL

for (i in 1:nrow(df_avg)){
  
  t0_ratio <- df_avg$m1_m0_r[which(df_avg$protein == df_avg$protein[i] & 
                                   df_avg$peptide == df_avg$peptide[i] &
                                   df_avg$sex == df_avg$sex[i] & 
                                   df_avg$leg == df_avg$leg[i] & 
                                   df_avg$week == 8)]
  
  which(df_avg$protein == df_avg$protein[i] & 
        df_avg$peptide == df_avg$peptide[i] &
        df_avg$sex == df_avg$sex[i] & 
        df_avg$leg == df_avg$leg[i] & 
        df_avg$week == 8)
  
  temp <- (df_avg$M_1[i] - (t0_ratio * df_avg$M_0[i])) / 
    (df_avg$M_0[i] + (df_avg$M_1[i] - (t0_ratio * df_avg$M_0[i])))
  
  ratios <- c(ratios, temp)

  }

df_avg$ratio_cor <- ratios

