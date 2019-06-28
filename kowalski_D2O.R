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

df_T0 <- df_tidy %>%
  filter(., week == 8.0 & isotope == "M_0" | week == 8 & isotope == "M_1") %>%
  group_by(protein, peptide, isotope, sex, leg)

df_messy <- df_T0 %>%
  spread(., "isotope", "abundance") %>%
  mutate(., m1_m0_r = M_1 / M_0)

# how do we do this in r in a tidy_df?? !!!!!!!!!

iso_ratio <- function(dat){
  result <- dat$abundance[which(dat$isotope == "M_1")] / 
    dat$abundance[which(dat$isotope == "M_0")]
  result <- as.data.frame(result)
  return(result)
}

df_T0 <- df_tidy %>%
  filter(., week == 8.0 & isotope == "M_0" | week == 8 & isotope == "M_1") %>%
  group_by(protein, peptide, isotope, sex, leg) %>%
  do(iso_ratio(.))







