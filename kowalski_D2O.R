# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(broom)

# source the functions in the proteomics_functions.R script!

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# Tidy the data ----------------------------------------------------------------

df_messy <- read.csv("kowalski_F_week_8_11_Myh_Hist_062819_isotopes.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df_tidy <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sex", "mouse", "leg", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week) - 8

# # average together the mice legs in the same week lines for each peptide
# df <- df_tidy %>%
#   group_by(protein, peptide, isotope, sex, leg, week) %>%
#   summarize(abundance = mean(abundance, na.rm = TRUE))

#     (M_1 - (M_1/M_0(T0)) * M_0))
# -------------------------------------
# (M_0 + (M_1 - ((M_1/M_0(T0)) * M_0)))


# more corrected
#     (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)/M_0(T0)) * M_0))
# -------------------------------------
# (M_0 + (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)/M_0(T0)) * M_0)))

# m1_m0_r ratio calculation using messy data

df_iso <- df_tidy %>%
  filter(., isotope == "M_0" | isotope == "M_1" | isotope == "M_2" | 
           isotope == "M_3") %>%
  group_by(protein, peptide, isotope, sex, leg)

df_messy <- df_iso %>%
  spread(., "isotope", "abundance") %>%
  as.data.frame(.)

df_avg <- df_messy %>%
  group_by(protein, peptide, sex, leg, week) %>%
  summarize(., M_0 = mean(M_0, na.rm = TRUE), M_1 = mean(M_1, na.rm = TRUE),
            M_2 = mean(M_2, na.rm = TRUE), M_3 = mean(M_3, na.rm = TRUE)) %>%
  mutate(., t0_r = (M_1 + M_2 + M_3) / M_0)

df_messy <- df_messy %>%
  mutate(., t0_r = (M_1 + M_2 + M_3) / M_0)

df_avg$ratio <- iso_ratio_calc(df_avg)

df_messy$ratio <- iso_ratio_calc(df_messy)

df_avg <- filter(df_avg, is.na(ratio) != TRUE)

df_messy <- filter(df_messy, is.na(ratio) != TRUE)

# plot and run a non-linear regression -----------------------------------------

# filter for learning!
df_g <- filter(df_messy, protein == "Shared Myosin" & peptide == "DTQLHLDDALR")

plot_iso_r(df_g, g_title = "Protein", g_subtitle = "PEPTIDE")

# plot all proteins and peptides

pros <- as.character(unique(df_messy$protein))

plot_list <- list()

for (pro in pros){
  pro_df <- filter(df_messy, protein == pro)
  peps <- as.character(unique(pro_df$peptide))
    for (pep in peps){
      pep_df <- filter(pro_df, peptide == pep)
      tryCatch({
        g <- plot_iso_r(pep_df, g_title = pro, g_subtitle = pep)
      }, error = function(error){})
      plot_list[[pep]] <- g
    }
}

pdf("test.pdf", width = 10.75, height = 6)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory
