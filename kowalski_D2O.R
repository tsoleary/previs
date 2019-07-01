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

df_messy <- df_messy %>%
  mutate(., m1_m0_r = M_1 / M_0)

df_avg$ratio <- iso_ratio_calc(df_avg)

df_messy$ratio <- iso_ratio_calc(df_messy)

# plot and run a non-linear regression -----------------------------------------

# filter for learning!
df_g <- filter(df_messy, protein == "Shared Myosin" & peptide == "DTQLHLDDALR")

plot_iso_r(df_g, g_title = "Protein", g_subtitle = "PEPTIDE")


# # fit out put from http://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
# df_g %>% 
#   group_by(protein, peptide, leg) %>% 
#   do(fit = nls(ratio ~ SSasymp(week, yf, y0, log_alpha), data = .)) %>% 
#   tidy(fit) %>% 
#   select(leg, term, estimate) %>% 
#   spread(term, estimate) %>% 
#   mutate(alpha = exp(log_alpha))

# plot all proteins and peptides

pros <- as.character(unique(df_messy$protein))

plot_list <- list()

for (pro in pros){
  pro_df <- filter(df_messy, protein == pro)
  peps <- as.character(unique(pro_df$peptide))
    for (pep in peps){
      temp_pep_df <- filter(pro_df, peptide == pep)
      g <- plot_iso_r(temp_pep_df, g_title = pro, g_subtitle = pep)
      plot_list[[pep]] <- g
    }
}

pdf("plot_F_w8_w11_iso_ratios_all.pdf", width = 10.75, height = 6)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory

