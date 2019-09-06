# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(broom)
library(RColorBrewer)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Kowalski")

# Tidy the data ----------------------------------------------------------------

df_raw <- read.csv("Kowalski_M_week_8_11 Super library with Myh and Cox.csv")

samples <- colnames(df_raw)[4:length(colnames(df_raw))]

df_tidy <- df_raw %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sex", "mouse", "leg", "week"), sep = "_")

df_tidy$week <- as.numeric(df_tidy$week) - 8

# # average together the mice legs in the same week lines for each peptide
# df_tidy <- df_tidy %>%
#   group_by(protein, peptide, isotope, sex, leg, week) %>%
#   summarize(abundance = mean(abundance, na.rm = TRUE))

# # As above, but with Peptides averaged into protein
df_tidy <- df_tidy %>%
  group_by(protein, isotope, sex, leg, week) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

df <- df_tidy %>%
  spread(., "isotope", "abundance") %>%
  as.data.frame(.) %>%
  norm_iso(.)

df <- df[which(!is.na(df$M_3_norm)), ]

df$ratio <- iso_sum_ratio(df, "Myh1")

# spread out by isotope and add normalization to sum columns
# df <- df_tidy %>%
#   spread(., "isotope", "abundance") %>%
#   as.data.frame(.) %>%
#   norm_iso(.)

# ratio calculation ------------------------------------------------------------

# Fractional M+0 abundance: M_0/(M_1+M_2+M_3)

# corrected ratio
#     (M_1 - (M_1(T0)/M_0(T0)) * M_0))
# -------------------------------------
# (M_0 + (M_1 - ((M_1(T0)/M_0(T0)) * M_0)))

# more corrected ratio
#     (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)(T0)/M_0(T0)) * M_0))
# -------------------------------------------------------------
# (M_0 + (M_1 + M_2 + M_3 - ((M_1 + M_2 + M_3)(T0)/M_0(T0)) * M_0)))


# # t0_r calculation (M+1 only)
# df <- df[which(!is.na(df$M_3_norm)), ] %>%
#   group_by(protein, peptide, sex, leg, week) %>%
#   mutate(., t0_r = (M_1_norm) / M_0_norm)
# 

# t0_r calculation (M+1, 2, 3)
df <- df[which(!is.na(df$M_3_norm)), ] %>%
  group_by(protein, peptide, sex, leg, week) %>%
  mutate(., t0_r = (M_1_norm + M_2_norm + M_3_norm) / M_0_norm)

## If curve is desired for peptides lacking M+3 peak
# df1 <- df[which(is.na(df$M_3_norm)), ] %>%
#   group_by(protein, peptide, sex, leg, week) %>%
#   mutate(., t0_r = (M_1_norm + M_2_norm) / M_0_norm)
# df1$peptide <- str_pad(df1$peptide, side = "both", width = (max(nchar(as.character(df1$peptide))) + 2), pad = "*")
# 
# df2 <- df[which(!is.na(df$M_3_norm)), ] %>%
#   group_by(protein, peptide, sex, leg, week) %>%
#   mutate(., t0_r = (M_1_norm + M_2_norm + M_3_norm) / M_0_norm)
# 
# df <- union(df1, df2)
# rm(df1, df2)
# corrected ratio calculation

# df$ratio <- iso_ratio_calc_one(df) #For M+1 only
# 
# df$ratio <- iso_ratio_calc(df) #For M+1, 2, 3

# # modified proteomixr function! Use if you *really* want to see corrected ratios for
# # peptides without M+3 data.
# iso_ratio_calc_more <- function(dat){
#   
#   ratios <- NULL
#   
#   for (i in 1:nrow(dat)){
#     
#     t0_ratio <- mean(dat$t0_r[which(dat$protein == dat$protein[i] & 
#                                       dat$peptide == dat$peptide[i] &
#                                       dat$sex == dat$sex[i] & 
#                                       dat$leg == dat$leg[i] & 
#                                       dat$week == 0)], na.rm = TRUE)
#     if(is.na(dat$M_3[i]) == FALSE) {
#       temp <- (dat$M_1[i] + dat$M_2[i] + dat$M_3[i] - (t0_ratio * dat$M_0[i])) / 
#         (dat$M_0[i] + (dat$M_1[i] + dat$M_2[i] + dat$M_3[i] - 
#                          (t0_ratio * dat$M_0[i])))
#     }
#     if(is.na(dat$M_3[i]) == TRUE) {
#       temp <- (dat$M_1[i] + dat$M_2[i] - (t0_ratio * dat$M_0[i])) / 
#         (dat$M_0[i] + (dat$M_1[i] + dat$M_2[i] - 
#                          (t0_ratio * dat$M_0[i])))
#     }
#     ratios <- c(ratios, temp)
#     
#   }
#   
#   return(ratios)
#   
# }
# 
# df$ratio <- iso_ratio_calc_more(df)
# df <- filter(df, is.na(ratio) != TRUE)

# df <- filter(df, ratio >= 0)

## Get sum ratio over time for given, time-, leg-, and sex-matched peptide divisor

df$ratio <- iso_sum_ratio(df, "NAYEESLDHLETLKR")


# plot and run a non-linear regression for the corrected isotope ratio ---------

# filter for learning!
df_g <- filter(df, protein == "Act" & peptide == "QEYDEAGPSIVHR")
# plot_iso_r(df_g, g_title = "Protein", g_subtitle = "PEPTIDE")
plot_m0(df_g, g_title = "Protein", g_subtitle = "PEPTIDE")
plot_sum_ratio(df_g, "PEPTIDE")

# # plot all proteins and peptides
# pros <- as.character(unique(df$protein))
# 
# plot_list <- list()
# 
# for (pro in pros){
#   pro_df <- filter(df, protein == pro)
#   peps <- as.character(unique(pro_df$peptide))
#     for (pep in peps){
#       pep_df <- filter(pro_df, peptide == pep)
#       grapherror <- tryCatch({
#       g <- plot_iso_r(pep_df, g_title = pro, g_subtitle = pep)
#       }, error = function(error){ print(error) })
# 
#       if(inherits(grapherror, "error")) {next}
#       plot_list[[pep]] <- g
#     }
# }
# 
# pdf("Male Isotopic Ratio Weeks 8-11 Individual 4.pdf", width = 10.75, height = 6)
# 
# for(i in 1:length(plot_list)){
#   print(plot_list[[i]])
# }
# 
# dev.off() # pdf file will appear in working directory

pros <- as.character(unique(df$protein))

plot_list <- list()

for (pro in pros){
  pro_df <- filter(df, protein == pro)
  grapherror <- tryCatch({
    g <- plot_sum_ratio(pro_df, g_title = pro)
  }, error = function(error){ print(error) })
    if(inherits(grapherror, "error")) {next}
    plot_list[[pep]] <- g
  }

pdf("test3 myh1.pdf", width = 10.75, height = 6)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory


# pros <- as.character(unique(df$protein))
# 
# plot_list <- list()
# 
# for (pro in pros){
#   pro_df <- filter(df, protein == pro)
#   peps <- as.character(unique(pro_df$peptide))
#   for (pep in peps){
#     pep_df <- filter(pro_df, peptide == pep)
#     grapherror <- tryCatch({
#       g <- plot_m0(pep_df, g_title = pro, g_subtitle = pep)
#     }, error = function(error){ print(error) })
#     
#     if(inherits(grapherror, "error")) {next}
#     plot_list[[pep]] <- g
#   }
# }
# 
# pdf("Male M0 Abundance Weeks 8-11 Individual Expanded.pdf", width = 10.75, height = 6)
# 
# for(i in 1:length(plot_list)){
#   print(plot_list[[i]])
# }
# 
# dev.off() # pdf file will appear in working directory

# plot the changing isotopic distribution over time ----------------------------

# tidy the data frame back by isotope
df_iso <- df %>%
  gather(., isotope, abundance, colnames(.)[grep("_norm", colnames(.))]) 

df_iso$isotope <- gsub("_norm", "", df_iso$isotope)
df_iso <- df_iso[!is.nan(df_iso$abundance), ]
df_iso$week <- as.character(df_iso$week)

#filter for learning
df_g <- filter(df_iso, protein == "Act", peptide == "QEYDEAGPSIVHR")
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

pdf("Male Isotopic Distribution Weeks 8-11 Individual Expanded.pdf", width = 10.75, height = 5)

for(i in 1:length(plot_list)){
  print(plot_list[[i]])
}

dev.off() # pdf file will appear in working directory

