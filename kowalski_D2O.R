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


ratios <- NULL

for (i in 1:nrow(df_avg)){
  
  t0_ratio <- df_avg$m1_m0_r[which(df_avg$protein == df_avg$protein[i] & 
                                   df_avg$peptide == df_avg$peptide[i] &
                                   df_avg$sex == df_avg$sex[i] & 
                                   df_avg$leg == df_avg$leg[i] & 
                                   df_avg$week == 0)]

  temp <- (df_avg$M_1[i] - (t0_ratio * df_avg$M_0[i])) / 
    (df_avg$M_0[i] + (df_avg$M_1[i] - (t0_ratio * df_avg$M_0[i])))
  
  ratios <- c(ratios, temp)

  }

df_avg$ratio <- ratios

# plot and run a non-linear regression

# plot_iso_ratios function
plot_iso_r <- function(dat, g_title, FUN = geom_jitter){
  
  g <- ggplot(dat, aes(x = week, y = ratio)) +
    FUN(mapping = aes(x = week, y = ratio, fill = leg), 
        alpha = 0.5, size = 3, pch = 21,  color = "black", width = 0.05) +
    geom_smooth(method = "nls", 
                method.args = list(formula = y ~ yf * x / (tf + x),
                                   start = list(tf = 3, yf = .6)),
                data = df_g,
                se = FALSE,
                aes(color = leg), show.legend = FALSE) +
    labs(title = g_title, x = "Week", y = "M0/(M0+M1)\nCorrected", 
         fill = "Leg") +
    expand_limits(x = 0, y = 0) +
    theme_classic()
  
  return(g)
  
}


ggplot(df_g)


df_g <- filter(df_avg, protein == "Myh4" & peptide == "LQDAEEHVEAVNSK")


ggplot(df_g, aes(week, ratio)) +
  geom_point() +   
  geom_smooth(method = "nls", 
              method.args = list(formula = y ~ yf * x / (tf + x),
                                 start = list(tf = 3, yf = .6)),
              data = df_g,
              se = FALSE,
              aes(color = leg))



plot_iso_r(df_g, "WORK!")



fit <- nls(ratio_cor ~ SSasymp(week, yf, y0, log_alpha), data = df_g)


df_g %>% 
  group_by(protein, peptide, leg) %>% 
  do(fit = nls(ratio_cor ~ SSasymp(week, yf, y0, log_alpha), data = .)) %>% 
  tidy(fit) %>% 
  select(sensor, term, estimate) %>% 
  spread(term, estimate) %>% 
  mutate(alpha = exp(log_alpha))
