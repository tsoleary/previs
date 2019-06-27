# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(plotly)
library(Rdisop)

# source the functions in the proteomics_functions.R script!

setwd("C:/Users/PrevBeast/Documents/R/Helms")

# Tidy the data ----------------------------------------------------------------

df_messy <- read.csv("All_isotopes_organized_select_reprex.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df_tidy <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("group", "duplicate", "hrs"), sep = "_")

df_tidy$hrs <- as.numeric(df_tidy$hrs)

# average together the duplicate lines
df_avg_3 <- df_tidy %>%
  group_by(protein, peptide, isotope, group, hrs) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

# spreads out to columns by isotope: M_0, M_3, & Sum
df <- as.data.frame(spread(df_avg, "isotope", "abundance"))

# plot by protein, peptide, and isotope ----------------------------------------

pro <- as.character(unique(df$protein))

for (i in 1:length(pro)){
  
  temp <- filter(df, protein == pro[i])
  
  pep <- as.character(unique(temp$peptide))
  
  for (pep_x in pep){
  
    temp_pep <- filter(temp, peptide == pep_x)
    
    cols <- c(grep(c("M_0"), colnames(temp_pep)),
              grep(c("M_3"), colnames(temp_pep)),
              grep(c("Sum"), colnames(temp_pep)))
    
    for (j in cols){
     
      g1 <- ggplot(data = temp_pep) +
              geom_point(mapping = aes(x = hrs, 
                                       y = as.numeric(temp_pep[, colnames(temp_pep)[j]]), 
                                       color = Group), 
                         size = 3, alpha = 0.8) + 
              labs(title = paste(pro[i], pep_x, sep = " -- "), 
                   subtitle = colnames(temp_pep)[j],
                   y = "Abundance", x = "Time (hrs)") +
              theme_classic()
      plot(g1)
    }
  }
}

# figure out how to fit a one phase decay to the M-0 values -- probably need 
# to do the corrected M-0 values
# maybe figure out how to put all three M_0, M_3, and Sum in the same sort of 
# space
# plotly interactive too
# make the rest of it nice

# working on aesthetics
# ggplot(data = temp_pep) +
#   geom_point(mapping = aes(x = hrs,
#                            y = as.numeric(temp_pep[, colnames(temp_pep)[7]]),
#                            fill = sample, alpha = 0.5),
#              shape = 21, size = 3, stroke = 2, alpha = 0.9) +
#   labs(title = paste(pro[5], pep_x, sep = " -- "),
#        subtitle = colnames(temp_pep)[7],
#        y = "Abundance", x = "Time (hrs)") +
#   theme_classic()

# ggsave("test.pdf", dpi = 320) could save the specific files in a .pdf using
# the paste0() as the argument of the ggsave() in the for loop.


# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/

# histone normalization of all isotopes ----------------------------------------

# average histones
df_hist <- df %>%
  group_by(protein, group, hrs) %>%
  summarize(avg_pro = mean(Sum)) %>%
  filter(protein == "Histones")

df <- cbind(df, hist_avg = df_hist$avg_pro)

# normailize all isotopes by histones 
isotopes <- c("M_0", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "Sum")

for (col in isotopes){
  df[, paste(col, "hist_norm", sep = "_")] <- df[, col] / df$hist_avg
}

# Mega Model -------------------------------------------------------------------

# initial conditions
peptide <- "SAMPLLLLLLER"

# initial total abundance at T0
t0_abun <- 1000

# fraction of labeled Leucine in the media 
per_lab <- 0.50

# rates
deg_old <- 0.0500
deg_new <- 0.0500
syn <- 50

# length of time
time <- 0:168

# model data frame -------------------------------------------------------------

mod <- mega_model(peptide, deg_old, deg_new, syn, t0_abun, per_lab, time)


# need to figure out how to work backwards to get the M_0 corrected values etc
