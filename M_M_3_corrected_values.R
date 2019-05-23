# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(plotly)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")

df_messy <- read.csv("All_isotopes_organized_select_reprex.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df_tidy <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("group", "duplicate", "hrs"), sep = "_")

df_tidy$hrs <- as.numeric(df_tidy$hrs)

# average together the duplicate lines
df_avg <- df_tidy %>%
  group_by(protein, isotope, peptide, group, hrs) %>%
  summarize(abundance = mean(abundance))

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

# percent distribution of each isotope -----------------------------------------

# count number of leucines in a peptide
df$num_L <- str_count(df$peptide, pattern = "L")

# function to determine the distribution 
# p <- percent labeled leucine in media
# L <- number of leucines in peptide
# mass isotopomer (M_0 = 0, M_3 = 1, M_6 = 2, ...)
# (factorial(L)/(factorial(L-i)*factorial(i)))*(p^i)*(1-p)^(L-i)


# calculates the isotopic distribution of 
iso_dist <- function (pep, p = 0.5){
  
  L <- str_count(pep, pattern = "L")
  isos <- 0:L
  dis <- NULL

  for (i in isos){
    temp <- (factorial(L)/(factorial(L-i)*factorial(i)))*(p^i)*(1-p)^(L-i)
    dis <- c(dis, temp)
  }
  return(dis)
}

iso_dist("SAMPLLLLER")


# M_0, M_3, M_6 correction -----------------------------------------------------

# M_0 correction 

for (i in length(unique(hrs))){
  
}


# M_3 and M_6 correction


# Mega Model -------------------------------------------------------------------

# initial conditions
peptide <- "SAMPLLER"
num_L <- str_count(peptide, pattern = "L")
t0_abun <- 1000
per_lab <- 0.50

# rates
deg_old <- 0.0500
deg_new <- 0.0500


# data frame 
time <- 0:168

# nat_iso distribution
iso <- 0:8
m_z <- c(916.49, 917.49, 918.50, 919.50, 920.50, 921.50, 922.50, 923.50, 924.50)
per_total <- c(56.83, 28.41, 10.87, 3.05, 0.69, 0.13, 0.02, 0.00, 0.00)
per_max <- c(100.00, 49.99, 19.12, 5.36, 1.21, 0.23, 0.04, 0.01, 0.00)

nat_iso <- data.frame("isotope" = iso, "m/z" = m_z, "percent_total" = per_total,
                      "percent_max" = per_max)

#Rdisop install (takes a long time)
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rdisop")

library("Rdisop")

aa_form <- read.csv("aa_molecular_formula.csv")

L_row <- grep("L", aa_form$letter)

Leu <- as.character(aa_form[L_row, "molecular_formula"])

Leu2 <- getMolecule(Leu, z=1)

getIsotope(Leu, 1)
getIsotope(Leu2, seq(1,5))
















# Model data frame -------------------------------------------------------------
mod <- data.frame("time" = time)

