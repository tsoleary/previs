# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(plotly)
library(Rdisop)

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

aa_for <- read.csv("aa_molecular_formula.csv")

# initial conditions -----------------------------------------------------------

# peptide specific
peptide <- "SAMPLLER"
D3_pep_dist <- iso_dist(peptide, p = 0.5)

# natural isotopic distribution of the peptide
nat_iso <- pep_iso(peptide, charge = 1)

# initial total abundance at T0
t0_abun <- 1000

# rates
deg_old <- 0.0500
deg_new <- 0.0500
syn <- 50

# length of time
time <- 0:168

# model data frame -------------------------------------------------------------

mod <- mega_model("SAMPLLER", deg_old = 0.05, deg_new = 0.05, syn = 50, 
                  t0_abun = 1000, per_lab = 0.5, time = c(0:168))


# need to create a function that combines all individual isotopes into one total

mutate(mod, sum = colnames(mod)[grep("M_0", colnames(mod))], sum)



df_test <- mod
colnames(df_test) <- gsub("^.*(M_[0-9]+)", "\\1", colnames(df_test))

mod_comb <- sapply(split.default(df_test, colnames(df_test)), 
                   rowSums, na.rm = TRUE)


y <- colnames(mod_comb)


# combine the individual isotopomers from different pools back to a total 
combine_iso <- function (dat){
  
  # removes the M_ from each column name?
  colnames(dat) <- gsub("^.*(M_[0-9]+)", "\\1", colnames(dat))
  
  # combine all columns that have the same name by their 
  df <- sapply(split.default(dat, colnames(dat)), 
                     rowSums, na.rm = TRUE)
  
  # add the word total to every colname with a M_0 pattern
  colnames(df) <- gsub("^.*(M_[0-9]+)", 
                             paste0("\\1", "_total"), colnames(df))  
  
  # remove pook and time from the df
  df <- df[, grep(c("M_"), colnames(df))]
  
  # pading with a 0 to allow sorting of the columns 
  colnames(df) <- str_pad(gsub("M_", "", colnames(df)), 8, pad = "0")
  
  df <- df[, sort(colnames(df))]
  
  colnames(df) <- c(paste0("M_", colnames(df)))
                    
  return(df)
  
}





# the columns come out of order because it is lexographic indexing not numerical
# need to pad the numbers with a 0

xy <- str_pad(gsub("M_", "", colnames(x)), 8, pad = "0")

sort(x)

# need to figure out how to work backwards to get the M_0 corrected values etc






