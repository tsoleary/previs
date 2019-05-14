# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

#require("plotly")

library(tidyverse)
library(plotly)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")

df_messy <- read.csv("pinpoint_m_m_3_for_r.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sample", "duplicate", "hrs"), sep = "_")

# Correct M_0 and M_3 ----------------------------------------------------------

df_avg <- df %>%
  group_by(protein, type, peptide, sample, hrs) %>%
  summarize(abundance = mean(abundance))
  
df_avg$hrs <- as.numeric(df_avg$hrs)

# spreads out by "type" M_0, M_3, Sum

df_spread <- as.data.frame(spread(df_avg, "type", "abundance"))

# creates different data frames for each one.
# df_sum <- filter(df_avg, type == "Sum")
# df_m0 <- filter(df_avg, type == "M_0")
# df_m3 <- filter(df_avg, type == "M_3")

# trying to make a bunch of plots at once with a loop

pro <- as.character(unique(df_spread$protein))

for (i in 1:length(pro)){
  
  temp <- filter(df_spread, protein == pro[i])
  
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
                                       color = sample), size = 3) + 
              labs(title = paste(pro[i], pep_x, sep = " -- "), 
                   subtitle = colnames(temp_pep)[j]) +
              ylab("abundance")
      plot(g1)
    }
  }
}

# so now a few things

# 1. work on the aesthetics of the graph (rounder points?)
# 2. figure out how to fit a one phase decay to the M-0 values -- probably need 
# to do the corrected M-0 values
# 3. maybe figure out how to put all three M_0, M_3, and Sum in the same sort of 
# space
# 4. plotly interactive too
# 5. make the rest of it nice
# 

# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/

ggplotly(graph, tooltip = c("protein"))

pros <- as.character(unique(df_sum$protein))



