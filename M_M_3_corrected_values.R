# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

#require("plotly")

library(tidyverse)
library(plotly)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")

df_messy <- read.csv("reprex.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sample", "duplicate", "hrs"), sep = "_")

# Correct M_0 and M_3 ----------------------------------------------------------

df_avg <- df %>%
  group_by(protein, type, peptide, sample, hrs) %>%
  summarize(abundance = mean(abundance))

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
  
  loop_vector <- c(grep(c("M_0"), colnames(temp)),
                   grep(c("M_3"), colnames(temp)),
                   grep(c("Sum"), colnames(temp)))
  
  for (j in loop_vector){
   
     ggplot(data = temp) +
       geom_point(mapping = aes(x = hrs, 
                                y = as.numeric(temp[, colnames(temp)[j]]), 
                                color = sample)) + 
       labs(title = pro[i]) +
       ylab(abundance)
  
  }
  
}


# Graph using ggplot2 ----------------------------------------------------------

ggplot(data = df_sum) +
  geom_point(mapping = aes(x = hrs, y = abundance, color = sample)) +
  labs(title = "data")


# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/

ggplotly(graph, tooltip = c("protein"))

pros <- as.character(unique(df_sum$protein))



