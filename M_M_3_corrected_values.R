# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

#require("plotly")

library(tidyverse)
library(plotly)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")
list.files()

df_messy <- read.csv("pinpoint_m_m_3_for_r.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("sample", "duplicate", "hrs"), sep = "_")

# Correct M_0 and M_3 ----------------------------------------------------------

df_avg <- df %>%
  group_by(protein, type, peptide, sample, hrs) %>%
  summarize(abundance = mean(abundance))

df_sum <- filter(df_avg, type == "Sum")
df_m0 <- filter(df_avg, type == "M_0")
df_m3 <- filter(df_avg, type == "M_3")


# Graph using ggplot2 ----------------------------------------------------------

ggplot(data = df_sum) +
  geom_point(mapping = aes(x = hrs, y = abundance, color = sample))

# Creating individual plots ---------------------------------------------------

df2 <- df %>%
  group_by(protein, peptide, sample) %>%
  ggplot + geom_point(mapping = aes(x = hrs, y = abundance, color = sample))

graph <- ggplot(data = df_sum) +
           geom_point(mapping = aes(x = hrs, y = abundance, color = sample))

ggplotly(graph, tooltip = c("protein"))




# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/

pros <- as.character(unique(df_sum$protein))



