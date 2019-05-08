# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

library(tidyverse)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")
list.files()

df <- read.csv("reprex.csv")

samples <- colnames(df)[4:length(colnames(df))]

df2 <- df %>% 
  gather(sample, abundance, samples, na.rm = TRUE) %>%
  separate("sample", c("sample", "duplicate", "time"), sep = "_")

# Correct M_0 and M_3 ----------------------------------------------------------

df_dup_avg <- df3 %>%
                group_by(protein, type, peptide, sample, time) %>%
                summarize(abundance = mean(abundance))

df5 <- filter(df_dup_avg, type == "Sum")

# Graph using ggplot2 ----------------------------------------------------------

ggplot(data = df5) +
  geom_point(mapping = aes(x = time, y = abundance, color = sample))



