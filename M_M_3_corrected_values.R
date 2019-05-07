# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

library(tidyverse)

# Tidy the data ----------------------------------------------------------------

getwd()
setwd("C:/Users/PrevBeast/Documents/R/Helms")
list.files()

df <- read.csv("reprex.csv")

samples <- c("DF_A_T0", "DF_A_T2", "DF_B_T0", "DF_B_T2", 
             "L2_A_T0", "L2_A_T2", "L2_B_T0", "L2_B_T2")

df2 <- df %>% 
  gather(sample, abundance, samples, na.rm = TRUE)

df3 <- separate(df2, "sample", c("sample", "duplicate", "time"), sep = "_")

# Correct M_0 and M_3 ----------------------------------------------------------

df4 <- df3 %>%
         group_by(protein, type, peptide, sample, time) %>%
         summarize(abundance = mean(abundance))

df5 <- filter(df4, type == "Sum")


#okay so the next things that I want to do are graph the sums over time

# so I need to convert the times from T0 T2 to 0 hours, 12 hours etc....


