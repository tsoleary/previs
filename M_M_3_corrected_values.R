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
  separate("sample", c("Group", "duplicate", "hrs"), sep = "_")

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

# so now a few things

# 2. figure out how to fit a one phase decay to the M-0 values -- probably need 
# to do the corrected M-0 values
# 3. maybe figure out how to put all three M_0, M_3, and Sum in the same sort of 
# space
# 4. plotly interactive too
# 5. make the rest of it nice
# 

# working on aesthetics
ggplot(data = temp_pep) +
  geom_point(mapping = aes(x = hrs,
                           y = as.numeric(temp_pep[, colnames(temp_pep)[7]]),
                           fill = sample, alpha = 0.5),
             shape = 21, size = 3, stroke = 2, alpha = 0.9) +
  labs(title = paste(pro[5], pep_x, sep = " -- "),
       subtitle = colnames(temp_pep)[7],
       y = "Abundance", x = "Time (hrs)") +
  theme_classic()

# ggsave("test.pdf", dpi = 320) could save the specific files in a .pdf using
# the paste0() as the argument of the ggsave() in the for loop.


# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/


# doing the M_0, M_3, Sum corrected math ---------------------------------------




