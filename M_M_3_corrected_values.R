# Take pinpoint M and M+3 values to calculate M corrected and M+3 corrected ----

library(tidyverse)
library(plotly)

# Tidy the data ----------------------------------------------------------------

setwd("C:/Users/PrevBeast/Documents/R/Helms")

df_messy <- read.csv("All_isotopes_organized.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("group", "duplicate", "hrs"), sep = "_")

df$hrs <- as.numeric(df$hrs)

# average together the duplicate lines
df_avg <- df %>%
  group_by(protein, isotope, peptide, group, hrs) %>%
  summarize(abundance = mean(abundance))

# spreads out to columns by isotope: M_0, M_3, & Sum
df_spread <- as.data.frame(spread(df_avg, "isotope", "abundance"))

# make a reprex small data frame
df_isos <- select(df_spread, c("protein", "peptide", "group", "hrs", "M_0", 
                               "M_3", "M_6", "Sum"))
df_reprex <- df_isos %>% 
  filter(df_isos$hrs == c(0, 12, 24, 36)) 

df_reprex <- filter(df_reprex, df_reprex$protein == "MYBPC3")
df_reprex <- filter(df_reprex, df_reprex$peptide == "LNFDLIQELSHEAR")

# plot by protein, peptide, and isotope ----------------------------------------

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

# doing the M_0, M_3, Sum corrected math ---------------------------------------




