setwd("C:/Users/PrevBeast/Documents/R/Padron")

data_raw <- read.csv("DMG Unmod.csv") %>%
  full_join(., read.csv("DMG Mg.csv"), by = "File") %>%
  full_join(., read.csv("DMG Oxi.csv"), by = "File") %>%
  full_join(., read.csv("File IDs.csv"), by = "File")

write.csv(data_raw, file = "DMG with File ID.csv")
