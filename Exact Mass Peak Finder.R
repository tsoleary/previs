library(tidyverse)
setwd("C:\\Users\\PrevBeast\\Documents\\R\\Previs")

## Import clipboard of Exact Masses from chosen scan(s) in CSV format
## (Code below omits header and sets appropriate colnames/)

mass_list <- read.csv("ALQEAH 5 LV.csv") %>%
    .[8:nrow(.), ]
colnames(mass_list) <- c("m/z", "Intensity")

## Decoy list to verify operations working as intended
# mass_list2 <- mass_list

## Coerces factor-class columns to numeric (character first because of factor-
## to-numeric weirdness.)
mass_list$`m/z` <- as.numeric(as.character(mass_list$`m/z`))
mass_list$Intensity <- as.numeric(as.character(mass_list$Intensity))
##  Diffdiffsign operation is key for finding maxima here. Theoretically does not
##  need to be a column. Might change later.
mass_list$diffsign <- c(sign(diff(mass_list$Intensity)), 0)
mass_list$diffdiffsign <- c(0, diff(sign(diff(mass_list$Intensity))), 0)
## Subset
maxima <- mass_list[mass_list$diffdiffsign == -2, 1:ncol(mass_list)] %>%
  .[.$`m/z` > 936.76 & .$`m/z` < (936.76 + (1/3 * 21)), ]
## Sets indices to '1:nrow(maxima)' to avoid base-R bug where subsets of a DF
## with indices higher than the number of rows introduces NAs.
## (Apparently could have avoided by using dplyr::filter - learning experience.)
rownames(maxima) <- 1:nrow(maxima)
  

## Find series - Before using, manually pick a m/z value from 'maxima' to be
## first ion in series. MAY NOT BE INDEX OF 1!! Depends highly on subset!
## May modify for cases in which selected starting-m/z is not M+0.
m0 <- 936.8176
envelope <- maxima[1,]

for (i in 1:21) {
  if (i == 1) {
  tar <- m0 + 0.33 } else {
  tar <- tar + 0.33
  }
  maxima$iondiff <- maxima$`m/z` - tar
  mindiff <- which.min(abs(maxima$iondiff))
  if (abs(maxima$iondiff[mindiff]) < 0.025) {
  envelope <- full_join(envelope, maxima[mindiff, ])
  tar <- maxima$'m/z'[mindiff] } else {
  next
  }
}

rm(tar)

write.csv(envelope, "test series ALQEAH LV 5.csv", row.names = FALSE)


## To do: Functionalize script. Inputs should be M0 mass, charge at the very least.
## Inputs for desired number of isotopomers and maximum deltamass might be good.
## In cases where unlabeled peaks disappear, could be useful to start from
## M+X peak.

#############################################

