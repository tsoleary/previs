library(tidyverse)
setwd("C:\\Users\\PrevBeast\\Documents\\R\\Previs\\GTLEDQ")

## Import clipboard of Exact Masses from chosen scan(s) in CSV format
## (Code below omits header and sets appropriate colnames/)

# pepcsv_list <- list.files()[grep("ALQEAH", list.files())]

for (i in c(1:22)) {
maxima <- find_maxima("GTLEDQ", "LV", i)

sample <- paste("LV", i)

envelope <- find_envelopes(maxima, 1100.56, 2, rep = sample)

if (i == 1) {
  envelopes <- envelope
} else {
  envelopes <- full_join(envelope, envelopes, by = "Isotopomer")
}
}

## Find series - Before using, manually pick a m/z value from 'maxima' to be
## first ion in series. MAY NOT BE INDEX OF 1!! Depends highly on subset!
## May modify for cases in which selected starting-m/z is not M+0.
# m0 <- 936.8044
# m0ind <- which.min(abs(maxima$`m/z` - m0))
# m0 <- maxima$`m/z`[m0ind]
# envelope <- maxima[m0ind,]
# NAcount <- 0
# 
# for (i in 1:21) {
#   if (i == 1) {
#   tar <- m0 + 0.33 } else {
#   tar <- tar + 0.33
#   }
#   maxima$iondiff <- maxima$`m/z` - tar
#   mindiff <- which.min(abs(maxima$iondiff))
#   if (abs(maxima$iondiff[mindiff]) < 0.02) {
#   envelope <- full_join(envelope, maxima[mindiff, ])
#   tar <- maxima$'m/z'[mindiff] } else {
#   NAcount <- NAcount + 1
#   if (NAcount > 1) {
#     break
#   }
#   blank <- c(NA, NA, NA, NA, NA)
#   names(blank) <- names(envelope)
#   envelope <- bind_rows(envelope, blank)
#   
#   next
#   }
# }
# 
# rm(tar)



write.csv(envelopes, "GTLEDQ all LV auto.csv", row.names = FALSE)


## To do: Functionalize script. Inputs should be M0 mass, charge at the very least.
## Inputs for desired number of isotopomers and maximum deltamass might be good.
## In cases where unlabeled peaks disappear, could be useful to start from
## M+X peak.

#############################################

# find_envelope <- function(data, mass, charge, deltamax = 0.02) {
# m0ind <- which.min(abs(data$`m/z` - mass))
# m0 <- data$`m/z`[m0ind]
# envelope <- data[m0ind,]
# zreciprocal <- 1/charge
# NAcount <- 0
# 
# for (i in 1:21) {
#   if (i == 1) {
#     tar <- m0 + zreciprocal } else {
#       tar <- tar + zreciprocal
#     }
#   data$iondiff <- data$`m/z` - tar
#   mindiff <- which.min(abs(data$iondiff))
#   if (abs(data$iondiff[mindiff]) < deltamax) {
#     envelope <- full_join(envelope, data[mindiff, ])
#     tar <- data$'m/z'[mindiff] } else {
#       NAcount <- NAcount + 1
#       if (NAcount > 1) {
#         break
#       }
#       blank <- c(NA, NA, NA, NA, NA)
#       names(blank) <- names(envelope)
#       envelope <- bind_rows(envelope, blank)
#       
#       next
#     }
# }
# return(envelope)
# }
# rm(tar)

find_maxima <- function(peptide, tissue, rep) {
  
  file <- paste0(peptide, " ", rep, " ", tissue, ".csv")  
  
  mass_list <- read.csv(file) %>%
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
  mass_list$diffdiffsign <- c(0, diff(sign(diff(mass_list$Intensity))), 0)
  ## Subset
  maxima <- mass_list[mass_list$diffdiffsign == -2, 1:ncol(mass_list)] # %>%
    # .[.$`m/z` > 936.76 & .$`m/z` < (936.76 + (1/3 * 21)), ]
  ## Sets indices to '1:nrow(maxima)' to avoid base-R bug where subsets of a DF
  ## with indices higher than the number of rows introduces NAs.
  ## (Apparently could have avoided by using dplyr::filter - learning experience.)
  rownames(maxima) <- 1:nrow(maxima)
  return(maxima)
}

find_envelopes <- function(data, mass, charge, rep, deltamax = 0.02) {
  m0ind <- which.min(abs(data$`m/z` - mass))
  m0 <- data$`m/z`[m0ind]
  envelope <- data[m0ind,]
  zreciprocal <- 1/charge
  NAcount <- 0
  # Isos <- "M+0"
  
  for (i in 1:21) {
    if (i == 1) {
      tar <- m0 + zreciprocal } else {
        tar <- tar + zreciprocal
        # tempiso <- paste0("M+", i-1)
        # Isos <- c(Isos, tempiso)
      }
    data$iondiff <- data$`m/z` - tar
    mindiff <- which.min(abs(data$iondiff))
    if (abs(data$iondiff[mindiff]) < deltamax) {
      envelope <- full_join(envelope, data[mindiff, ])
      tar <- data$'m/z'[mindiff] } else {
        NAcount <- NAcount + 1
        if (NAcount > 1) {
          break
        }
        blank <- c(NA, NA, NA, NA)
        names(blank) <- names(envelope)
        envelope <- bind_rows(envelope, blank)
        
        next
      }
  }
  envelope$Isotopomer <- paste0("M+", 0:(nrow(envelope) - 1))
  envelope <- select(envelope, Isotopomer, Intensity)
  colnames(envelope) <- c("Isotopomer", rep)
  return(envelope)
}
rm(tar)


