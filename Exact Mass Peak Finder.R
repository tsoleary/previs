library(tidyverse)
setwd("C:\\Users\\PrevBeast\\Documents\\R\\Previs\\D3L mouse study\\Shared Myh Peptides\\VQLL\\Extraction")

## In case you're on your second or third peptide..
rm(envelope, envelopes, maxima, maxtemp, compiled_maxima)

## Pep: Peptide identifier as used in source filenames
pep <- "VQLL"
## mz: Peptide m/z value
mz <- 861.98
## z: Peptide charge-state
z <- 2
## Tis: Source-tissue as used in source filenames
tis <- "LV"
## Samps: If not all sample indices are used (as in extraction test), this should
## be a vector containing all sample indices of interest
Samps <- c(7, 8, 19, 20)
## Fracs: Vector of fraction-identifiers in extraction expt.
Fracs <- c("P", "A", "B")


## For-loop runs find_maxima to find all local maxima in each clipboard file.
## After each file's list of maxima is found, find_envelope attempts to find
## the isotopic envelope of interest given the maxima df, m/z value, and charge
## of the peptide. 
## Given the vastly different x-axis values between sample clipboards, I append
## a bunch of NA values to the compiled_maxima df so all columns have the same
## length. This and 'envelopes' are the exports of note. The former mainly
## being useful for manual validation.


## For whole-tissue D3 samples:
# for (i in Samps) {
# maxima <- find_maxima(pep, tis, i)
# sample <- paste(tis, i)
# maxtemp <- maxima[ , 1:2]
# colnames(maxtemp) <- paste(colnames(maxtemp), sample)
# 
# envelope <- find_envelope(maxima, mz, z, rep = sample)
# 
# if (i == Samps[1]) {
#   envelopes <- envelope
#   compiled_maxima <- maxtemp
# } else {
#   envelopes <- full_join(envelope, envelopes, by = "Isotopomer")
#     if (nrow(maxtemp) > nrow(compiled_maxima)) {
#       compiled_maxima[(nrow(compiled_maxima)+1):nrow(maxtemp), ] <- rep(NA,
#       times = (nrow(maxtemp) - nrow(compiled_maxima)))
#     }
#     if (nrow(compiled_maxima) > nrow(maxtemp)) {
#       maxtemp[(nrow(maxtemp)+1):nrow(compiled_maxima), ] <- rep(NA,
#       times = (nrow(compiled_maxima) - nrow(maxtemp)))
#   }
#   compiled_maxima <- bind_cols(maxtemp, compiled_maxima)
# }
# }

## For Extractions: 
for (i in Samps) {
  for (j in Fracs) {
    maxima <- find_maxima(pep, tis, i, extraction = TRUE, fraction = j)
    sample <- paste(tis, i, j)
    maxtemp <- maxima[ , 1:2]
    colnames(maxtemp) <- paste(colnames(maxtemp), sample)
    
    envelope <- find_envelope(maxima, mz, z, rep = sample)
    
    if (i == Samps[1] & j == Fracs[1]) {
      envelopes <- envelope
      compiled_maxima <- maxtemp
    } else {
      envelopes <- full_join(envelope, envelopes, by = "Isotopomer")
      if (nrow(maxtemp) > nrow(compiled_maxima)) {
        compiled_maxima[(nrow(compiled_maxima)+1):nrow(maxtemp), ] <- rep(NA,
                                                                          times = (nrow(maxtemp) - nrow(compiled_maxima)))
      }
      if (nrow(compiled_maxima) > nrow(maxtemp)) {
        maxtemp[(nrow(maxtemp)+1):nrow(compiled_maxima), ] <- rep(NA,
                                                                  times = (nrow(compiled_maxima) - nrow(maxtemp)))
      }
      compiled_maxima <- bind_cols(maxtemp, compiled_maxima)
    }
  }
}


## Exports two separate CSVs of 'envelopes' and 'compiled_maxima' respectively

write.csv(envelopes, "VQLL LV Extraction test auto.csv", row.names = FALSE)
write.csv(compiled_maxima, "VQLL LV Extraction test maxima.csv", row.names = FALSE)


## Run these first! This is the quick-and-dirty heights script (as opposed to
## using actual peak-areas, which would be more rigorous), so I'm hestitant
## to put these in the actual proteomics_functions.R script.

find_maxima <- function(peptide, tissue, rep, extraction = FALSE, fraction = NULL) {
  if (extraction == TRUE){
    file <- paste0(peptide, " ", rep, " ", tissue, " ", fraction, ".csv")
  } else {
    file <- paste0(peptide, " ", rep, " ", tissue, ".csv")  
  }
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

find_envelope <- function(data, mass, charge, rep, deltamax = 0.02, NAcounter = TRUE) {
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
          if (NAcounter == TRUE) {
          NAcount <- NAcount + 1
          if (NAcount > 1) {
            if ((i < 4) | (i < 7 & any(!is.na(envelope[, 2])))) {
            NAcount <- 0
            blank <- rep(NA, times = ncol(envelope))
            names(blank) <- names(envelope)
            envelope <- bind_rows(envelope, blank)
            next
            }
            break
          }
          blank <- rep(NA, times = ncol(envelope))
          names(blank) <- names(envelope)
          envelope <- bind_rows(envelope, blank)
          
          next
        }
      }
  }
  envelope$Isotopomer <- paste0("M+", 0:(nrow(envelope) - 1))
  envelope <- select(envelope, Isotopomer, Intensity)
  colnames(envelope) <- c("Isotopomer", rep)
  return(envelope)
}
rm(tar)

