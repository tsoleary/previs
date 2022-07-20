library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/Esser/Titin 1.5E7 Threshold")

# Add linkers to Uniprot domain list
# NEEDED INPUT: Pasted domain list from Uniprot, with residue-ranges split into
# two separate columns (Start = Starting residue number, End = Ending residue
# number). Each row should be a domain.

dom_raw <- read.csv("Exons NCBI only.csv")

# Start and End columns must be numeric for the operations to follow
dom_raw$Start <- as.numeric(dom_raw$Start)
dom_raw$End <- as.numeric(dom_raw$End)

# Now we add some linkers
for (i in rownames(dom_raw)) {
  
  thisrow <- as.numeric(i)
  nextrow <- (thisrow + 1)
  nextrow2 <- (thisrow + 2)
  lastrow <- (thisrow - 1)
  
  nextstart <- dom_raw[nextrow, 1]
  nextstart2 <- dom_raw[nextrow2, 1]
  
  thisend <- dom_raw[i, 2]
  nextend <- dom_raw[nextrow, 2]
  nextend2 <- dom_raw[nextrow2, 2]
  lastend <- dom_raw[lastrow, 2]
  
  if (i == 1) {
    count <- 0
    linkerlist <- data.frame(Start = numeric(), End = numeric(), 
                             Description = character())
  }
  
  if (i == length(rownames(dom_raw))){
    print(paste("here's the end, it's", thisend))
    break
  }
  
  if ((nextstart - thisend) == 1) {
  print(paste("no linker after", dom_raw[i, 3]))
  } 
  if (i != 1) {
    if ((thisend - lastend) < 0) {
    print(paste("Nested domain at", thisrow))
    next
    }
    
  }
  if (nextstart - thisend > 1) {
    print(paste("Linker", (count + 1)))
    count <- (count + 1)
    linkname <- paste("Linker", count)
    linktemp <- data.frame(Start = (thisend + 1), End = (nextstart - 1),
                           Description = linkname)
    linkerlist <- full_join(linkerlist, linktemp)
  }
  if (thisrow != nrow(dom_raw) & nextstart - thisend < 0 & nextstart2 - thisend > 1
      & nextend - thisend < 0) {
    print(paste("Linker", (count + 1)))
    count <- (count + 1)
    linkname <- paste("Linker", count)
    linktemp <- data.frame(Start = (thisend + 1), End = (nextstart2 - 1),
                           Description = linkname)
    linkerlist <- full_join(linkerlist, linktemp)
  }
}

dom_linked <- full_join(dom_raw, linkerlist) %>%
arrange(., Start)


# For averaging later, we need to add a way of distinguishing the end of a domain
# name (i.e., when pulling up all rows where "domain" contains "Ig-like 1",
# we don't want Ig-like 10, Ig-like 11, Ig-like 100, etc. showing up.

dom_linked$termdomain <- paste0(dom_linked$Description, "xx")

# List overlapping domains while we're at it

dom_linked$overlaps <- character(length = nrow(dom_linked))
for (i in 1:nrow(dom_linked)) {
  temp <- which((dom_linked$End[i] > dom_linked$Start & dom_linked$End[i] < dom_linked$End) |
                  (dom_linked$Start[i] > dom_linked$Start & dom_linked$Start[i] < dom_linked$End) |
                  dom_linked$Start[i] < dom_linked$Start & dom_linked$End[i] > dom_linked$End, TRUE)
  dom_linked$overlaps[i] <- paste(dom_linked$Description[temp], collapse = "; ")
}

# Now we have a data.frame with all Uniprot domains, and linkers interspersed
# appropriately. Next, import a .csv file with Start, End, and abundances for
# peptides detected via PD

peps <- read.csv("Mus Ttn 203 TA All Peptides threshold.csv")

# Try to define a useful function for per-peptide domain assignment here

peps$domain <- character(length = nrow(peps))
for (i in 1:nrow(peps)) {
  temp <- which((peps$PepEnd[i] > dom_linked$Start & peps$PepEnd[i] < dom_linked$End) |
              (peps$PepStart[i] > dom_linked$Start & peps$PepStart[i] < dom_linked$End) |
                peps$PepStart[i] < dom_linked$Start & peps$PepEnd[i] > dom_linked$End, TRUE)
  peps$domain[i] <- paste(dom_linked$termdomain[temp], collapse = "; ")
}

# Take some group averages
controls <- grep("Control", colnames(peps))
KOs <- grep("KO", colnames(peps))

peps$avg_control <- by_group(peps, controls, FUN = median)
peps$avg_ko <- by_group(peps, KOs, FUN = median)

# Look at the scatter plot for a flat section..

plot(x = peps$PepStart, peps$avg_control, ylim = c(0, 2.5E+8))
points(peps$PepStart, peps$avg_ko)

write.csv(peps, file = "TA Titin 203 annotated peptides.csv")

# Taking in manually normalized values with summary ratio pre-calculated

normval <- read.csv("Renorm Aband Top3 WRT Cronos.csv")

## Mean ratio for each domain

normval <- normval[which(!is.na(normval$KO.Control)), ]

# for (i in 1:nrow(dom_linked)) {
#   if (i == 1) {
#     avrat <- data.frame(domain = character(), start = numeric(), end = numeric(),
#                         overlaps = character(), AVG.KO = numeric(), AVG.Control = numeric(),
#                         "KO.Control" = numeric(), peptides = numeric())
#   }
#   dommean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "KO.Control"]))
#   ctrlmean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_control"]))
#   komean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_ko"]))
#   temp <- data.frame(domain = dom_linked$Description[i], start =
#                        dom_linked$Start[i], end = dom_linked$End[i],
#                      overlaps = dom_linked$overlaps[i], AVG.KO = komean, 
#                      AVG.Control = ctrlmean, "KO.Control" = dommean,
#                      peptides = length(grep(dom_linked$termdomain[i], normval$domain)))
#   if (is.nan(temp$KO.Control)) {
#     temp$KO.Control <- NA
#     }
#   avrat <- full_join(avrat, temp)
# }

for (i in 1:length(unique(normval$domain))) {
  if (i == 1) {
    avrat <- data.frame(domain = character(), start = numeric(), end = numeric(),
                        overlaps = character(), AVG.KO = numeric(), AVG.Control = numeric(),
                        AVG.KO.Norm = numeric(), AVG.Control.Norm = numeric(),
                        "KO.Control" = numeric(), peptides = numeric())
  }
  dommean <- mean(as.numeric(normval[which(normval$domain == unique(normval$domain)[i]), "KO.Control"]))
  ctrlmean <- mean(as.numeric(normval[which(normval$domain == unique(normval$domain)[i]), "avg_control"]))
  komean <- mean(as.numeric(normval[which(normval$domain == unique(normval$domain)[i]), "avg_ko"]))
  ctrlnormmean <- mean(as.numeric(normval[which(normval$domain == unique(normval$domain)[i]), "avg_control_norm"]))
  konormmean <- mean(as.numeric(normval[which(normval$domain == unique(normval$domain)[i]), "avg_ko_norm"]))
  if(str_detect(unique(normval$domain)[i], ";") == FALSE) {
    domind <- which(dom_linked$termdomain == unique(normval$domain)[i])
    temp <- data.frame(domain = dom_linked$termdomain[domind], start =
                         dom_linked$Start[domind], end = dom_linked$End[domind],
                       overlaps = dom_linked$overlaps[domind], AVG.KO = komean, 
                       AVG.Control = ctrlmean,
                       AVG.KO.Norm = konormmean, AVG.Control.Norm = ctrlnormmean, 
                       "KO.Control" = dommean,
                       peptides = length(which(normval$domain == unique(normval$domain)[i])))
  }
  if(str_detect(unique(normval$domain)[i], ";") == TRUE) {
    pepind <- which(normval$domain == unique(normval$domain)[i])
    temp <- data.frame(domain = unique(normval$domain)[i], start =
                         min(normval$PepStart[pepind]), end = max(normval$PepEnd[pepind]),
                       overlaps = unique(normval$domain)[i], AVG.KO = komean, 
                       AVG.Control = ctrlmean, "KO.Control" = dommean,
                       peptides = length(pepind))
  }
  if (is.nan(temp$KO.Control)) {
    temp$KO.Control <- NA
  }
  avrat <- full_join(avrat, temp)
}

# Now, stitch in domains not detected for complete list (NA values for Ratio)

for (i in 1:nrow(dom_linked)) {
  if(is.element(dom_linked$termdomain[i], avrat$domain)) {
    next
  }
  if(!is.element(dom_linked$termdomain[i], avrat$domain)) {
    temp <- data.frame(domain = dom_linked$termdomain[i], start =
                         dom_linked$Start[i], end = dom_linked$End[i],
                       overlaps = dom_linked$overlaps[i], AVG.KO = NA, 
                       AVG.Control = NA, "KO.Control" = NA,
                       peptides = 0)
    if(i == 1) {
      avratfull <- full_join(avrat, temp)
    } else {
      avratfull <- full_join(avratfull, temp)
    }
    
  }
}

write.csv(avratfull, "Titin Aband Top3 norm Exon view norms.csv")
