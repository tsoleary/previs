library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/HCM/Titin Splice")

# Add linkers to Uniprot domain list
# NEEDED INPUT: Pasted domain list from Uniprot, with residue-ranges split into
# two separate columns (Start = Starting residue number, End = Ending residue
# number). Each row should be a domain.

dom_raw <- read.csv("Titin Domains Uniprot Q8WZ42.csv")

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

peps <- read.csv("HCM_duplicate_1_samples_061019 Titin 8E6 Threshold.csv")

# Try to define a useful function for per-peptide domain assignment here

peps$domain <- character(length = nrow(peps))
for (i in 1:nrow(peps)) {
  temp <- which((peps$PepEnd[i] > dom_linked$Start & peps$PepEnd[i] < dom_linked$End) |
              (peps$PepStart[i] > dom_linked$Start & peps$PepStart[i] < dom_linked$End) |
                peps$PepStart[i] < dom_linked$Start & peps$PepEnd[i] > dom_linked$End, TRUE)
  peps$domain[i] <- paste(dom_linked$termdomain[temp], collapse = "; ")
}

# Take some group averages
controls <- grep("GOL", colnames(peps))
DIS1 <- grep("DIS_1", colnames(peps))
DIS2 <- grep("DIS_2", colnames(peps))
DIS3 <- grep("DIS_3", colnames(peps))

peps$avg_control <- by_group(peps, controls, FUN = mean)
peps$avg_dis1 <- by_group(peps, DIS1, FUN = mean)
peps$avg_dis2 <- by_group(peps, DIS2, FUN = mean)
peps$avg_dis3 <- by_group(peps, DIS3, FUN = mean)

# Look at the scatter plot for a flat section..

plot(x = peps$PepStart, peps$avg_control, ylim = c(0, 2E+8))
points(peps$PepStart, peps$avg_ko)

write.csv(peps, file = "Titin rep1 annotated peptides 8E6 Threshold.csv")

# Taking in manually normalized values with summary ratio pre-calculated

normval <- read.csv("Titin rep1 annotated peptides 8E6 Threshold Norm.csv")

## Mean ratio for each domain

normval <- normval[which(!is.na(normval$DIS1.GOL) | !is.na(normval$DIS2.GOL) |
                     !is.na(normval$DIS3.GOL)), ]

for (i in 1:nrow(dom_linked)) {
  if (i == 1) {
    avrat <- data.frame(domain = character(), start = numeric(), end = numeric(),
                        overlaps = character(), AVG.DIS1 = numeric(), AVG.DIS2 = numeric(),
                        AVG.DIS3 = numeric(), AVG.Control = numeric(),
                        "DIS1.GOL" = numeric(), "DIS2.GOL" = numeric(),
                        "DIS3.GOL" = numeric(), peptides = numeric())
  }
  dommean1 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS1.GOL"]), na.rm = TRUE)
  dommean2 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS2.GOL"]), na.rm = TRUE)
  dommean3 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS3.GOL"]), na.rm = TRUE)
  ctrlmean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_control"]), na.rm = TRUE)
  dis1mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis1"]), na.rm = TRUE)
  dis2mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis2"]), na.rm = TRUE)
  dis3mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis3"]), na.rm = TRUE)
  temp <- data.frame(domain = dom_linked$Description[i], start =
                       dom_linked$Start[i], end = dom_linked$End[i],
                     overlaps = dom_linked$overlaps[i], AVG.DIS1 = dis1mean,
                     AVG.DIS2 = dis2mean, AVG.DIS3 = dis3mean,
                     AVG.Control = ctrlmean, "DIS1.GOL" = dommean1,
                     "DIS2.GOL" = dommean2, "DIS3.GOL" = dommean3,
                     peptides = length(grep(dom_linked$termdomain[i], normval$domain)))
  # if (is.nan(temp$KO.Control)) {
  #   temp$KO.Control <- NA
  #   }
  avrat <- full_join(avrat, temp)
}

write.csv(avrat, "Titin Domain DIS1 2 3 Ratio 8E6 Threshold.csv")
