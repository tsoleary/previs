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

# # Take some group averages
controls <- grep("Control", colnames(peps))
g1 <- grep("KO", colnames(peps))

# RELIC: Might not need again
peps$avg_control <- by_group(peps, controls, FUN = mean)
peps$avg_g1 <- by_group(peps, g1, FUN = mean)

# Look at the scatter plot for a flat section..

plot(x = peps$PepStart, peps$avg_control, ylim = c(0, 2E+8))
points(peps$PepStart, peps$avg_ko)

write.csv(peps, file = "Titin rep1 annotated peptides 8E6 Threshold rev.csv")

# Taking in manually normalized values with summary ratio pre-calculated

normval <- read.csv("Renorm Aband WRT Cronos.csv")

# Median, sd, & ratio of peptides ----------------------------------------------

normcols <- which(str_detect(colnames(normval), "_norm"))


group1 <- grep("KO", colnames(normval)) %>%
  intersect(., normcols)
ctrl <- grep("Control", colnames(normval)) %>%
  intersect(., normcols)

normval$group1_med <- by_group(normval, group1)
normval$ctrl_med <- by_group(normval, ctrl)

# Standard deviation between samples for each peptide
normval$group1_sd <- by_group(normval, group1, FUN = sd)
normval$ctrl_sd <- by_group(normval, ctrl, FUN = sd)

normval$KO.Control <- abun_ratio(normval, "group1_med")

normval$group1_df <- count_df(normval, group1)
normval$ctrl_df <- count_df(normval, ctrl)

# Removing rows with NA values for the all ratios
normval <- normval[which(!is.na(normval$KO.Control)),]

# Remove outliers --------------------------------------------------------------

pro_out1 <- by_variable(normval, "KO.Control", VAR = "domain") %>% as.data.frame %>%
  rownames_to_column("domain")

normval <- rm_outliers_dom(normval, pro_out1, "KO.Control")


# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(normval) %>%
  group_by(domain) %>%
  # top_n(n = pep_top, wt = ctrl_med) %>%
  as.data.frame

# Regional Averages -------------------------------------------------------------

group_names <- c("group1_med", "ctrl_med")

regional <- by_variable(data_top, group_names, VAR = "domain") %>%
  as.data.frame %>%
  rownames_to_column("domain")

regional$domain <-
  regional$domain %>%
  as.character

data_top$group1_sd_df <- square_x_df(data_top, "group1_sd", "group1_df")
data_top$ctrl_sd_df <- square_x_df(data_top, "ctrl_sd", "ctrl_df")

# Creating a temp data frame to calculate the ratio1 sd
temp_df <- NULL 
sd_sum <- c("group1_sd_df", "ctrl_sd_df", "group1_df", "ctrl_df")
temp_df <- by_variable(data_top, sd_sum, FUN = sum, VAR = "domain") %>% as.data.frame

group1_pooled_sd <- temp_df$group1_sd_df / temp_df$group1_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
regional <- cbind(regional, group1_pooled_sd)
regional <- cbind(regional, ctrl_pooled_sd)

# Grouped relative protein abundance ratio
regional$KO.Control <- as.numeric(by_variable(normval, "KO.Control", VAR = "domain"))

# # Adding grouped ratio if regions with NA ratio values are removed
# withnan <- as.numeric(by_variable(normval, "KO.Control", VAR = "domain"))
# 
# regional$KO.Control <- withnan[!is.nan(withnan)]

# Standard deviation relative protein abundance ratio
normval$KO.Control_sd <- ratio_sd(normval, "KO.Control", "group1_sd", "group1_med", 
                           "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio1
normval$KO.Control_df <- normval$group1_df + normval$ctrl_df
normval$KO.Control_sd_df <- square_x_df(normval, "KO.Control_sd", "KO.Control_df")

# ratio1 sd_calc
sd_calc <- c("KO.Control_df", "KO.Control_sd_df")
temp_df <- by_variable(normval, sd_calc, FUN = sum, VAR = "domain") %>% as.data.frame

KO.Control_sd <- temp_df$KO.Control_sd_df / temp_df$KO.Control_df
regional <- cbind(regional,  KO.Control_sd)
# NOTE !!!! sd > 0 condition is only to make this work with NA exons (I.E., ratio undefined)
# regional <- cbind(regional,  KO.Control_sd[KO.Control_sd > 0])


# Statistics -------------------------------------------------------------------

log_norm <- log(normval[, c(group1, ctrl)])
colnames(log_norm) <- paste(colnames(log_norm), sep = "_", "log")
normval <- cbind(normval, log_norm)

log_cols <- grep("log", colnames(normval))
group1_log_cols <- intersect(log_cols, grep("KO", colnames(normval)))
ctrl_log_cols <- intersect(log_cols, grep("Control", colnames(normval)))


# Problem if different number of samples be in each group
length(group1_log_cols)
length(ctrl_log_cols)

# add 1 dummy columns to Control
normval$dummy1_Control_norm_log <- NA
# normval$dummy2_Control_norm_log <- NA
# normval$dummy3_Control_norm_log <- NA
# normval$dummy4_Control_norm_log <- NA
# normval$dummy5_Control_norm_log <- NA
# normval$dummy6_Control_norm_log <- NA
# normval$dummy7_Control_norm_log <- NA
# normval$dummy8_Control_norm_log <- NA
# normval$dummy9_Control_norm_log <- NA
# normval$dummy10_Control_norm_log <- NA
# normval$dummy11_Control_norm_log <- NA
# normval$dummy12_Control_norm_log <- NA
# normval$dummy13_Control_norm_log <- NA
# normval$dummy14_Control_norm_log <- NA
# normval$dummy15_Control_norm_log <- NA
# normval$dummy16_Control_norm_log <- NA
# normval$dummy17_Control_norm_log <- NA
# normval$dummy18_Control_norm_log <- NA
# normval$dummy19_Control_norm_log <- NA
# normval$dummy20_Control_norm_log <- NA

log_cols <- grep("log", colnames(normval))
group1_log_cols <- intersect(log_cols, grep("KO", colnames(normval)))
ctrl_log_cols <- intersect(log_cols, grep("Control", colnames(normval)))

stacked <- data.frame(normval[, "domain"], 
                      stack(normval[, group1_log_cols]),
                      stack(normval[, ctrl_log_cols]))


colnames(stacked) <- c("domain", "group1_log", "samp1",
                       "ctrl_log", "ctrl")

# t-test for group1
p_vals1 <- pval_ttest(stacked, "group1_log", "ctrl_log", col = "domain")
colnames(p_vals1) <- c("domain", "pvalue1")

p_vals1$domain <-
  p_vals1$domain %>% 
  as.character

regional <- full_join(regional, p_vals1, by = "domain")


# Minimum number of peptides for each protein group ----------------------------
# min_pep <- 5 
regional$peptides <- table(normval$domain)
# regional <- filter(regional, regional$peptides >= min_pep)

write.csv(regional, "Titin Exon abundance ratios3.csv")
write.csv(normval, "Titin Exon peptide abundance ratios3.csv")

# ## Mean ratio for each domain
# 
# # normval <- normval[which(!is.na(normval$KO.Control) | !is.na(normval$DIS2.GOL) |
# #                      !is.na(normval$DIS3.GOL)), ]
# 
# for (i in 1:nrow(dom_linked)) {
#   if (i == 1) {
#     avrat <- data.frame(domain = character(), start = numeric(), end = numeric(),
#                         overlaps = character(), AVG.DIS1 = numeric(), AVG.DIS2 = numeric(),
#                         AVG.DIS3 = numeric(), AVG.Control = numeric(),
#                         "KO.Control" = numeric(), "DIS2.GOL" = numeric(),
#                         "DIS3.GOL" = numeric(), peptides = numeric())
#   }
#   dommean1 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "KO.Control"]), na.rm = TRUE)
#   dommean2 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS2.GOL"]), na.rm = TRUE)
#   dommean3 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS3.GOL"]), na.rm = TRUE)
#   ctrlmean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_control"]), na.rm = TRUE)
#   dis1mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis1"]), na.rm = TRUE)
#   dis2mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis2"]), na.rm = TRUE)
#   dis3mean <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "avg_dis3"]), na.rm = TRUE)
#   temp <- data.frame(domain = dom_linked$Description[i], start =
#                        dom_linked$Start[i], end = dom_linked$End[i],
#                      overlaps = dom_linked$overlaps[i], AVG.DIS1 = dis1mean,
#                      AVG.DIS2 = dis2mean, AVG.DIS3 = dis3mean,
#                      AVG.Control = ctrlmean, "KO.Control" = dommean1,
#                      "DIS2.GOL" = dommean2, "DIS3.GOL" = dommean3,
#                      peptides = length(grep(dom_linked$termdomain[i], normval$domain)))
#   # if (is.nan(temp$KO.Control)) {
#   #   temp$KO.Control <- NA
#   #   }
#   avrat <- full_join(avrat, temp)
# }
# 
# write.csv(avrat, "Titin Domain DIS1 2 3 Ratio 8E6 Threshold.csv")
