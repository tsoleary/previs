library(tidyverse)
source("C:/Users/PrevBeast/Documents/GitHub/Previs/proteomics_functions.R")

setwd("C:/Users/PrevBeast/Documents/R/HCM/Titin Splice/Revised ratio and stats")

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

# # Take some group averages
# controls <- grep("GOL", colnames(peps))
# DIS1 <- grep("DIS_1", colnames(peps))
# DIS2 <- grep("DIS_2", colnames(peps))
# DIS3 <- grep("DIS_3", colnames(peps))

## RELIC: Might not need again
# peps$avg_control <- by_group(peps, controls, FUN = mean)
# peps$avg_dis1 <- by_group(peps, DIS1, FUN = mean)
# peps$avg_dis2 <- by_group(peps, DIS2, FUN = mean)
# peps$avg_dis3 <- by_group(peps, DIS3, FUN = mean)

# # Look at the scatter plot for a flat section..
# 
# plot(x = peps$PepStart, peps$avg_control, ylim = c(0, 2E+8))
# points(peps$PepStart, peps$avg_ko)

write.csv(peps, file = "Titin rep1 annotated peptides 8E6 Threshold rev.csv")

# Taking in manually normalized values with summary ratio pre-calculated

normval <- read.csv("Titin rep1 annotated peptides 8E6 Threshold rev Norm PEVK subdivision.csv")

# Median, sd, & ratio of peptides ----------------------------------------------

normcols <- which(str_detect(colnames(normval), "_norm"))


group1 <- grep("DIS_1", colnames(normval)) %>%
  intersect(., normcols)
group2 <- grep("DIS_2", colnames(normval)) %>%
  intersect(., normcols)
group3 <- grep("DIS_3", colnames(normval)) %>%
  intersect(., normcols)
ctrl <- grep("GOL", colnames(normval)) %>%
  intersect(., normcols)

normval$group1_med <- by_group(normval, group1)
normval$group2_med <- by_group(normval, group2)
normval$group3_med <- by_group(normval, group3)
normval$ctrl_med <- by_group(normval, ctrl)

# Standard deviation between samples for each peptide
normval$group1_sd <- by_group(normval, group1, FUN = sd)
normval$group2_sd <- by_group(normval, group2, FUN = sd)
normval$group3_sd <- by_group(normval, group3, FUN = sd)
normval$ctrl_sd <- by_group(normval, ctrl, FUN = sd)

normval$DIS1.GOL <- abun_ratio(normval, "group1_med")
normval$DIS2.GOL <- abun_ratio(normval, "group2_med")
normval$DIS3.GOL <- abun_ratio(normval, "group3_med")

normval$group1_df <- count_df(normval, group1)
normval$group2_df <- count_df(normval, group2)
normval$group3_df <- count_df(normval, group3)
normval$ctrl_df <- count_df(normval, ctrl)

# Removing rows with NA values for the all ratios
normval <- normval[which(!is.na(normval$DIS1.GOL) | !is.na(normval$DIS2.GOL) |
                     !is.na(normval$DIS3.GOL)) ,]

# Remove outliers --------------------------------------------------------------

pro_out1 <- by_variable(normval, "DIS1.GOL", VAR = "region") %>% as.data.frame %>%
  rownames_to_column("region")

normval <- rm_outliers_reg(normval, pro_out1, "DIS1.GOL")

pro_out2 <- by_variable(normval, "DIS2.GOL", VAR = "region") %>% as.data.frame %>%
  rownames_to_column("region")

normval <- rm_outliers_reg(normval, pro_out2, "DIS2.GOL")

pro_out3 <- by_variable(normval, "DIS3.GOL", VAR = "region") %>% as.data.frame %>%
  rownames_to_column("region")

normval <- rm_outliers_reg(normval, pro_out3, "DIS3.GOL")

# Data frame with only top few ionizing peptides -------------------------------

pep_top <- 3
data_top <-
  tbl_df(normval) %>%
  group_by(region) %>%
  top_n(n = pep_top, wt = ctrl_med) %>%
  as.data.frame

# Regional Averages -------------------------------------------------------------

group_names <- c("group1_med", "group2_med", "group3_med", "ctrl_med")

regional <- by_variable(data_top, group_names, VAR = "region") %>%
  as.data.frame %>%
  rownames_to_column("region")

regional$region <-
  regional$region %>%
  as.character

data_top$group1_sd_df <- square_x_df(data_top, "group1_sd", "group1_df")
data_top$group2_sd_df <- square_x_df(data_top, "group2_sd", "group2_df")
data_top$group3_sd_df <- square_x_df(data_top, "group3_sd", "group3_df")
data_top$ctrl_sd_df <- square_x_df(data_top, "ctrl_sd", "ctrl_df")

# Creating a temp data frame to calculate the ratio1 sd
temp_df <- NULL 
sd_sum <- c("group1_sd_df", "ctrl_sd_df", "group1_df", "ctrl_df")
temp_df <- by_variable(data_top, sd_sum, FUN = sum, VAR = "region") %>% as.data.frame

group1_pooled_sd <- temp_df$group1_sd_df / temp_df$group1_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
regional <- cbind(regional, group1_pooled_sd)

# Creating a temp data frame to calculate the DIS2.GOL sd
temp_df <- NULL 
sd_sum <- c("group2_sd_df", "ctrl_sd_df", "group2_df", "ctrl_df")
temp_df <- by_variable(data_top, sd_sum, FUN = sum, VAR = "region") %>% as.data.frame

group2_pooled_sd <- temp_df$group2_sd_df / temp_df$group2_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
regional <- cbind(regional, group2_pooled_sd)

# Creating a temp data frame to calculate the DIS3.GOL sd
temp_df <- NULL 
sd_sum <- c("group3_sd_df", "ctrl_sd_df", "group3_df", "ctrl_df")
temp_df <- by_variable(data_top, sd_sum, FUN = sum, VAR = "region") %>% as.data.frame

group3_pooled_sd <- temp_df$group3_sd_df / temp_df$group3_df
ctrl_pooled_sd <- temp_df$ctrl_sd_df / temp_df$ctrl_df
regional <- cbind(regional, group3_pooled_sd, ctrl_pooled_sd)

# Grouped relative protein abundance ratio
regional$DIS1.GOL <- by_variable(normval, "DIS1.GOL") %>% as.numeric
regional$DIS2.GOL <- by_variable(normval, "DIS2.GOL") %>% as.numeric
regional$DIS3.GOL <- by_variable(normval, "DIS3.GOL") %>% as.numeric

# Standard deviation relative protein abundance ratio
normval$DIS1.GOL_sd <- ratio_sd(normval, "DIS1.GOL", "group1_sd", "group1_med", 
                           "ctrl_sd", "ctrl_med")
normval$DIS2.GOL_sd <- ratio_sd(normval, "DIS2.GOL", "group2_sd", "group2_med", 
                           "ctrl_sd", "ctrl_med")
normval$DIS3.GOL_sd <- ratio_sd(normval, "DIS3.GOL", "group3_sd", "group3_med", 
                           "ctrl_sd", "ctrl_med")

# Standard deviation grouped relative protein abundance ratio1
normval$DIS1.GOL_df <- normval$group1_df + normval$ctrl_df
normval$DIS1.GOL_sd_df <- square_x_df(normval, "DIS1.GOL_sd", "DIS1.GOL_df")

# Standard deviation grouped relative protein abundance DIS2.GOL
normval$DIS2.GOL_df <- normval$group2_df + normval$ctrl_df
normval$DIS2.GOL_sd_df <- square_x_df(normval, "DIS2.GOL_sd", "DIS2.GOL_df")

# Standard deviation grouped relative protein abundance DIS3.GOL
normval$DIS3.GOL_df <- normval$group2_df + normval$ctrl_df
normval$DIS3.GOL_sd_df <- square_x_df(normval, "DIS3.GOL_sd", "DIS3.GOL_df")

# ratio1 sd_calc
sd_calc <- c("DIS1.GOL_df", "DIS1.GOL_sd_df")
temp_df <- by_variable(normval, sd_calc, FUN = sum, VAR = "region") %>% as.data.frame

DIS1.GOL_sd <- temp_df$DIS1.GOL_sd_df / temp_df$DIS1.GOL_df
regional <- cbind(regional, DIS1.GOL_sd)

# DIS2.GOL sd_calc
sd_calc <- c("DIS2.GOL_df", "DIS2.GOL_sd_df")
temp_df <- by_variable(normval, sd_calc, FUN = sum, VAR = "region") %>% as.data.frame

DIS2.GOL_sd <- temp_df$DIS2.GOL_sd_df / temp_df$DIS2.GOL_df
regional <- cbind(regional, DIS2.GOL_sd)

# DIS3.GOL sd_calc
sd_calc <- c("DIS3.GOL_df", "DIS3.GOL_sd_df")
temp_df <- by_variable(normval, sd_calc, FUN = sum, VAR = "region") %>% as.data.frame

DIS3.GOL_sd <- temp_df$DIS3.GOL_sd_df / temp_df$DIS3.GOL_df
regional <- cbind(regional, DIS3.GOL_sd)

# Statistics -------------------------------------------------------------------

log_norm <- log(normval[, c(group1, group2, group3, ctrl)])
colnames(log_norm) <- paste(colnames(log_norm), sep = "_", "log")
normval <- cbind(normval, log_norm)

log_cols <- grep("log", colnames(normval))
group1_log_cols <- grep("DIS_1_norm_log", colnames(normval))
group2_log_cols <- grep("DIS_2_norm_log", colnames(normval))
group3_log_cols <- grep("DIS_3_norm_log", colnames(normval))
ctrl_log_cols <- grep("GOL_norm_log", colnames(normval))


# Problem if different number of samples be in each group
length(group1_log_cols)
length(group2_log_cols)
length(group3_log_cols)
length(ctrl_log_cols)

# add 12 dummy columns to DIS_2
normval$dummy1_DIS_2_norm_log <- NA
normval$dummy2_DIS_2_norm_log <- NA
normval$dummy3_DIS_2_norm_log <- NA
normval$dummy4_DIS_2_norm_log <- NA
normval$dummy5_DIS_2_norm_log <- NA
normval$dummy6_DIS_2_norm_log <- NA
# normval$dummy7_DIS_2_norm_log <- NA
# normval$dummy8_DIS_2_norm_log <- NA
# normval$dummy9_DIS_2_norm_log <- NA
# normval$dummy10_DIS_2_norm_log <- NA
# normval$dummy11_DIS_2_norm_log <- NA
# normval$dummy12_DIS_2_norm_log <- NA

# add 9 dummy columns to DIS_3
normval$dummy1_DIS_3_norm_log <- NA
normval$dummy2_DIS_3_norm_log <- NA
normval$dummy3_DIS_3_norm_log <- NA
normval$dummy4_DIS_3_norm_log <- NA
normval$dummy5_DIS_3_norm_log <- NA
normval$dummy6_DIS_3_norm_log <- NA
normval$dummy7_DIS_3_norm_log <- NA
normval$dummy8_DIS_3_norm_log <- NA
normval$dummy9_DIS_3_norm_log <- NA
# normval$dummy10_DIS_3_norm_log <- NA
# normval$dummy11_DIS_3_norm_log <- NA
# normval$dummy12_DIS_3_norm_log <- NA
# normval$dummy13_DIS_3_norm_log <- NA
# normval$dummy14_DIS_3_norm_log <- NA
# normval$dummy15_DIS_3_norm_log <- NA
# normval$dummy16_DIS_3_norm_log <- NA
# normval$dummy17_DIS_3_norm_log <- NA
# normval$dummy18_DIS_3_norm_log <- NA

# add 10 dummy columns to GOL
normval$dummy1_GOL_norm_log <- NA
normval$dummy2_GOL_norm_log <- NA
normval$dummy3_GOL_norm_log <- NA
normval$dummy4_GOL_norm_log <- NA
normval$dummy5_GOL_norm_log <- NA
normval$dummy6_GOL_norm_log <- NA
normval$dummy7_GOL_norm_log <- NA
normval$dummy8_GOL_norm_log <- NA
normval$dummy9_GOL_norm_log <- NA
normval$dummy10_GOL_norm_log <- NA
# normval$dummy11_GOL_norm_log <- NA
# normval$dummy12_GOL_norm_log <- NA
# normval$dummy13_GOL_norm_log <- NA
# normval$dummy14_GOL_norm_log <- NA
# normval$dummy15_GOL_norm_log <- NA
# normval$dummy16_GOL_norm_log <- NA
# normval$dummy17_GOL_norm_log <- NA
# normval$dummy18_GOL_norm_log <- NA
# normval$dummy19_GOL_norm_log <- NA
# normval$dummy20_GOL_norm_log <- NA

group1_log_cols <- grep("DIS_1_norm_log", colnames(normval))
group2_log_cols <- grep("DIS_2_norm_log", colnames(normval))
group3_log_cols <- grep("DIS_3_norm_log", colnames(normval))
ctrl_log_cols <- grep("GOL_norm_log", colnames(normval))

stacked <- data.frame(normval[, "region"], 
                      stack(normval[, group1_log_cols]),
                      stack(normval[, group2_log_cols]),
                      stack(normval[, group3_log_cols]),
                      stack(normval[, ctrl_log_cols]))


colnames(stacked) <- c("region", "group1_log", "samp1", 
                       "group2_log", "samp2","group3_log", "samp3",
                       "ctrl_log", "ctrl")

# t-test for group1
p_vals1 <- pval_ttest(stacked, "group1_log", "ctrl_log", col = "region")
colnames(p_vals1) <- c("region", "pvalue1")

p_vals1$region <-
  p_vals1$region %>% 
  as.character

regional <- full_join(regional, p_vals1, by = "region")

# t-test for group2
p_vals2 <- pval_ttest(stacked, "group2_log", "ctrl_log", col = "region")
colnames(p_vals2) <- c("region", "pvalue2")

p_vals2$region <-
  p_vals2$region %>% 
  as.character

regional <- full_join(regional, p_vals2, by = "region")

# t-test for group3
p_vals3 <- pval_ttest(stacked, "group3_log", "ctrl_log", col = "region")
colnames(p_vals3) <- c("region", "pvalue3")

p_vals3$region <-
  p_vals3$region %>% 
  as.character

regional <- full_join(regional, p_vals3, by = "region")

# Minimum number of peptides for each protein group ----------------------------
# min_pep <- 5 
regional$peptides <- table(normval$region)
# regional <- filter(regional, regional$peptides >= min_pep)

write.csv(regional, "Titin regional abundance ratios with outlier subPEVK.csv")
write.csv(normval, "Titin peptide abundance ratios with outlier subPEVK.csv")

# ## Mean ratio for each domain
# 
# # normval <- normval[which(!is.na(normval$DIS1.GOL) | !is.na(normval$DIS2.GOL) |
# #                      !is.na(normval$DIS3.GOL)), ]
# 
# for (i in 1:nrow(dom_linked)) {
#   if (i == 1) {
#     avrat <- data.frame(domain = character(), start = numeric(), end = numeric(),
#                         overlaps = character(), AVG.DIS1 = numeric(), AVG.DIS2 = numeric(),
#                         AVG.DIS3 = numeric(), AVG.Control = numeric(),
#                         "DIS1.GOL" = numeric(), "DIS2.GOL" = numeric(),
#                         "DIS3.GOL" = numeric(), peptides = numeric())
#   }
#   dommean1 <- mean(as.numeric(normval[grep(dom_linked$termdomain[i], normval$domain), "DIS1.GOL"]), na.rm = TRUE)
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
#                      AVG.Control = ctrlmean, "DIS1.GOL" = dommean1,
#                      "DIS2.GOL" = dommean2, "DIS3.GOL" = dommean3,
#                      peptides = length(grep(dom_linked$termdomain[i], normval$domain)))
#   # if (is.nan(temp$KO.Control)) {
#   #   temp$KO.Control <- NA
#   #   }
#   avrat <- full_join(avrat, temp)
# }
# 
# write.csv(avrat, "Titin Domain DIS1 2 3 Ratio 8E6 Threshold.csv")
