# Functions used in Proteomic analysis for data imported from PD 2.2 -----------

# Function applied to each peptide for a group of data with median as default
as_group <- function (dat, col, FUN = median){
  list <- NULL
  for (i in 1:nrow(dat)){
    temp <- FUN(as.numeric(dat[i, col], na.rm = TRUE))
    list <- c(list, temp)
  }
  return(list)
}

# Function to count degrees of freedom in a group, defaults without duplicates
count_df <- function (dat, col, dup = 1){
  list <- NULL
  for (i in 1:nrow(dat)){
    temp <- length(which(!is.na(dat[i, col])))
    temp <- temp / dup
    if(temp >= 1){
      temp <- temp - 1
    }
    list <- c(list, temp)
  }
  return(list)
}

# Relative abundance of each protein using top ionizers
protein_group <- function (dat, groups, FUN = mean){
  tab <- NULL
  for (col in groups){
    temp <- tapply(dat[, col],
                   dat$Master.Protein.Accessions,
                   FUN,
                   na.rm = TRUE)
    tab <- cbind(tab, temp)
  }
  colnames(tab) <- groups
  return(tab)
}

# Standard deviation of grouped relative abundance using top ionizers
square_x_df <- function (dat, group_sd, group_df){
  (dat[, group_sd])^2 * (dat[, group_df])
}

# Taylor Expansion to get the sd of the ratio of two means
ratio_sd <- function (dat, ratio, group1_sd, group1_med, ctrl_sd, ctrl_med){
  ratio_sd_pep <- NULL
  for (i in 1:nrow(dat)){
    temp1 <- (dat[i, group1_sd] / dat[i, group1_med])^2
    temp2 <- (dat[i, ctrl_sd] / dat[i, ctrl_med])^2
    result <- dat[i, ratio] * sqrt(temp1 + temp2)
    ratio_sd_pep <- c(ratio_sd_pep, result)
  }
  return(ratio_sd_pep)
}

# t-test function
pval_ttest <- function (dat, group, ctrl){
  Master.Protein.Accessions <- NULL
  p_value <- NULL
  for (pro in unique(dat$Master.Protein.Accessions)){
    temp <- filter(dat, dat$Master.Protein.Accessions == pro)
    pval_temp <- t.test(temp[, group], temp[, ctrl])$p.value
    Master.Protein.Accessions <- c(Master.Protein.Accessions, pro)
    p_value <- c(p_value, pval_temp)
  }
  return(as.data.frame(cbind(Master.Protein.Accessions, p_value)))
}


# Convert Master.Protein.Accession to gene symbol
mpa_to_gene <- function (dat, gene_dat){
  dat$gene <- dat$Master.Protein.Accessions
  for (i in 1:nrow(dat)){
    temp <- which(dat$gene[i] == gene_dat$Accession, TRUE)
    if (length(temp) == 1){
      dat$gene <- gsub(dat$gene[i], gene_dat$Gen[temp], dat$gene)
    }
  }
  return(dat$gene)
}

# Remove outliers in peptide data
rm_outliers <- function (dat, ratio){
  
  sd_ratio_temp_df <- protein_group(dat, ratio, FUN = sd) %>% 
    as.data.frame %>% 
    rownames_to_column("Master.Protein.Accessions") %>%
    'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))
  
  pro_temp <- full_join(protein, sd_ratio_temp_df, 
                        by = "Master.Protein.Accessions")
  
  pro_temp$max_ratio <- pro_temp$ratio + 2 * pro_temp$sd_ratio
  pro_temp$min_ratio <- pro_temp$ratio - 2 * pro_temp$sd_ratio
  
  data_rm_out <- NULL
  
  for (pro in unique(dat$Master.Protein.Accessions)){
    temp <- filter(dat, dat$Master.Protein.Accessions == pro)
    
    rm_high <- which(temp$ratio > pro_temp$max_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)])
    
    rm_low <- which(temp$ratio < pro_temp$min_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)]) 
    
    rm <- c(rm_high, rm_low)
    temp_rm <- temp[-rm, ]
    data_rm_out <- bind_rows(data_rm_out, temp_rm)
  }
  return(data_rm_out)
}

# Abundance ratio function
abun_ratio <- function (dat, group, ctrl = "ctrl_med"){
  dat[, group] / dat[, ctrl]
}





