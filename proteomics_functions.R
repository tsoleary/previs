# proteomic analysis top 3 abun & pairwise ratios ------------------------------

# function applied to each peptide for a group of data with median as default
by_group <- function (dat, col, FUN = median){
  list <- NULL
  for (i in 1:nrow(dat)) {
    temp <- FUN(as.numeric(dat[i, col]), na.rm = TRUE)
    list <- c(list, temp)
  }
  return(list)
}

# count the degrees of freedom in a group, defaults without duplicates
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

# relative abundance of each protein using top ionizers
by_protein <- function (dat, groups, FUN = mean){
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

# standard deviation of grouped relative abundance using top ionizers
square_x_df <- function (dat, group_sd, group_df){
  (dat[, group_sd])^2 * (dat[, group_df])
}

# Taylor expansion to get the sd of the ratio of two means
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
pval_ttest <- function (dat, group, ctrl, col = "Master.Protein.Accessions"){
  name <- NULL
  p_value <- NULL
  for (pro in unique(dat[, col])) {
    temp <- dplyr::filter(dat, dat[, col] == pro)
    pval_temp <- tryCatch(t.test(temp[, group], temp[, ctrl])$p.value, 
                          error=function(err) NA)
    name <- c(name, pro)
    p_value <- c(p_value, pval_temp)
  }
  return(as.data.frame(cbind(name, p_value)))
}

# convert Master.Protein.Accession to gene symbol
mpa_to_gene <- function (dat, gene_dat){
  dat$gene <- dat$Master.Protein.Accessions
  for (i in 1:nrow(dat)){
    temp <- which(dat$gene[i] == gene_dat$Accession, TRUE)
    if (length(temp) == 1){
      dat$gene <- gsub(dat$gene[i], gene_dat$Gene[temp], dat$gene)
    }
  }
  return(dat$gene)
}

# remove outliers in peptide data
rm_outliers <- function (dat, pro_df, ratio, mult = 2){
  
  sd_ratio_temp_df <- by_protein(dat, ratio, FUN = sd) %>%
    as.data.frame %>%
    rownames_to_column("Master.Protein.Accessions") %>%
    'colnames<-' (c("Master.Protein.Accessions", "sd_ratios"))
  
  pro_temp <- dplyr::full_join(pro_df, sd_ratio_temp_df,
                               by = "Master.Protein.Accessions")
  
  pro_temp$max_ratio <- pro_temp$ratio + mult * pro_temp$sd_ratio
  pro_temp$min_ratio <- pro_temp$ratio - mult * pro_temp$sd_ratio
  
  data_rm_out <- NULL
  
  for (pro in unique(dat$Master.Protein.Accessions)){
    temp <- dplyr::filter(dat, dat$Master.Protein.Accessions == pro)
    
    rm_high <- which(temp[, ratio] > pro_temp$max_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)])
    
    rm_low <- which(temp[, ratio] < pro_temp$min_ratio[which(
      pro_temp$Master.Protein.Accessions == pro)])
    
    rm <- c(rm_high, rm_low)
    if (length(rm) > 0){
      temp_rm <- temp[-rm, ]
      data_rm_out <- dplyr::bind_rows(data_rm_out, temp_rm)
    } else {
      data_rm_out <- dplyr::bind_rows(data_rm_out, temp)
    }
  }
  return(data_rm_out)
}

# Abundance ratio function
abun_ratio <- function (dat, group, ctrl = "ctrl_med"){
  dat[, group] / dat[, ctrl]
}

# whole proteome ---------------------------------------------------------------
proteome <- function(file_name, organism,
                     wt_samps = "all", group = FALSE, group_names = NULL,
                     top_pep = 3, min_pep = 3,
                     norm = FALSE, norm_method = NULL, norm_pro = NULL,
                     csv = TRUE){
  
  # import data frame
  dat <- read.csv(file_name)
  
  # data frame format check
  if (colnames(dat)[1] != "Annotated.Sequence" &
      colnames(dat)[2] != "Modifications" &
      colnames(dat)[3] != "Master.Protein.Accessions")
    stop ("invalid data frame: must have Annotated.Sequence, Modifications,
          Master.Protein.Accessions")
  
  # import gene accession
  if (organism == "mouse") {
    gene_df <- read.csv(
      "C:/Users/PrevBeast/Documents/Fasta Files/Mouse/mouse_fasta_accession_gene.csv")
  }
  if (organism == "human"){
    gene_df <- read.csv(
      "C:/Users/PrevBeast/Documents/Fasta Files/Human/human_fasta_accession_gene.csv")
  }
  if (organism != "mouse" & organism != "human")
    stop ("only human or mouse organism currently supported")
  
  
  # determine the top ionizers and sort by them taking only the top_pep
  if (wt_samps == "all"){
    wt_grp <- 4:ncol(dat)
  } else {
    if (sum(grepl(wt_samps, colnames(dat))) < 1)
      stop ("wt_samps must be a pattern that is present in the sample
            (column) names that are being used to weight the top ionizers")
    wt_grp <- grep(wt_samps, colnames(dat))
  }
  
  dat$wt_grp <- by_group(dat, wt_grp)
  
  df <- tbl_df(dat) %>%
    group_by(Master.Protein.Accessions) %>%
    top_n(n = top_pep, wt = wt_grp)
  
  
  # do normalization if indicated --------
  if (norm == TRUE){
    
    if (norm_method == "sum_total") {
      norm_value <- colSums(df[, wt_grp], na.rm = TRUE)
      
      raw_abun_mat <- as.matrix(df[, wt_grp])
      
      norm_abun <- t(t(raw_abun_mat)/norm_value)
      colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
      norm_test <- as.data.frame(norm_abun)
      df <- tibble::as_tibble(df)
      df <- cbind(df, norm_test)
      
    }
    
    if (norm_method == "protein"){
      
      norm_pep <- subset(df, df$Master.Protein.Accessions == norm_pro)
      numeric_cols <- which(sapply(norm_pep, is.numeric) == TRUE)
      raw_abun <- numeric_cols[-length(numeric_cols)]
      norm_value <- sapply(norm_pep[, raw_abun], mean, na.rm = TRUE)
      raw_abun_mat <- as.matrix(df[, raw_abun])
      norm_abun <- t(t(raw_abun_mat)/norm_value)
      colnames(norm_abun) <- paste(colnames(norm_abun), sep = "_", "norm")
      norm_test <- as.data.frame(norm_abun)
      df <- tibble::as_tibble(df)
      df <- cbind(df, norm_test)
    }
  }
  
  df <- as.data.frame(df)
  
  if (group == TRUE){
    for (i in 1:length(group_names)){
      
      if (norm == TRUE){
        pat <- paste0(group_names[i], "_norm")
      } else {
        pat <- group_names[i]
      }
      samp_col_group <- grep(pat, colnames(df))
      df[, paste0(group_names[i], "_med")] <- by_group(df, samp_col_group)
    }
    
    pro_df <- by_protein(df, colnames(df)[grep("_med", colnames(df))]) %>%
      as.data.frame %>%
      rownames_to_column("Master.Protein.Accessions")
    
    colnames(pro_df) <- sub("_med", "", colnames(pro_df))
    
  } 
  
  
  if (group == FALSE) {
    if (norm == FALSE) {
    all <- colnames(dat)[-c(1,2,3,length(dat))]
    
    pro_df <- by_protein(df, all) %>%
      as.data.frame %>%
      rownames_to_column("Master.Protein.Accessions")
    }
    if (norm == TRUE){
      pro_df <- by_protein(df, colnames(df)[grep("_norm", colnames(df))]) %>%
        as.data.frame %>%
        rownames_to_column("Master.Protein.Accessions")
      
      colnames(pro_df) <- sub("_norm", "", colnames(pro_df))
    }
    
  }
 
  
  pro_df$gene <- mpa_to_gene(pro_df, gene_df)
  
  pro_df$peptides <- table(dat$Master.Protein.Accessions)
  pro_df <- dplyr::filter(pro_df, pro_df$peptides >= min_pep)
  
  col_order <- c("gene", "Master.Protein.Accessions", "peptides",
                 colnames(pro_df)[2:(length(colnames(pro_df))-2)])
  
  pro_df <- pro_df[, col_order]
  
  colnames(pro_df)[1] <- "Gene"
  colnames(pro_df)[2] <- "Accession"
  colnames(pro_df)[3] <- "# of peptides"
  
  pro_df <- dplyr::filter(pro_df, Accession != "")
  
  file_output <- paste0(gsub(".csv", "", file_name), "_r_",
                        format(Sys.time(), "%d %b %Y %H:%M:%S"),
                        ".csv")
  
  if (csv == TRUE){
    write.csv(pro_df, file_output, row.names = FALSE)
    print(paste("result file", file_output ,"is located in", getwd()))
  }
  pro_df <- tibble::as_tibble(pro_df)
  return(pro_df)
}


# proteins changing over time --------------------------------------------------

# get the slope, intercept r_sq, and pval of each with do(lin_fit(.))
lin_fit <- function(dat){
  the_fit <- lm(dat$abundance ~ dat$week, dat)
  p_val <- anova(the_fit)$'Pr(>F)'[1]
  slo_int <- data.frame(t(coef(the_fit)))
  r_sq <- summary(the_fit)$r.squared
  result <- cbind(slo_int, r_sq, p_val)
  colnames(result) <- c("intercept", "slope", "r_squared", "p_value")
  return(result)
}

# get the regression eqn, r_sq, & pval for plotting from a lm_object 
lm_eqn <- function(lm_object) {
  eq <-
    substitute(
      italic(y) == m ~ italic(x) + b,
      list(
        b = as.character(signif(coef(lm_object)[1], digits = 2)),
        m = as.character(signif(coef(lm_object)[2], digits = 2))
      )
    )
  r <- 
    substitute(
      italic(r) ^ 2 ~ "=" ~ r2,
      list(
        r2 = as.character(signif(summary(lm_object)$r.squared, digits = 3))
      )
    )
  p <-  
    substitute(
      italic("p-val") ~ "=" ~ pval,
      list(
        pval = as.character(signif(anova(lm_object)$'Pr(>F)'[1], digits = 3))
      )
    )
  result <- c(as.character(as.expression(eq)), as.character(as.expression(r)),
              as.character(as.expression(p)))
  return(result)
}

# plot_pro function
plot_pro <- function(dat, g_title, FUN = geom_jitter, x_pos = 9, y1_pos = 0.83,
                     y2_pos = 0.30){
  
  df_l <- filter(dat, leg == "L")
  df_r <- filter(dat, leg == "R")
  
  lm_l <- lm(abundance ~ week, df_l)
  lm_r <- lm(abundance ~ week, df_r)
  
  eqn_l <- lm_eqn(lm_l)
  eqn_r <- lm_eqn(lm_r)
  
  g <- ggplot(dat, aes(x = week, y = abundance)) +
    FUN(mapping = aes(x = week, y = abundance, fill = leg), 
        alpha = 0.5, size = 3, pch = 21,  color = "black", width = 0.05) +
    labs(title = g_title, x = "Week", y = "Raw Abundance", fill = "Leg") +
    expand_limits(x = 0, y = 0) +
    theme_classic() +  
    expand_limits(x = 0, y = 0) +
    geom_smooth(mapping = aes(color = leg), method = 'lm', se = FALSE, 
                size = 1.1, show.legend = FALSE, linetype = "dotted") +
    annotate("text", x = x_pos, y = (y1_pos)*(max(dat$abundance)), 
             label = eqn_l[1], parse = TRUE, color = "#F98B86") +
    annotate("text", x = x_pos, y = (y1_pos - 0.05)*(max(dat$abundance)), 
             label = eqn_l[2], parse = TRUE, color = "#F98B86") +
    annotate("text", x = x_pos, y = (y1_pos - 0.11)*(max(dat$abundance)), 
             label = eqn_l[3], parse = TRUE, color = "#F98B86") +
    annotate("text", x = x_pos, y = (y2_pos)*(max(dat$abundance)), 
             label = eqn_r[1], parse = TRUE, color = "#53D3D7") +
    annotate("text", x = x_pos, y = (y2_pos - 0.05)*(max(dat$abundance)), 
             label = eqn_r[2], parse = TRUE, color = "#53D3D7") +
    annotate("text", x = x_pos, y = (y2_pos - 0.11)*(max(dat$abundance)), 
             label = eqn_r[3], parse = TRUE, color = "#53D3D7") +
    coord_cartesian(xlim = c(1, 8), clip = 'off') +
    theme(plot.margin = unit(c(1, 5, 1, 1), "lines"))
  return(g)
}

# convert an individual accession to gene for each graph 
indiv_mpa_to_gene <- function (Acc_pro, gene_dat){
  temp <- which(Acc_pro == gene_dat$Master.Protein.Accessions, TRUE)
  gene <- gsub(Acc_pro, gene_dat$gene[temp], Acc_pro)
  return(gene)
}

# D3_Leu labeled turnover study ------------------------------------------------

# calculates the distribution of newly synthesized D3-Leu peptides
iso_dist <- function (pep, p = 0.5){
  
  L <- str_count(pep, pattern = "L")
  isos <- 0:L
  dis <- NULL
  
  for (i in isos){
    temp <- (factorial(L)/(factorial(L-i)*factorial(i)))*(p^i)*(1-p)^(L-i)
    dis <- c(dis, temp)
  }
  return(dis)
}

# isotope distribution for a peptide - data frame output
pep_iso <- function (pep, max_iso = 9, charge = 1){
  pep_length <- nchar(pep)
  pep_bond <- pep_length - 1
  
  aa_sep <- NULL
  for (i in 1:pep_length){
    aa <- substr(pep, i, i)
    aa_sep <- c(aa_sep, aa)
  }
  
  # reference data frame with molecular formulas of all amino acids
  aa_form <- data.frame("letter" = c("A", "R", "N", "D", "C", "Q", "E", "G", 
                                     "H", "I", "L", "K", "M", "F", "P", "S", 
                                     "T", "W", "Y", "V"),
                        "molecular_formula" = c("C3H7NO2", "C6H14N4O2", 
                                                "C4H8N2O3", "C4H7NO4", 
                                                "C3H7NO2S", "C5H10N2O3", 
                                                "C5H9NO4", "C2H5NO2", 
                                                "C6H9N3O2", "C6H13NO2", 
                                                "C6H13NO2", "C6H14N2O2", 
                                                "C5H11NO2S", "C9H11NO2", 
                                                "C5H9NO2", "C3H7NO3", 
                                                "C4H9NO3", "C11H12N2O2", 
                                                "C9H11NO3", "C5H11NO2"))
  
  mols <- NULL
  for (j in 1:length(aa_sep)){
    row <- grep(aa_sep[j], aa_form$letter)
    mol <- as.character(aa_form[row, "molecular_formula"])
    mols <- paste0(mols, mol)
  }
  
  pep_tot <- getMolecule(mols)$formula
  
  element <- function(formula){
    # pattern to match the initial element assumes element starts with an upper 
    # case and optional lower case followed by zero or more digits.
    first <- "^([[:upper:]][[:lower:]]?)([0-9]*).*"
    # inverse of above to remove the initial element
    last <- "^[[:upper:]][[:lower:]]?[0-9]*(.*)"
    result <- list()
    extract <- formula
    # repeat as long as there is data
    while ((start <- nchar(extract)) > 0){
      chem <- sub(first, '\\1 \\2', extract)
      extract <- sub(last, '\\1', extract)
      # if the number of characters is the same, then there was an error
      if (nchar(extract) == start){
        warning("Invalid formula:", formula)
        return(NULL)
      }
      # append to the list
      result[[length(result) + 1L]] <- strsplit(chem, ' ')[[1]]
    }
    result
  }
  
  ans <- unlist(element(pep_tot))
  
  if (ans[length(ans)] == "S"){
    ans <- c(ans, "0")
  }
  
  elem <- ans[rep(seq(from = 1, to = length(ans), by = 2), 1)]
  num <- ans[rep(seq(from = 2, to = length(ans), by = 2), 1)]
  
  df <- as.data.frame(cbind(elem, num))
  df$num <- as.numeric(as.character(df$num))
  
  # correcting for loss of water in the peptide bond
  df[elem == "H", "num"] <- (df[elem == "H", "num"] - (2 * pep_bond) + charge)
  df[elem == "O", "num"] <- (df[elem == "O", "num"] - 1 * pep_bond)
  
  # convert from data frame back to molecular fomula string
  df_mat <- as.matrix(df)
  df_list <- NULL
  for (k in seq(1:nrow(df_mat))){
    df_list <- c(df_list, df_mat[k, ])
  }
  pep_mol <- gsub(" ", "", 
                  paste(as.character(df_list), sep = "", collapse = ""))
  
  # get isotope on final molecule
  pep_molecule <- getMolecule(pep_mol, z = charge)
  pep_dist <- as.data.frame(getIsotope(pep_molecule, seq(1, max_iso)), 
                            row.names = c("m_z", "per_total"))
  colnames(pep_dist) <- paste0("M_", 0:(max_iso-1))
  pep_dist <- as.data.frame(t(pep_dist))
  pep_dist$m_z <- (pep_dist$m_z / charge)
  pep_dist <- round(pep_dist, 3)
  
  return(pep_dist)
  
}

# deg_syn function
deg_syn <- function (dat, deg, syn, initial, D3_pep, D3_i){
  
  temp <- NULL
  list <- NULL
  D3_syn <- syn * D3_pep[D3_i + 1]
  
  for (i in 1:((nrow(dat) - 1))){
    if (is.null(temp) == TRUE){
      temp <- initial
      list<- c(temp)
    }
    temp <- temp - (temp * deg) + D3_syn
    list <- c(list, temp)
  }
  return(list)
}

# make_pools function 
make_pools <- function(dat, pool_names, syn, deg_old, deg_new, t0_abun, D3_pep){
  
  temp <- NULL
  df <- NULL
  
  for (pool in pool_names){
    if (grepl("old", pool) == TRUE){
      f_syn <- 0
      f_deg <- deg_old
      f_initial <- t0_abun
    } else {
      f_syn <- syn
      f_deg <- deg_new
      f_initial <- 0
    }
    
    temp <- deg_syn(dat, f_deg, f_syn, f_initial, D3_pep,
                    D3_i = as.numeric(gsub("^D3_([0-9]+)_.*", "\\1", pool)))
    df <- cbind(df, temp)
  }
  colnames(df) <- pool_names
  return(df)
}

# function for the distribution of all the isotopes within a pool 
pool_iso_dist <- function (dat, pool, peptide) {
  
  isos <- paste0("M_", seq(((as.numeric(gsub("^D3_([0-9]+)_.*", 
                                             "\\1", pool)))*3), 
                           by = 1, length = 9))
  
  nat_iso <- pep_iso(peptide)
  
  temp <- NULL
  df <- NULL
  for (i in 1:length(isos)){
    temp <- dat[, pool] * nat_iso[i, "per_total"]
    df <- cbind(df, temp)
  }
  
  colnames(df) <- paste(str_replace(pool, "_pool", ""), isos, sep = "_")
  
  return(df)
}

# function to get make an individual isotope distribution for each pool
all_isos <- function(dat, pool_cols, peptide){
  
  temp <- NULL
  df <- NULL
  
  for (pools in pool_cols){
    temp <- pool_iso_dist(dat, pools, peptide)
    df <- cbind(df, temp)
  }
  return(df)
}

# coombine the individual isotopomers from different pools back to a total 
combine_iso <- function (dat){
  
  # removes the M_ from each column name?
  colnames(dat) <- gsub("^.*(M_[0-9]+)", "\\1", colnames(dat))
  
  # combine all columns that have the same name by their 
  df <- sapply(split.default(dat, colnames(dat)), 
               rowSums, na.rm = TRUE)
  
  # add the word total to every colname with a M_0 pattern
  colnames(df) <- gsub("^.*(M_[0-9]+)", 
                       paste0("\\1", "_total"), colnames(df))  
  
  # remove pook and time from the df
  df <- df[, grep(c("M_"), colnames(df))]
  
  # pading with a 0 to allow sorting of the columns 
  colnames(df) <- str_pad(gsub("M_", "", colnames(df)), 8, pad = "0")
  
  df <- df[, sort(colnames(df))]
  
  colnames(df) <- c(paste0("M_", colnames(df)))
  
  return(df)
  
}

# mega model creating the whole thing so it is one line of code
mega_model <- function (peptide, deg_old, deg_new, syn, t0_abun, per_lab, time){
  # create a data frame with a specified length of time 
  mod <- data.frame("time" = time)
  
  # define the number of D3 Leucines that will be in the peptide and their 
  # proportional syntheis rates based on the percent labeling 
  D3_pep_dist <- iso_dist(peptide, p = per_lab)
  
  # define the pool names to be used as an argument in the make pools function
  pool_names <- c("D3_0_old_pool", "D3_0_new_pool",
                  paste("D3", 1:(length(D3_pep_dist) - 1), "pool", sep = "_"))
  
  # make pools for all 
  pool_df <- make_pools(mod, pool_names, syn, deg_old, deg_new, t0_abun, 
                        D3_pep_dist)
  
  # cbind to the data frame with the time vector
  mod <- cbind(mod, pool_df)
  
  # create a data frame with all isotopes for each pool
  df_all_isos <- all_isos(mod, pool_names, peptide)
  
  # cbind all columns together
  mod <- cbind(mod, df_all_isos)
  
  # sum all the individual isotopes to get the totals
  mod <- cbind(mod, combine_iso(mod))
  
  return(mod)
}

# D2O Turnover functions -------------------------------------------------------

# plot isotopic ratios 
plot_iso_r <- function(dat, g_title, g_subtitle, FUN = geom_jitter, 
                       x_pos = 3.35, y_pos1 = 0.83, y_pos2 = 0.30){
  
  df <- dat %>%
    group_by(leg) %>%
    do(fit = nls(ratio ~ SSasymp(week, yf, y0, log_k), data = .)) %>%
    tidy(fit) %>%
    select(leg, term, estimate) %>%
    spread(term, estimate) %>%
    mutate(k = exp(log_k))
  
  L_k_yf <- k_yf(df, "L")
  R_k_yf <- k_yf(df, "R")
  
  g <- ggplot(dat, aes(x = week, y = ratio)) +
    FUN(mapping = aes(x = week, y = ratio, fill = leg), 
        alpha = 0.5, size = 3, pch = 21,  color = "black", width = 0.05) +
    geom_smooth(method = "nls", formula = y ~ SSasymp(x, yf, y0, log_k),
                data = dat,
                se = FALSE,
                aes(color = leg), show.legend = FALSE) +
    annotate("text", x = x_pos, y = (y_pos1)*(max(dat$ratio)), 
             label = L_k_yf[1], parse = TRUE, color = "#F98B86") +
    annotate("text", x = x_pos, y = (y_pos1 - 0.07)*(max(dat$ratio)), 
             label = L_k_yf[2], parse = TRUE, color = "#F98B86") +
    annotate("text", x = x_pos, y = (y_pos2)*(max(dat$ratio)), 
             label = R_k_yf[1], parse = TRUE, color = "#53D3D7") +
    annotate("text", x = x_pos, y = (y_pos2 - 0.07)*(max(dat$ratio)), 
             label = R_k_yf[2], parse = TRUE, color = "#53D3D7") +
    labs(title = g_title, subtitle = g_subtitle, 
         x = "Week", y = "Isotopic Ratio\nCorrected", # change later!
         fill = "Leg") +
    expand_limits(x = 0, y = 0) +
    theme_classic() +    
    coord_cartesian(xlim = c(0, 3), clip = 'off') +
    theme(plot.margin = unit(c(1, 2.25, 1, 1), "lines"))
  return(g)
  
}

# corrected m0/(m0+m1)
iso_ratio_calc <- function(dat){
  
  ratios <- NULL
  
  for (i in 1:nrow(dat)){
    
    t0_ratio <- mean(dat$t0_r[which(dat$protein == dat$protein[i] & 
                                         dat$peptide == dat$peptide[i] &
                                         dat$sex == dat$sex[i] & 
                                         dat$leg == dat$leg[i] & 
                                         dat$week == 0)], na.rm = TRUE)
    
    temp <- (dat$M_1[i] + dat$M_2[i] + dat$M_3[i] - (t0_ratio * dat$M_0[i])) / 
      (dat$M_0[i] + (dat$M_1[i] + dat$M_2[i] + dat$M_3[i] - 
                       (t0_ratio * dat$M_0[i])))
    
    ratios <- c(ratios, temp)
    
  }
  
  return(ratios)
  
}

# get the k and yf value into an equation to be put on graph 
k_yf <- function(nls_df, leg) {
  k <-
    substitute(
      italic(k) ~ "=" ~ k_val,
      list(
        k_val = as.character(signif(nls_df$k[which(nls_df$leg == leg)], 
                                    digits = 3))
      )
    )
  yf <- 
    substitute(
      italic(yf) ~ "=" ~ yf_val,
      list(
        yf_val = as.character(signif(nls_df$yf[which(nls_df$leg == leg)], 
                                    digits = 2))
      )
    )
  result <- c(as.character(as.expression(k)), as.character(as.expression(yf)))
  return(result)
}

# plot the isotopic distribution changing over time
plot_all_isos <- function(dat, g_title, g_subtitle, FUN = geom_jitter){
  g <- ggplot(dat, aes(x = isotope, y = abundance)) + 
    FUN(mapping = aes(x = isotope, y = abundance, fill = week), 
        alpha = 0.5, size = 3, pch = 21, color = "black", width = 0.01) +
    scale_fill_brewer(palette = "Spectral") +
    labs(title = g_title, subtitle = g_subtitle, 
         x = "Isotopomer", y = "Abundance\n(% Total)",
         fill = "Week") + 
    theme_bw() +
    facet_wrap( ~ dat$leg) +
    theme(strip.background = element_rect(color = "black", fill = "#53D3D7"))
  return(g)
}

# normalize isotope columns
norm_iso <- function(dat){
  isotopes <- colnames(dat)[grep("M_", colnames(dat))]
  df <- NULL
  for(iso in isotopes){
    temp <- dat[iso]/dat["Sum"]
    df <- c(df, temp)
  }
  df <- as.data.frame(df)
  colnames(df) <- paste0(isotopes, "_norm")
  dat <- cbind(dat, df)
}


