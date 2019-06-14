# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(plotly)
library(Rdisop)

setwd("C:/Users/PrevBeast/Documents/R/Helms")

# Tidy the data ----------------------------------------------------------------

df_messy <- read.csv("All_isotopes_organized_select_reprex.csv")

samples <- colnames(df_messy)[4:length(colnames(df_messy))]

df_tidy <- df_messy %>% 
  gather(sample, abundance, samples) %>%
  separate("sample", c("group", "duplicate", "hrs"), sep = "_")

df_tidy$hrs <- as.numeric(df_tidy$hrs)

# average together the duplicate lines
df_avg_3 <- df_tidy %>%
  group_by(protein, peptide, isotope, group, hrs) %>%
  summarize(abundance = mean(abundance, na.rm = TRUE))

# spreads out to columns by isotope: M_0, M_3, & Sum
df <- as.data.frame(spread(df_avg, "isotope", "abundance"))

# plot by protein, peptide, and isotope ----------------------------------------

pro <- as.character(unique(df$protein))

for (i in 1:length(pro)){
  
  temp <- filter(df, protein == pro[i])
  
  pep <- as.character(unique(temp$peptide))
  
  for (pep_x in pep){
  
    temp_pep <- filter(temp, peptide == pep_x)
    
    cols <- c(grep(c("M_0"), colnames(temp_pep)),
              grep(c("M_3"), colnames(temp_pep)),
              grep(c("Sum"), colnames(temp_pep)))
    
    for (j in cols){
     
      g1 <- ggplot(data = temp_pep) +
              geom_point(mapping = aes(x = hrs, 
                                       y = as.numeric(temp_pep[, colnames(temp_pep)[j]]), 
                                       color = Group), 
                         size = 3, alpha = 0.8) + 
              labs(title = paste(pro[i], pep_x, sep = " -- "), 
                   subtitle = colnames(temp_pep)[j],
                   y = "Abundance", x = "Time (hrs)") +
              theme_classic()
      plot(g1)
    }
  }
}

# figure out how to fit a one phase decay to the M-0 values -- probably need 
# to do the corrected M-0 values
# maybe figure out how to put all three M_0, M_3, and Sum in the same sort of 
# space
# plotly interactive too
# make the rest of it nice

# working on aesthetics
# ggplot(data = temp_pep) +
#   geom_point(mapping = aes(x = hrs,
#                            y = as.numeric(temp_pep[, colnames(temp_pep)[7]]),
#                            fill = sample, alpha = 0.5),
#              shape = 21, size = 3, stroke = 2, alpha = 0.9) +
#   labs(title = paste(pro[5], pep_x, sep = " -- "),
#        subtitle = colnames(temp_pep)[7],
#        y = "Abundance", x = "Time (hrs)") +
#   theme_classic()

# ggsave("test.pdf", dpi = 320) could save the specific files in a .pdf using
# the paste0() as the argument of the ggsave() in the for loop.


# Using plotly package ---------------------------------------------------------
# http://www.rebeccabarter.com/blog/2017-04-20-interactive/

# histone normalization of all isotopes ----------------------------------------

# average histones
df_hist <- df %>%
  group_by(protein, group, hrs) %>%
  summarize(avg_pro = mean(Sum)) %>%
  filter(protein == "Histones")

df <- cbind(df, hist_avg = df_hist$avg_pro)

# normailize all isotopes by histones 
isotopes <- c("M_0", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "Sum")

for (col in isotopes){
  df[, paste(col, "hist_norm", sep = "_")] <- df[, col] / df$hist_avg
}

# Mega Model -------------------------------------------------------------------

aa_for <- read.csv("aa_molecular_formula.csv")

# initial conditions -----------------------------------------------------------

# # peptide specific
# peptide <- "SAMPLLER"
# D3_pep_dist <- iso_dist(peptide, p = 0.5)
# 
# # natural isotopic distribution of the peptide
# nat_iso <- pep_iso(peptide, charge = 1)
# 
# # initial total abundance at T0
# t0_abun <- 1000
# 
# # rates
# deg_old <- 0.0500
# deg_new <- 0.0500
# syn <- 50
# 
# # length of time
# time <- 0:168

# model data frame -------------------------------------------------------------

# mod <- data.frame("time" = time)
# 
# pool_names <- c("D3_0_old_pool", "D3_0_new_pool", 
#                 paste("D3", 1:(length(D3_pep_dist) - 1), "pool", sep = "_"))
# 
# pool_df <- make_pools(pool_names, syn = 50, deg_old = 0.05, 
#                       deg_new = 0.05, t0_abun = 1000)
# 
# mod <- cbind(mod, pool_df)
# 
# mod <- cbind(mod, pool_iso_dist("D3_0_old_pool"))
# 
# pool_cols <- colnames(mod)[grep("pool", colnames(mod))]
# 
# df_all_isos <- all_isos(pool_cols)
# 
# mod <- cbind(mod, df_all_isos)

# total function 
model_turnover <- function (peptide, deg_old, deg_new, syn, t0_abun, per_lab, 
                            time = 0:168){
  
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
    
  # distribution of pools based on the number of labeled leucines
  L <- str_count(peptide, pattern = "L")
  isos <- 0:L
  p <- per_lab
  D3_pep_dist <- NULL
  
  for (i in isos){
    temp <- (factorial(L)/(factorial(L-i)*factorial(i)))*(p^i)*(1-p)^(L-i)
    D3_pep_dist <- c(D3_pep_dist, temp)
  }
  
  # coverting the primary sequence of the peptide into a total molecular formula
  pep_length <- nchar(peptide)
  pep_bond <- pep_length - 1
  
  aa_sep <- NULL
  for (i in 1:pep_length){
    aa <- substr(peptide, i, i)
    aa_sep <- c(aa_sep, aa)
  }
  
  mols <- NULL
  for (j in 1:length(aa_sep)){
    row <- grep(aa_sep[j], aa_form$letter)
    mol <- as.character(aa_form[row, "molecular_formula"])
    mols <- paste0(mols, mol)
  }
  
  pep_tot <- getMolecule(mols)$formula
  # pattern to match the initial element assumes element starts with an upper 
  # case and optional lower case followed by zero or more digits.
  first <- "^([[:upper:]][[:lower:]]?)([0-9]*).*"
  # inverse of above to remove the initial element
  last <- "^[[:upper:]][[:lower:]]?[0-9]*(.*)"
  result <- list()
  extract <- pep_tot
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
  nat_iso <- round(pep_dist, 3)
  
  pool_names <- c("D3_0_old_pool", "D3_0_new_pool", 
                  paste("D3", 1:(length(D3_pep_dist) - 1), "pool", sep = "_"))
  
  mod <- data.frame("time" = time)
  
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
    
    temp <- deg_syn(df, f_deg, f_syn, f_initial, 
                    D3_i = as.numeric(gsub("^D3_([0-9]+)_.*", "\\1", pool)))
    
    deg_syn <- function (df, deg, syn, initial, D3_i){
      
      temp <- NULL
      list <- NULL
      D3_syn <- syn * D3_pep_dist[D3_i + 1]
      
      for (i in 1:((nrow(df) - 1))){
        if (is.null(temp) == TRUE){
          temp <- initial
          list<- c(temp)
        }
        temp <- temp - (temp * deg) + D3_syn
        list <- c(list, temp)
      }
      return(list)
    }
    
    df <- cbind(df, temp)
  }
  
  colnames(df) <- pool_names
  
  mod <- cbind(mod, make_pools(pool_names, syn, deg_old, deg_new, t0_abun))
  
  all_isos <- function(dat, pool_cols){
    
    temp <- NULL
    df <- NULL
    
    for (pools in pool_cols){
      temp <- pool_iso_dist(dat, pools)
      
      pool_iso_dist <- function (df, pool) {
        
        isos <- paste0("M_", seq(((as.numeric(gsub("^D3_([0-9]+)_.*", 
                                                   "\\1", pool)))*3), 
                                 by = 1, length = 9))
        
        temp <- NULL
        df <- NULL
        for (i in 1:length(isos)){
          temp <- df[, pool] * nat_iso[i, "per_total"]
          df <- cbind(df, temp)
        }
        
        
        colnames(df) <- paste(str_replace(pool, "_pool", ""), isos, sep = "_")
        
        return(df)
      }
      
      
      
      
      df <- cbind(df, temp)
    }
    return(df)
  }
  
  
  
  mod <- cbind(mod, all_isos(colnames(mod)[grep("pool", colnames(mod))]))
  
  return(mod)
  
}

x <- model_turnover("SAMPLLER", 0.05, 0.05, 50, 1000, 0.50)



# need to create a function that combines all individual isotopes into one total

mutate(mod, sum = colnames(mod)[grep("M_0", colnames(mod))], sum)

df_test <- mod
colnames(df_test) <- gsub("^.*(M_[0-9]+)", "\\1", colnames(df_test))

mod[, paste("total", iso, sep = "_")] <- mod[, 6] + mod[, 15]

mod_comb <- sapply(split.default(df_test, colnames(df_test)), 
                   rowSums, na.rm = TRUE)


y <- colnames(mod_comb)

# the columns come out of order because it is lexographic indexing not numerical
# need to pad the numbers with a 0

x <- str_pad(gsub("M_", "", y), 2, pad = "0")

sort(x)

# need to figure out how to work backwards to get the M_0 corrected values etc






