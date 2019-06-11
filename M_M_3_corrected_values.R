# Take pinpoint isotopes plot and correct for synthesis and decay calc ---------

library(tidyverse)
library(plotly)
library("Rdisop")

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

aa_form <- read.csv("aa_molecular_formula.csv")

# Functions --------------------------------------------------------------------

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
  pep_mol <- str_replace(paste(as.character(df_list), sep = "", collapse = ""),
                         " ", "")
  pep_mol <- str_replace(pep_mol, " ", "")
  
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
deg_syn <- function (deg, syn, initial, D3_i){
  
  temp <- NULL
  list <- NULL
  D3_syn <- syn * D3_pep_dist[D3_i + 1]
  
  for (i in 1:((nrow(mod) - 1))){
    if (is.null(temp) == TRUE){
      temp <- initial
      list<- c(temp)
    }
    temp <- temp - (temp * deg) + D3_syn
    list <- c(list, temp)
  }
  return(list)
}

# initial conditions -----------------------------------------------------------

# peptide specific
peptide <- "SAMPLLER"
D3_pep_dist <- iso_dist(peptide, p = 0.5)

# natural isotopic distribution of the peptide
nat_iso <- pep_iso(peptide, charge = 1)

# initial total abundance at T0
t0_abun <- 1000

# rates
deg_old <- 0.0500
deg_new <- 0.0500
syn <- 50

# length of time
time <- 0:168

# model data frame -------------------------------------------------------------

mod <- data.frame("time" = time)

# D3_0_old
mod$D3_0_old_pool <- deg_syn(deg_old, syn = 0, initial = t0_abun, D3_i = 0)

# D3_0_new
mod$D3_0_new_pool <- deg_syn(deg_new, syn, initial = 0, D3_i = 0)

# D3_1_new
mod$D3_1_new_pool <- deg_syn(deg_new, syn, initial = 0, D3_i = 1)

# D3_2_new
mod$D3_2_new_pool <- deg_syn(deg_new, syn, initial = 0, D3_i = 2)


pool_names <- c("D3_0_old_pool", "D3_0_new_pool", 
                paste("D3", 1:(length(D3_pep_dist) - 1), "new_pool", sep = "_"))

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
  
  temp <- deg_syn(f_deg, f_syn, f_initial, 
                  D3_i = as.numeric(gsub("^D3_([0-9]+)_.*", "\\1", pool)))
  
  df <- cbind(df, temp)
  
}

colnames(df) <- pool_names


mod <- cbind(mod, df)




# create a function that now makes these columns automatically for each pool

for (i in 1:length(D3_pep_dist)){
  
}


# function for the distribution of all the isotopes within a pool --------------
pool_iso_dist <- function (pool, isos = paste0("M_", 0:9)){
  df <- NULL
  for (i in isos){
    temp <- mod[, pool] * nat_iso[i, "per_total"]
    df <- cbind(df, temp)
  }
  colnames(df) <- paste(pool, isos, sep = "_")
  return(df)
}

df_test <- cbind(mod, pool_iso_dist("D3_0_old_pool"))

pool_cols <- grep("pool", colnames(mod))

for (i in pool_cols){
  
}
