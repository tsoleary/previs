# D3 Leucine turnover modeling functions ---------------------------------------

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

# make_pools function 
make_pools <- function(pool_names, syn, deg_old, deg_new, t0_abun){
  
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
  
  return(df)
  
}

# function for the distribution of all the isotopes within a pool 
pool_iso_dist <- function (pool) {
  
  isos <- paste0("M_", seq(((as.numeric(gsub("^D3_([0-9]+)_.*", 
                                             "\\1", pool)))*3), 
                           by = 1, length = 9))
  
  temp <- NULL
  df <- NULL
  for (i in 1:length(isos)){
    temp <- mod[, pool] * nat_iso[i, "per_total"]
    df <- cbind(df, temp)
  }
  
  
  colnames(df) <- paste(str_replace(pool, "_pool", ""), isos, sep = "_")
  
  return(df)
}

# function to get make an individual isotope distribution for each pool
all_isos <- function(pool_cols){
  
  temp <- NULL
  df <- NULL
  
  for (pools in pool_cols){
    temp <- pool_iso_dist(pools)
    df <- cbind(df, temp)
  }
  return(df)
}

