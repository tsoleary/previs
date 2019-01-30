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
