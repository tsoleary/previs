#####Set Up#####
getwd()   #See what folder the working directory is set to
setwd("C:/Users/xyz/Documents/R/Previs")   #Set the working directory to the specfic folder that the .csv file is located in

mouse_raw <- read.csv("WT vs KO_pep.csv") #Store the csv as a data.frame with the reference name
#Use the function read.csv() to read the file. Only need to use the "file_name.csv" if you have the correct working directory location. Otherwise you need the whole file path.

#View data set with environment dashboard and functions below
#head(mouse_raw) #See first six lines of the data.frame with this function
#summary(mouse_raw)
#str(mouse_raw) #Gives the class of each variable and the factor levels etc. 
################

#######################
#####Normalization#####
#######################
#B2RQQ1 = Myh6 & Q91Z83 = Myh7

#####Groups##### 
#KO_raw <- c(4:12)
WT_raw <- c(13:21) #Only WT needed for normalization

shared_myo <- "B2RQQ1; Q91Z83"   #the shared myosin peptides are used for normalization

mouse_shared_myo <- mouse_raw[mouse_raw$Master.Protein.Accessions == shared_myo, ]  #creates a d.f. with just shared myosin peptides

WT_raw_med <- NULL
for (i in 1:nrow(mouse_shared_myo)){
  WT_raw_med_temp <- median(as.numeric(mouse_shared_myo[i, WT_raw], na.rm = TRUE))
  WT_raw_med <- c(WT_raw_med, WT_raw_med_temp)
}
mouse_shared_myo$WT_raw_med <- WT_raw_med  #creates column with median values of raw WT for each peptide
myo_pep_norm = 15   #the number of peptides used for normalization
top_shared_myo <- mouse_shared_myo[order(mouse_shared_myo$WT_raw_med, decreasing = TRUE)[1:myo_pep_norm], ]

#sapply() to average normalization peptides
norm_value <- sapply(top_shared_myo[ ,4:ncol(top_shared_myo)], mean) #gets the average for each sample (column) for the normalization peptides selected


#Bind data.frame with raw and norm values 
mouse_norm <- mouse_raw[,4:ncol(mouse_raw)]/norm_value   #divide all columns by a specific number
colnames(mouse_norm) <- paste("norm", colnames(mouse_norm), sep = "_")   #create a prefix for all normalized colnames 
mouse <- cbind(mouse_raw, mouse_norm)

#######################
#######################
#######################


#####Normalized Abundance Groups#####
#Define the columns that correspond to the group
KO <- c(22:30)
WT <- c(31:39)
#####################################


######Median Columns#####
KO_med <- NULL # initialize the vector
WT_med <- NULL
for (i in 1:nrow(mouse)){
  KO_med_temp <- median(as.numeric(mouse[i, KO], na.rm = TRUE)) 
  WT_med_temp <- median(as.numeric(mouse[i, WT], na.rm = TRUE))
  KO_med <- c(KO_med, KO_med_temp)
  WT_med <- c(WT_med, WT_med_temp)
}
mouse$KO_med <- KO_med
mouse$WT_med <- WT_med
#########################


#####Standard Deviation#####
KO_sd <- NULL
WT_sd <- NULL
for (i in 1:nrow(mouse)){
  KO_sd_temp <- sd(as.numeric(mouse[i, KO], na.rm = TRUE)) 
  WT_sd_temp <- sd(as.numeric(mouse[i, WT], na.rm = TRUE))
  KO_sd <- c(KO_sd, KO_sd_temp)
  WT_sd <- c(WT_sd, WT_sd_temp)
}
mouse$KO_sd <- KO_sd
mouse$WT_sd <- WT_sd
############################


#####Ratio#####
KO.WT_ratio <- NULL
KO.WT_ratio <- KO_med / WT_med
mouse$KO.WT_ratio <- KO.WT_ratio
###############


#####Average across proteins#####
#?aggregate() or #?tapply()
#result_df <- aggregate(mouse$KO.WT_ratio ~ mouse$Master.Protein.Accessions, FUN = mean)
result_vector <- tapply(mouse$KO.WT_ratio, mouse$Master.Protein.Accessions, mean, na.rm = TRUE)

max_pep = 15   #the max number of peptides incuded in the average for the ratio
mouse_max15_pep <- mouse[order(mouse$WT_med, decreasing = TRUE)[1:max_pep], ]


#Ideas##
#Need to reduce the max number of peptides used for quantification to 15
#and restrict the protein groups to 5 or more peptides
pep_count <- table(mouse$Master.Protein.Accessions)   #gives how many peptides in each protein accension category
five <- pep_count > 5   #gives the proteins that have more than 5 peptides in each category
tab$Master.Protein.Accessions[five == TRUE] #not working...

####################
#####Statistics#####
####################

#create new columns with natural log of the 
ln_norm_abund <- mouse_norm[,1:ncol(mouse_norm)]   #divide all columns by a specific number
colnames(mouse_norm) <- paste("norm", colnames(mouse_norm), sep = "_")   #create a prefix for all normalized colnames 
mouse <- cbind(mouse_raw, mouse_norm)
norm_cols <- c()

