library(tidyverse)
setwd("C:/Users/PrevBeast/Documents/R/Esser/Exon Organization and Annotation")

features <- read_file("Human Titin Gene Features.txt")

str_locate_all(features, "exon")
test <- str_extract_all(features, "\\b[:alnum:]+\\b") %>%
  unlist(.)

# To-do: Exon Lengths

for (i in 1:length(test)) {
  notegrab <- FALSE
  notetemp <- NA
  if (i == 1){

    weirdcount <- 0
    flipouts2 <- c("start")
    params <- data.frame(start = numeric(), end = numeric(), number = numeric(),
                         note = character())
  }
  if (test[i] == "exon"){
    if (is.na(as.numeric(test[i+1]))|is.na(as.numeric(test[i+2]))) {
      # print(i)
      # print(test[i-1])
      # flipouts2 <- c(flipouts2, test[i-1])
      # weirdcount <- (weirdcount + 1)
      next
    } else {
      print(paste("here at", i))
      
      for (j in i:length(test)){
        if (test[j] == "note") {
          notetemp <- character()
          notegrab <- TRUE
          next
        }
        if (notegrab == TRUE & test[j] != "number") {
          notetemp <- paste(notetemp, test[j])
        }
        if (test[j] == "number") {
          print(j)
          notegrab <- FALSE
          print(notetemp)
          numtemp <- test[j+1]
          break
        }
      }
      paramtemp <- data.frame(start = as.numeric(test[i+1]), end = as.numeric(test[i+2]),
                              number = as.numeric(numtemp), note = as.character(notetemp))
      params <- full_join(params, paramtemp)
    }
    
    # print(test[i+1])
    # print(test[i+2])
  }
  # if (i == 4){
  #   break
  # }
  
}

write.csv(params, file = "Human exon positions test.csv")

# exons <- read.csv("Ensemble Exon numbers.csv")
# seq <- read_file("Titin Lance.txt")

# splice <- str_sub(seq, exons[1,1], exons[1,2])
# splice <- paste0(splice, str_sub(seq, exons[2,1], exons[2,2]))


for (i in 1:nrow(exons)) {
  if (i == 1){
    splice <- str_sub(seq, exons$Start[i], exons$End[i])
  }
  if (i > 1){
    splice <- paste0(splice, str_sub(seq, exons$Start[i], exons$End[i]))
  }
}

str_sub(splice, str_length(splice)-8, str_length(splice))

write_file(splice, "splice_seq2.txt")

ex2 <- str_sub(seq, exons$Start[3], exons$End[3])


cool <- "cool and beans 81"
str_extract_all(cool, "\\b[:alnum:]+\\b")
str_sub(seq, 76158, 82385)



shopping_list <- c("apples x4", "bag of flour", "bag of sugar", "milk x2")
str_extract(shopping_list, "\\d")
str_extract(shopping_list, "[a-z]+")
str_extract(shopping_list, "[a-z]{1,4}")
str_extract(shopping_list, "\\b[a-z]{1,4}\\b")

# Extract all matches
str_extract_all(shopping_list, "[a-z]+")
str_extract_all(shopping_list, "\\b[:alnum:]+\\b")
str_extract_all(shopping_list, "\\d")

# Simplify results into character matrix
str_extract_all(shopping_list, "\\b[a-z]+\\b", simplify = TRUE)
str_extract_all(shopping_list, "\\d", simplify = TRUE)

# Extract all words
str_extract_all("This is, suprisingly, a sentence.", boundary("word"))


regexpr("a", shopping_list[1])
