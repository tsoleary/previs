library(seqinr)

fasta <- read.fasta("Homo_sapiens_Uniprot (taxonomy 9606).fasta", 
                    seqtype = "AA", as.string = TRUE)

pro_des <- str_extract(attr(fasta[[1]], "Annot"), "HUMAN\\s[[:print:]]+\\sOS=") %>%
  gsub("\\sOS=", "", .) %>% gsub("MOUSE\\s", "", .)
                 
Accession <- NULL
Gene <- NULL
Description <- NULL
                 
for (i in 1:length(fasta)){
  pro_acc <- str_extract(attr(fasta[[i]], "name"), "\\|[[:alnum:]]+\\|") %>%
               gsub("\\|", "", .)
  pro_gene <- str_extract(attr(fasta[[i]], "Annot"), "GN=[[:alnum:]]+") %>%
                gsub("GN=", "", .)
  pro_des <- str_extract(attr(fasta[[i]], "Annot"), "HUMAN\\s[[:print:]]+\\sOS=") %>%
    gsub("\\sOS=", "", .) %>% gsub("HUMAN\\s", "", .)  
  
  if (is.na(pro_gene) == TRUE){
    pro_gene <- pro_des
  }
  
  Accession <- c(Accession, pro_acc)
  Gene <- c(Gene, pro_gene)
  Description <- c(Description, pro_des)
  
  gene_df <- cbind(Accession, Gene, Description)
  
}                

write.csv(gene_df, "human_fasta_accession_gene.csv")                 
                 
                 
                 
