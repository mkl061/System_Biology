library(tidyverse)

setwd("C:/Users/Legion/Downloads")

if ("Swiss_Prot.tsv" %in% dir()) {
  swiss_prot <- read.delim("Swiss_Prot.tsv")  
} else {
  cat("OBS! The file is not in working directory. Fetching it from web, but it takes a while")
  swiss_prot <- read.delim("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29")
}


data1 <- read.csv("nodes.csv")
data2 <- read.csv("edges.csv")
data3 <- read.csv("Marius_edges_test.csv")

# str_rev <- function(string) {
#   splt_str <- strsplit(string, NULL)[[1]]
#   rev_str <- paste(rev(splt_str), collapse = "")
#   
#   return(rev_str)
# }











# Test
a <- read_BioGateway(data3) 

b <- a %>% 
  filter(
    target %in% swiss_prot$Entry
    )
  


