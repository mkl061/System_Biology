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

# Function to reverse all elements in a character column:
str_rev_in_vec <- function(vec) {
  vec <- vec %>% 
    strsplit(NULL) %>% 
    lapply(rev) %>% 
    lapply(paste, collapse = "")
  
  return(vec)
}



## Reads a BioGateway network's edge table. 
## Returns a df of only source, target and interaction:
read_BioGateway <- function(df) {
  org <- df$Edge.Id
  int <- df$Source.Graph
  
  # Source:
  cut1 <- org %>% 
    str_locate(.,
               str_c(";",
                     int)) %>% 
    .[,1]-1
  
  new <- org %>% str_sub(end = cut1)
  
  cut2 <- new %>%
    str_rev_in_vec() %>% 
    str_locate("/") %>% 
    .[,1]-1
  
  source <- new %>% 
    str_rev_in_vec() %>% 
    str_sub(end = cut2) %>% 
    str_rev_in_vec()
  
  
  # Target:
  
  cut3 <- org %>% 
    str_rev_in_vec() %>% 
    str_locate("/") %>% 
    .[,1]-1
  
  target <- org %>%
    str_rev_in_vec() %>% 
    str_sub(end = cut3) %>% 
    str_rev_in_vec()
  
  
  df <- cbind.data.frame(
    "source"=unlist(source),
    "target"=unlist(target),
    "int"=unlist(int)
  )
  
  return(df)
}






# Test
a <- read_BioGateway(data3) 

b <- a %>% 
  filter(
    target %in% swiss_prot$Entry
    )
  

b <- a %>% 
  filter(int == "gene") %>% 
  filter(source == "IRAK1") %>%
  filter(target %in% swiss_prot$Entry)


