library(tidyverse)

setwd("C:/Users/Legion/Downloads")

data1 <- read.csv("nodes.csv")
data2 <- read.csv("edges.csv")
data3 <- read.csv("Marius_edges_test.csv")

str_rev <- function(string) {
  splt_str <- strsplit(string, NULL)[[1]]
  rev_str <- paste(rev(splt_str), collapse = "")
  
  return(rev_str)
}

str_rev_in_vec <- function(vec) {
  vec <- vec %>% 
    strsplit(NULL) %>% 
    lapply(rev) %>% 
    lapply(paste, collapse = "")
  
  return(vec)
}



source_target_int <- function(row, df=data3) {
  org <- df$Edge.Id[row]
  int <- df$Source.Graph[row]
  
  # Source:
  cut1 <- org %>% 
    str_locate(., 
               str_c(";", 
                     int)) %>% 
    .[1]-1
  
  new <- org %>% str_sub(end = cut1)
  
  cut2 <- new %>%
    str_rev() %>% 
    str_locate("/") %>% 
    .[1]-1
  
  source <- new %>% 
    str_rev() %>% 
    str_sub(end = cut2) %>% 
    str_rev()
  
  
  # Target:
  
  cut3 <- org %>% 
    str_rev() %>% 
    str_locate("/") %>% 
    .[1]-1
  
  df$target <- org %>%
    str_rev() %>% 
    str_sub(end = cut3) %>% 
    str_rev()
  
  return(c(source, df$target, int))
  
}

a <- source_target_int(1)
b <- data1 %>% 
  select(id, UniprotKB) %>% 
  filter(id == a[1]) %>% 
  select(UniprotKB) %>%
  .[1]

if (a[2] == b) {
  cat("TRUE")
} else {
  cat("FALSE")
}


df_source_target_int <- function(df=data3) {
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
  
  
  return(cbind.data.frame(unlist(source), unlist(target), unlist(int)))
}


a <- df_source_target_int() %>% filter(source == "IL1R1")




cut1 <- data3$Edge.Id %>% 
  str_locate(.,
             str_c(";",
                   data3$Source.Graph)) %>% 
  .[,1]-1
cut1

new <- data3$Edge.Id %>% str_sub(end = cut1)
new

cut2 <- new %>%
  str_rev_in_vec() %>% 
  str_locate("/") %>% 
  .[,1]-1
cut2


source <- new %>% 
  str_rev_in_vec() %>% 
  str_sub(end = cut2) %>% 
  str_rev_in_vec()
source


cut3 <- data3$Edge.Id %>% 
  str_rev_in_vec() %>% 
  str_locate("/") %>% 
  .[,1]-1
cut3

target <- data3$Edge.Id %>%
  str_rev_in_vec() %>% 
  str_sub(end = cut3) %>% 
  str_rev_in_vec()
target

a <- unlist(source)
