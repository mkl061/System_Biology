library(tidyverse)


#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)


setwd("C:/Users/Legion/Downloads")
# API URL for ALL human proteins in the Swiss-Prot database:
# https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29

swiss_prot <- read.delim("Swiss.tsv")

## Your spreadsheet:
#dat <- read.csv("original_set.csv")
dat <- read.csv("nodes.csv")
dat$Common_name <- dat$Common_name %>% str_replace("\xa0", "") # Fixing a recurring problem

## Trim all whitespaces in all cells with strings:
for (col_ind in 1:ncol(dat)) {
  if (is.character(dat[,col_ind])) {
    dat[,col_ind] <- sapply(dat[,col_ind], str_trim)
  }
}


ID_mapping <- function(from, to, df, keys_column, from_column) {
  ## Uniprot IDs given when mapping by "gene symbol":
  uni_first <- mapIds(org.Hs.eg.db,
                      #keys = df$id,
                      keys = unlist(df[keys_column]),
                      column = to,
                      keytype = from,
                      multiVals = "list") # A gene symbol may refer to multiple IDs (return a vector of these)
  
  
  # We usually start by finding Uniprot IDs!
  if (to == "UNIPROT") {
    swiss_prot <- read.delim("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29")
    
    
    for (lst_ind in 1:length(uni_first)) { # List index val
      new_vec <- vector() # Empty vector
      for (vec_elm in uni_first[[lst_ind]]) { # Elements of the list-element
        if (vec_elm %in% swiss_prot$Entry == T) { # Checks if element (i.e. Uniprot ID) is in Swiss-Prot (i.e. is it reviewed by expert)
          new_vec <- append(new_vec, vec_elm) # Adds the element to the vector
        }
      }
      
      if (length(new_vec) == 0) { # If the vector is empty, set NA
        new_vec <- NA
      }
      
      uni_first[[lst_ind]] <- new_vec
      
    }
  }
  
  
  ## Unlist the list:
  uni_first <- unlist(uni_first)
  
  
  tabl_check <- df %>% 
    dplyr::select("id"=!!sym(keys_column), "Our_ID"=!!sym(from_column)) %>% 
    add_column("Mapping_ID" = uni_first) %>% 
    mutate(final = ifelse(Our_ID == "", # Is Our_ID an empty string?
                          Mapping_ID, # Yes; Set as the ID gotten from ID mapping
                          ifelse(is.na(Our_ID),  # No; Initiate another ifelse(). Is Our_ID a NA?
                                 Mapping_ID,  # Yes; Set as the ID gotten from ID mapping
                                 Our_ID # No; keep Our_ID as it was
                                 ) 
                          )
          )
  
  df$UniprotKB <- tabl_check$final
  
  return(df)
}


View(ID_mapping("SYMBOL", "UNIPROT", dat, keys_column = "id", from_column = "UniprotKB"))
ny <- ID_mapping("SYMBOL", "UNIPROT", dat, keys_column = "id", from_column = "UniprotKB")
liste <- ID_mapping("UNIPROT", "ENSEMBL", ny, keys_column = "UniprotKB", from_column = "Ensembl")
View(ID_mapping("UNIPROT", "SYMBOL", dat, keys_column = "UniprotKB", from_column = "Ensembl"))



## Uniprot IDs given when mapping by "gene symbol":
uni_first <- mapIds(org.Hs.eg.db,
              keys = dat$id,
              column = "UNIPROT",
              keytype = "SYMBOL",
              multiVals = "list") # A gene symbol may refer to multiple IDs (return a vector of these)


## Given that uni_first is a list, we need to find a single identifier for
## each protein before we can turn the list into a vector:
uni_sec <- uni_first


for (lst_ind in 1:length(uni_sec)) {
  new_vec <- vector()
  for (vec_elm in uni_sec[[lst_ind]]) {
    if (vec_elm %in% swiss_prot$Entry == T) {
      new_vec <- append(new_vec, vec_elm)
    }
  }
  
  if (length(new_vec) == 0) {
    new_vec <- NA
  }
  
  uni_sec[[lst_ind]] <- new_vec
  
  df <- df %>% 
    mutate(UniprotKB = ifelse(
      tabl_check$Our_ID == "" && (tabl_check$Mapping_ID != "" || is.na(tabl_check$Mapping_ID)),
      tabl_check$final <- tabl_check$Mapping_ID,
      tabl_check$final <- tabl_check$Our_ID
    ))
  
  return(df)
}

## Unlist the list:
uni_sec <- unlist(uni_sec)

# if (all(lengths(uni_sec))) {
#   uni_sec <- unlist(uni_sec)
# } else {
#   cat("OBS! Some proteins have multiple identifiers in Swiss-Prot")
#   view(
#     uni_sec[which(lengths(uni_sec) != 1)]  
#   )
# }

#dat.copy <- dat
#dat$UniprotKB[1:5] <- ""


tabl_check <- cbind.data.frame(dat$id,
                               dat$UniprotKB,
                               uni_sec)
base::colnames(tabl_check) <- c("id", "Our_ID", "Mapping_ID")
tabl_check$identical <- tabl_check$Our_ID == tabl_check$Mapping_ID

# OBS!!! These are the non-identical ones:
tabl_check %>% filter(identical == F)

tabl_check$final <- NA
for (i in 1:nrow(tabl_check)) {
  if (tabl_check$Our_ID[i] == "" && (tabl_check$Mapping_ID[i] != "" || is.na(tabl_check$Mapping_ID[i]))) {
    tabl_check$final[i] <- tabl_check$Mapping_ID[i]
  } else {
    tabl_check$final[i] <- tabl_check$Our_ID[i]
  }
}

merg_table <- cbind.data.frame(tabl_check$id, tabl_check$final)
colnames(merg_table) <- c("id", "UniprotKB")

new_dat <- dat
new_dat$UniprotKB <- tabl_check$final

hei <- hei %>% 
  mutate(UniprotKB.x = ifelse((UniprotKB.y=="" || is.na(UniprotKB)),))
#hei <- merge(dat, merg_table, by="id", all.x = T)
#dat2 <- left_join(dat, merg_table)
