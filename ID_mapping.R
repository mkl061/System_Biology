### Package handling ----

# Bioconductor (i.e. BiocManager):
if (!require("BiocManager", quietly = T)) {
  install.packages("BiocManager", quiet = T)
}
BiocManager::install(update = T, ask = F, ... = c(quiet = T))


# # Biomartr:
# if (!require("biomartr", quietly = T)) {
#   BiocManager::install("Biostrings", update = T, ask = F, force = T)
#   BiocManager::install("biomaRt", update = T, ask = F)
#   
#   install.packages("biomartr", dependencies = T, quiet = T)
# }

if (!require("biomaRt", quietly = T)) {
  install.packages("biomaRt", quiet = T)
}

# Human data base:
BiocManager::install("org.Hs.eg.db", force = T)


# Tidyverse:
if (!require("tidyverse", quietly = T)) {
  install.packages("tidyverse", quiet = T)
}


### Working directory ----
if (Sys.info()[1] == "Linux") {
  setwd("/mnt/chromeos/MyFiles/Downloads")
} else if (Sys.info()[1] == "Windows") {
  setwd("C:/Users/Legion/Downloads")
}


### Custom function ----
# Splits apart a string at "|", removes duplicates, and splice it together again:
str_rm_duplicates <- function(string) {
  string <- string %>% 
    str_split(pattern = "\\|") %>% # Splitts the string at "|", into list
    unlist() %>% # Unlist to convert to vector
    base::unique() %>%  # Keep only the unique elements (i.e. removes duplicates)
    paste0(collapse = "|") # Collapse all elements (seperated by "|") into single string 
  return(string)
}

str_rm_NA <- function(string) {
  string <- string %>% 
    str_split(pattern = "\\|") %>%
    unlist() %>% 
    .[. != "NA"] %>% 
    paste0(collapse = "|")
  return(string)
}

# Function to reverse all elements in a character column:
str_rev_in_vec <- function(vec) {
  vec <- vec %>% 
    strsplit(NULL) %>% 
    lapply(rev) %>% 
    lapply(paste, collapse = "")
  
  return(vec)
}



### Input data ----
# Read files:
input_data <- read.csv("nodes.csv")

input_data$Common_name <- input_data$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

# Tidy up the file by removing unwanted white spaces:
for (col in 1:ncol(input_data)) {
  if (is.character(input_data[,col])) { # Checks if column is of character type
    input_data[,col] <- sapply(input_data[,col], 
                                   str_squish) # Removes white spaces in front, back and reduce multiple spaces to a single one
    #input_data[,col] <- sapply(input_data[,col], str_trim)
  }
}



### IDs with biomaRt ----
# listEnsembl()
# ensembl <- useEnsembl(biomart = "genes")
# datasets <- listDatasets(ensembl)

# Use only entries from Homo sapiens:
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attributes <- listAttributes(ensembl.con)
# filters <- listFilters(ensembl.con)

# Get IDs from other data bases, using the biomaRt package:
# OBS! This data frame will contain duplications, because
# some proteins may have multiple IDs within the same data bases.
multiple_ids <- getBM(attributes = c(
  "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
  "entrezgene_id", # I.e. NCBI gene ID
  "external_gene_name", 
  "external_synonym",
  "hgnc_id"), # HGCN ID
      filters = "uniprotswissprot",
      values = input_data$UniprotKB,
      mart = ensembl.con) #%>% 
  #add_column("dupl"=duplicated(.[,1])) %>% 
  #filter(dupl != T) #%>% View()


# Create a new df with removal of duplicates:
new_df <- input_data %>% dplyr::select(UniprotKB)
new_df$new <- NA
for (row in new_df$UniprotKB) {
  vec <- multiple_ids %>% 
    filter(uniprotswissprot == row) %>%
    dplyr::select(external_synonym) %>% 
    pull() %>% 
    paste0(collapse = "|")
  
  new_df$new[new_df$UniprotKB == row] <- vec
}
new_df <- new_df %>% 
  filter(UniprotKB != "") # Remove rows with empty cells

a <- left_join(input_data, new_df, by=c("UniprotKB" = "UniprotKB"))


### Find protein name and synonyms from Swiss-Prot db ----
# API URL for ALL human proteins in the Swiss-Prot database:
# https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29

# Checks if "Swiss_Prot.tsv" is in the working directory,
# if not it downloads it from the Uniprots webpage:
if ("Swiss_Prot.tsv" %in% dir()) {
  swiss_prot <- read.delim("Swiss_Prot.tsv")  
} else {
  cat("OBS! The Swiss_Prot file is not in working directory. Fetching it from web, but it takes a while")
  swiss_prot <- read.delim("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29")
}


# All the Uniprot IDs from the spreadsheet:
Uni_ids <- input_data %>% 
  dplyr::select(UniprotKB) %>% 
  pull() # Results in character vector

# Subset of swiss_prot df:
sub_df <- swiss_prot %>% .[.$Entry %in% Uni_ids, ]

# Create new columns:
sub_df$protein_name <- NA
sub_df$synonyms <- NA

# Filling the new columns with protein name and synonyms:
for (row in 1:nrow(sub_df)) { # Iterate through all the rows
  # Find the cell containing the preferred protein name, and synonyms:
  string <- sub_df[row,] %>% 
    dplyr::select(Protein.names) %>% # Focus on only the Protein.names column
    pull() # Results in a string
  
  # Take the string and cuts out the preferred protein name:
  p_name <- string %>% # Take the string
    str_sub(
      1, # First character 
      str_locate(., pattern = " \\(")[1] # Ends at first " ("
    ) 
  
  # Take the string not containing the preferred protein name,
  # cut it up into the synonyms, which are covered with "( )",
  # remove the parentheses, and return a characther vector:
  syn <- string %>% 
    str_sub(
      str_locate(., pattern = "\\(")[1]+1, # starts after the first "("
      nchar(string) # Ends at end of string
    ) %>% 
    str_replace_all(., pattern = "\\)", replacement = "") %>% # Replace ")"
    str_split(., pattern = " \\(") %>% unlist()
  
  # Remove possible EC-number from the synonyms: 
  if (!is.na(syn[1])) {
    flag <- numeric(0)
    for (i in 1:length(syn)) {
      if (str_starts(syn[i], "EC")) {
        flag <- base::append(flag, i)
      }
    }
    
    
    if (length(flag > 0)) {
      syn <- syn[-c(flag)]
    }
  }
  
  # Collapse the synonyms with the "|" between:
  syn <- paste0(syn, collapse = "|")
  
  # Take the entry name (e.g. TNFA_HUMAN) and remove "_HUMAN" 
  ent_name <- sub_df[row, "Entry.Name"] %>% 
    str_replace(pattern = "_HUMAN", replacement = "")
  
  # Add the entry name to the other synonyms:
  syn <- str_c(ent_name, syn, sep = "|")
  
  # Add the protein name and the synonyms to their columns:
  sub_df[row, "protein_name"] <- p_name
  sub_df[row, "synonyms"] <- syn
}



### Result from BioGateway query ----
## Reads a BioGateway network's edge table. 
## Returns a df of only source, target and interaction:
read_BioGateway <- function(df) {
  org <- df$Edge.Id # http://rdf.biogateway.eu/gene/9606/IL1R1;gene::http://semanticscience.org/resource/SIO_010078;http://rdf.biogateway.eu/prot/9606/B8ZZ73
  int <- df$Source.Graph # gene, prot2bp, or tfac2gene
  
  ## Source:
  # Index of first ";":
  cut1 <- org %>% 
    str_locate(.,
               str_c(";",
                     int)) %>% 
    .[,1]-1
  
  # Create a new string which drops all unwanted information
  # until cut1 (i.e. the index number):
  new <- org %>% 
    str_sub(end = cut1)
  
  # Reveres the new string, locate the first "/" which marks
  # the next wanted information, stores the index value:
  cut2 <- new %>%
    str_rev_in_vec() %>% 
    str_locate("/") %>% 
    .[,1]-1
  
  # Reverse the new string,cut all the information after cut2,
  # and reverse the string back to normal. Results in the
  # source node:
  source <- new %>% 
    str_rev_in_vec() %>% 
    str_sub(end = cut2) %>% 
    str_rev_in_vec()
  
  
  ## Target:
  # Reverse the original string, locate first "/", and store 
  # the index value:
  cut3 <- org %>% 
    str_rev_in_vec() %>% 
    str_locate("/") %>% 
    .[,1]-1
  
  # Reverse the original string, remove the information to cut3,
  # then reverse it back to normal:
  target <- org %>%
    str_rev_in_vec() %>% 
    str_sub(end = cut3) %>% 
    str_rev_in_vec()
  
  # Create a data frame with the source, tatget, and interaction type:
  df <- cbind.data.frame(
    "source"=unlist(source),
    "target"=unlist(target),
    "int"=unlist(int)
  )
  
  return(df)
}


BG_TFs_NFKB1_and_RELA <- read.csv("Marius_edges_test.csv")

res <- read_BioGateway(BG_TFs_NFKB1_and_RELA) %>% 
  filter(int != "prot2bp") # Don't need the GO-term information

genes <- res %>% 
  filter(int == "gene") %>% 
  filter(target %in% swiss_prot$Entry)


a <- left_join(genes, input_data[, c("UniprotKB", "id")], 
               by=c("target" = "UniprotKB"))


### Other ----

## Your spreadsheet:
#dat <- read.csv("original_set.csv")
dat <- read.csv("nodes.csv")
dat$Common_name <- dat$Common_name %>% str_replace("\xa0", "") # Fixing a recurring problem

## Trim all whitespaces in all cells with strings:
for (col in 1:ncol(dat)) {
  if (is.character(dat[,col])) {
    dat[,col] <- sapply(dat[,col], str_trim)
  }
}





ID_mapping <- function(from, to, df, keys_column, from_column) {
  ## Uniprot IDs given when mapping by "gene symbol":
  ids_from_db <- mapIds(org.Hs.eg.db,
                      #keys = df$id,
                      keys = unlist(df[keys_column]),
                      column = to,
                      keytype = from,
                      multiVals = "list") # A gene symbol may refer to multiple IDs (return a vector of these)
  
  
  # We usually start by finding Uniprot IDs!
  if (to == "UNIPROT") {
    # Checks if the Swiss_Prot.tsv file is in the directory.
    # If not, gets it from the web via API URL:
    if ("Swiss_Prot.tsv" %in% dir()) {
      swiss_prot <- read.delim("Swiss_Prot.tsv")  
    } else {
      cat("OBS! The Swiss_Prot file is not in working directory. Fetching it from web, but it takes a while")
      swiss_prot <- read.delim("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29")
    }
    
    
    for (lst_ind in 1:length(ids_from_db)) { # List index val
      new_vec <- vector() # Empty vector
      for (vec_elm in ids_from_db[[lst_ind]]) { # Elements of the list-element
        if (vec_elm %in% swiss_prot$Entry == T) { # Checks if element (i.e. Uniprot ID) is in Swiss-Prot (i.e. is it reviewed by expert)
          new_vec <- append(new_vec, vec_elm) # Adds the element to the vector
        }
      }
      
      if (length(new_vec) == 0) { # If the vector is empty, set NA
        new_vec <- NA
      }
      
      ids_from_db[[lst_ind]] <- new_vec
      
    }
  }
  
  ##############################################################################
  # OBS! It is missing some lines of code HERE!
  # Because the next line of code only works if all the list elements only
  # containt 1 element each! Which it seems to do for Swiss-Prot
  ##############################################################################
  
  
  ## Unlist the list:
  ids_from_db <- unlist(ids_from_db)
  
  
  tabl_check <- df %>% 
    dplyr::select("id"=!!sym(keys_column), "Our_ID"=!!sym(from_column)) %>% 
    add_column("Mapping_ID" = ids_from_db) %>% 
    mutate(final = ifelse(Our_ID == "", # Is Our_ID an empty string?
                          Mapping_ID, # Yes; Set as the ID gotten from ID mapping
                          ifelse(is.na(Our_ID),  # No; Initiate another ifelse(). Is Our_ID a NA?
                                 Mapping_ID,  # Yes; Set as the ID gotten from ID mapping
                                 Our_ID # No; keep Our_ID as it was
                                 ) 
                          )
          )
  
  df[from_column] <- tabl_check$final
  
  return(df)
}


View(ID_mapping("SYMBOL", "UNIPROT", dat, keys_column = "id", from_column = "UniprotKB"))
ny <- ID_mapping("SYMBOL", "UNIPROT", dat, keys_column = "id", from_column = "UniprotKB")
liste <- ID_mapping("UNIPROT", "ENSEMBL", ny, keys_column = "UniprotKB", from_column = "Ensembl")
View(ID_mapping("UNIPROT", "SYMBOL", dat, keys_column = "UniprotKB", from_column = "Ensembl"))



## Uniprot IDs given when mapping by "gene symbol":
ids_from_db <- mapIds(org.Hs.eg.db,
              keys = dat$id,
              column = "UNIPROT",
              keytype = "SYMBOL",
              multiVals = "list") # A gene symbol may refer to multiple IDs (return a vector of these)


## Given that ids_from_db is a list, we need to find a single identifier for
## each protein before we can turn the list into a vector:
uni_sec <- ids_from_db


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
