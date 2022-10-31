### Package handling ----

# Bioconductor (i.e. BiocManager):
if (!require("BiocManager", quietly = T)) {
  install.packages("BiocManager", quiet = T)
}
#BiocManager::install(update = T, ask = F, ... = c(quiet = T))
# OBS! The command above usually give a warning, but apparantly this
# can be ignored. See: https://support.bioconductor.org/p/128825/ 



if (!require("biomaRt", quietly = T)) {
  #install.packages("biomaRt", quiet = T)
  BiocManager::install("biomaRt")
}



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
# # Splits apart a string at "|", removes duplicates, and splice it together again:
# str_rm_duplicates <- function(string) {
#   string <- string %>% 
#     str_split(pattern = "\\|") %>% # Splits the string at "|", into list
#     unlist() %>% # Unlist to convert to vector
#     base::unique() %>%  # Keep only the unique elements (i.e. removes duplicates)
#     paste0(collapse = "|") # Collapse all elements (seperated by "|") into single string 
#   return(string)
# }
# 
# str_rm_NA <- function(string) {
#   string <- string %>% 
#     str_split(pattern = "\\|") %>%
#     unlist() %>% 
#     .[. != "NA"] %>% 
#     paste0(collapse = "|")
#   return(string)
# }



# Fucntion that looks throught all characther columns, and removes
# white spaces at front and end, and convert multiple spaces to single spaces:
remove_whitespaces <- function(df) {
  for (col in 1:ncol(df)) {
    if (is.character(df[,col])) { # Checks if column is of character type
      df[,col] <- sapply(
        df[,col],
        str_squish) # Removes white spaces in front, back and reduce multiple spaces to a single one
    }
  }
  return(df)
}



# Function to reverse all elements in a character column:
str_rev_in_vec <- function(vec) {
  vec <- vec %>% 
    strsplit(NULL) %>% 
    lapply(rev) %>% 
    lapply(paste, collapse = "")
  
  return(vec)
}




# Function to collapse strings of two different columns (same row) with
# specified separator:
col_splice_to_string <- function(df, row, col1, col2, from_sep, to_sep) {
  first <- df[row, col1] 
  second <- df[row, col2] 
  
  first_NA <- is.na(first)
  second_NA <- is.na(second)
  
  if (first_NA == F && second_NA == F) { # None are NA
    # None of the variables are equal to NA.
    # However, they may be equal to "" (i.e. they may be empty strings)
    
    # Create variables, holding either TRUE or FALSE, depending on
    # if the input variables are equal to "" or not:
    first_empty <- first == ""
    second_empty <- second == ""
    
    if (first_empty == F && second_empty == F) { # None are empty (i.e. "")
      # None of them are empty, then we can proceed to combining them:
      final <- str_flatten(c(first, second), collapse = from_sep)
      # Check if "|" is used as separator. If so, given that "|" is
      # a special character, we need to putt "\\" in front:
      if (from_sep == "|") {
        final <- final %>% str_split(pattern = "\\|") # splits (at "|") into list of vectors
      } else { # The separator is something else:
        final <- final %>% str_split(pattern = from_sep) # splits (at from_sep) into list of vectors
      }
      final <- final %>%
        base::unlist() %>% # Unlist to get the single vector
        base::unique() %>% # Keep only the unique elements of the vector
        str_flatten(collapse = to_sep) # Flattens the vector to a single string, with each element separated by the to_sep
      
      return(final) # Return the string
    } else if (first_empty == T && second_empty == T) {
      return(NA)
    } else if (first_empty == T) { # First is empty, but not Second
      return(second)
    } else { # Second is empty, but first is not
      return(first)
    }
    
    
  } else if (first_NA == T && second_NA == T) { # Both are NA
    return(NA)
  } else if (first_NA == T) { # First is NA, Second is not
    if (second == "") { # Checks if Second could be empty
      return(NA) # If so, return NA
    } else { 
      # Second is neither NA nor "", but we need to check if it 
      # consists of multiple terms or not:
      if (str_detect(str_squish(second), from_sep) == F) { # Is from_sep in Second
        return(second)
      } else { # from_sep is in Second
        # Check if "|" is used as separator. If so, given that "|" is
        # a special character, we need to putt "\\" in front:
        if (from_sep == "|") {
          second <- second %>% str_split(pattern = "\\|") # splits (at "|") into list of vectors
        } else { # The separator is something else:
          second <- second %>% str_split(pattern = from_sep) # splits (at from_sep) into list of vectors
        }
        second <- second %>%
          base::unlist() %>% # Unlist to get the single vector
          base::unique() %>% # Keep only the unique elements of the vector
          str_flatten(collapse = to_sep) # Flattens the vector to a single string, with each element separated by the to_sep
        
        return(second) # Return the string
      }
    }
  } else { # Second is NA, First is not
    if (first == "") { # Checks if first could be empty
      return(NA) # If so, return NA
    } else { 
      # first is neither NA nor "", but we need to check if it 
      # consists of multiple terms or not:
      if (str_detect(str_squish(first), from_sep) == F) { # Is from_sep in first
        return(first)
      } else { # from_sep is in first
        # Check if "|" is used as separator. If so, given that "|" is
        # a special character, we need to putt "\\" in front:
        if (from_sep == "|") {
          first <- first %>% str_split(pattern = "\\|") # splits (at "|") into list of vectors
        } else { # The separator is something else:
          first <- first %>% str_split(pattern = from_sep) # splits (at from_sep) into list of vectors
        }
        first <- first %>%
          base::unlist() %>% # Unlist to get the single vector
          base::unique() %>% # Keep only the unique elements of the vector
          str_flatten(collapse = to_sep) # Flattens the vector to a single string, with each element separated by the to_sep
        
        return(first) # Return the string
      }
    }
  }
}




# Function that will be used to fill cells of one column if the cell contain
# either NA or "": 
fill_cells <- function(to_col, from_col) {
  ret_col <- ifelse(
    is.na(to_col) == T,
    from_col,
    ifelse(
      to_col == "",
      from_col,
      to_col
    )
  )
  
  return(ret_col)
}



### Input data ----
## Nodes:

# Read files:
input_data_nodes <- read.csv("nodes.csv")

input_data_nodes$Common_name <- input_data_nodes$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

# Tidy up the file by removing unwanted white spaces:
input_data_nodes <- remove_whitespaces(input_data_nodes)



## Edges:
input_data_edges <- read.csv("edges.csv")

input_data_edges <- remove_whitespaces(input_data_edges)

# Fix capital letter in name:
input_data_edges <- input_data_edges %>% 
  mutate(Curator = str_to_title(Curator))


### IDs with biomaRt ----
# listEnsembl()
# ensembl <- useEnsembl(biomart = "genes")
# datasets <- listDatasets(ensembl)

# Use only entries from Homo sapiens:
# NB! This step takes a while, so I put it in an if-statement so that
# it doesn't have to run each time if its already in the global environment.
if ("ensembl.con" %in% ls() == F) {
  ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}

# attributes <- listAttributes(ensembl.con)
# filters <- listFilters(ensembl.con)

# Get IDs from other data bases, using the biomaRt package:
# OBS! This data frame will contain duplications, because
# some proteins may have multiple IDs within the same data bases.
multiple_ids <- getBM(attributes = c(
  "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
  "entrezgene_id", # I.e. NCBI gene ID
  "external_gene_name", # I.e. Gene symbol 
  "ensembl_gene_id", # Ensembl IDs
  "hgnc_id"), # I.e. HGCN ID
      filters = "uniprotswissprot",
      values = input_data_nodes$UniprotKB,
      mart = ensembl.con) %>% 
  mutate(dup = duplicated(uniprotswissprot)) %>% # OBS! Because Ensemble have multiple IDs per protein, we choose to only keep one
  filter(dup != T) %>% 
  dplyr::select(-dup)

# Creates a data frame with all the synonyms:
# Each synonym has its own row, meaning that there are multiple rows
# with the same id/value in the "uniprotswissprot" column.
raw_gene_syn_df <- getBM(attributes = c(
  "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
  "external_synonym"), # I.e. gene synonyms
  filters = "uniprotswissprot",
  values = input_data_nodes$UniprotKB,
  mart = ensembl.con)

# Takes the raw_gene_syn_df, and combine all the synonym for every Uniprot ID
# to a single row:
# "fin_gene_syn_df stands for "finished gene synonym data frame"
fin_gene_syn_df <- input_data_nodes %>% 
  dplyr::select(UniprotKB) # Uniprot IDs from input data

fin_gene_syn_df$gene_synonyms_BioMart <- NA # Create a new column

for (U_ID in fin_gene_syn_df$UniprotKB) { # Iterate trough all Uniprot IDs
  vec <- raw_gene_syn_df %>% # Get raw data
    filter(uniprotswissprot == U_ID) %>% # Filter to only rows with current Uniprot ID
    dplyr::select(external_synonym) %>% # Select only synonym column
    pull() %>% # Converts data frame column to character vector 
    paste0(collapse = "|") # Collapse all elements of vector to single string
  
  # Assigns the string as the cell-content of the column (i.e. gene synonym)
  # and the row (i.e. the one with column UniprotKB == the Uniprot ID):
  fin_gene_syn_df$gene_synonyms_BioMart[fin_gene_syn_df$UniprotKB == U_ID] <- vec
}

# Remove rows with empty Uniprot ID:
fin_gene_syn_df <- fin_gene_syn_df %>% 
  filter(UniprotKB != "") 

# Combine the data frame with different IDs with the the data frame
# with all the gene synonyms:
BioMart_df <- left_join(
  multiple_ids, fin_gene_syn_df, 
  by=c("uniprotswissprot" = "UniprotKB"))


## Tidying up:
rm(multiple_ids, 
   raw_gene_syn_df, 
   fin_gene_syn_df,
   U_ID,
   vec
   )



### Find protein name and synonyms from Swiss-Prot db ----
# API URL for ALL human proteins in the Swiss-Prot database:
# https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29

# Checks if "Swiss_Prot.tsv" is in the working directory,
# if not it downloads it from the Uniprot's webpage:
if ("Swiss_Prot.tsv" %in% dir()) {
  swiss_prot <- read.delim("Swiss_Prot.tsv")  
} else {
  cat("OBS! The Swiss_Prot file is not in working directory. Fetching it from web, but it takes a while")
  swiss_prot <- read.delim("https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cgene_synonym&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29")
}


# All the Uniprot IDs from the spreadsheet:
Uniprot_IDs <- input_data_nodes %>%
  filter(UniprotKB != "") %>% 
  dplyr::select(UniprotKB) %>% 
  pull() # Results in character vector

# Subset of swiss_prot df, with only IDs equal to those in our input data:
sub_df <- swiss_prot %>% .[.$Entry %in% Uniprot_IDs, ]

# Create new columns:
sub_df$protein_name <- NA
sub_df$protein_synonyms_swiss <- NA
sub_df$gene_synonyms_swiss <- NA

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
  # remove the parentheses, and return a character vector:
  syn <- string %>% 
    str_sub(
      str_locate(., pattern = "\\(")[1]+1, # starts after the first "("
      nchar(string) # Ends at end of string
    ) %>% 
    str_replace_all(., pattern = "\\)", replacement = "") %>% # Replace ")"
    str_split(., pattern = " \\(") %>% unlist()
  
  # Remove possible EC-number from the synonyms: 
  if (!is.na(syn[1])) {
    flag <- numeric(0) # Vector to hold index values
    for (i in 1:length(syn)) { # Iterate over all elements of "syn" vector
      if (str_starts(syn[i], "EC")) { # If element starts with "EC"
        flag <- base::append(flag, i) # Append index value for the element, to flag 
      }
    }
    
    # Checks if any elements of "syn" started with "EC" (i.e. is length of flag > 0)
    if (length(flag > 0)) {
      syn <- syn[-c(flag)] # Removes the element with index values in flag
    }
  }
  
  # Remove possible NAs that have been included:
  syn <- syn[!is.na(syn)]
  
  # Collapse the synonyms with the "|" between:
  syn <- paste0(syn, collapse = "|")
  
  # Take the entry name (e.g. TNFA_HUMAN) and remove "_HUMAN" 
  ent_name <- sub_df[row, "Entry.Name"] %>% 
    str_replace(pattern = "_HUMAN", replacement = "")
  
  # Add the entry name to the other synonyms:
  # There may not be any synonyms in the data base, and then we only
  # ad the entry name:
  if (syn != "") {
    syn <- str_c(ent_name, syn, sep = "|")
  } else {
    syn <- ent_name
  }
  
  # Use the custom function to collapse the gene name and all its synonyms
  # together in a single string:
  g_syn <- col_splice_to_string(df = sub_df,
              row = row,
              col1 = "Gene.Names", 
              col2 = "Gene.Names..synonym.",
              from_sep = " ", 
              to_sep = "|")
  
  
  # Assign the values to the respective columns
  sub_df[row, "protein_name"] <- p_name
  sub_df[row, "protein_synonyms_swiss"] <- syn
  sub_df[row, "gene_synonyms_swiss"] <- g_syn
  
}


# Create a new data frame with only the interesting to us:
Swiss_df <- sub_df %>% 
  dplyr::select(Entry, 
                protein_name, 
                protein_synonyms_swiss,
                gene_synonyms_swiss)


## Tidy up:
rm(
  Uniprot_IDs,
  sub_df,
  string, 
  p_name, 
  syn, 
  ent_name, 
  g_syn,
  i,
  flag,
  row
)

## Double merge, the BioMart_df and Swiss_df, to the input_data_nodes df:
merge_df <- input_data_nodes %>% 
  left_join(., BioMart_df, # First merging the input_data_nodes with BioMart_df
            by=c("UniprotKB" = "uniprotswissprot")) %>% 
  left_join(., Swiss_df, # Then merging the previosly created df with Swiss_df
            by=c("UniprotKB" = "Entry"))


# Create a new column for the final version of gene synonyms,
# and fill it with the function col_splice_to_string():
merge_df$final_gene_syn <- NA
for (i in 1:nrow(merge_df)) {
  if (merge_df[i, "UniprotKB"] != "") {
    merge_df[i, "final_gene_syn"] <- col_splice_to_string(
      df = merge_df,
      row = i,
      col1 = "gene_synonyms_BioMart",
      col2 = "gene_synonyms_swiss",
      from_sep = "|",
      to_sep = "|"
    )
  }
}



#
finished_df <- merge_df %>% 
  mutate(
    Common_name = fill_cells(Common_name, protein_name),
    NCBI_gene = fill_cells(NCBI_gene, entrezgene_id),
    Ensembl = fill_cells(Ensembl, ensembl_gene_id),
    HGNC = fill_cells(HGNC, hgnc_id),
    Origin = ifelse(is.na(Origin), "added", ifelse(Origin == "", "added", Origin)),
    Curator = str_to_title(Curator)
    ) %>% 
  select(
    id,
    Common_name,
    UniprotKB,
    molecule_type,
    "gene_synonyms" = final_gene_syn,
    "protein_synonyms" = protein_synonyms_swiss,
    Receptor_or_Ligand,
    Interleukine.1.signaling,
    TNF.alpha.Signaling,
    JNK.signaling,
    Reactome_ID,
    Signor_ID,
    Ensembl,
    HGNC,
    NCBI_gene,
    Added_from,
    Origin,
    Curator
  )
  



## Tidy up:
rm(merge_df)



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
  
  # Create a data frame with the source, target, and interaction type:
  df <- cbind.data.frame(
    "source"=unlist(source),
    "target"=unlist(target),
    "int"=unlist(int)
  )
  
  return(df)
}



## Takes the information from the BioGateway file and adds genes that are
## not currently in the network:
add_genes_BioGateway <- function(file, updated_df, add_to) {
  input_df <- read.csv(file) %>% 
    filter(Source.Graph != "prot2bp") # Remove the GO-term information
  
  interpreted_input <- read_BioGateway(input_df)
  
  if (add_to == "nodes") {
    genes <- interpreted_input %>% 
      filter(int == "gene") %>% # Only the genes
      filter(target %in% swiss_prot$Entry) %>% # Removes those duplicates (i.e. multiple IDs for a single protein)
      mutate(in_network = ifelse(target %in% updated_df$UniprotKB, T, F)) %>% # Create a new column with T or F depending on protein in network
      filter(in_network == T) %>% # Only keep those with proteins in network
      dplyr::select(-in_network) # Remove the column created earlier
    
    sub_df_from_updated_df <- updated_df %>% 
      filter(UniprotKB %in% genes$target) %>%
      filter(molecule_type != "gene") %>% # 
      mutate(molecule_type = "gene",
             id = str_c(id, "_gene"),
             Curator = "Marius",
             Receptor_or_Ligand = "")
    
    return_df <- bind_rows(updated_df, sub_df_from_updated_df) %>% 
      distinct(.keep_all = T)
    
    return(return_df)
    
  } else if (add_to == "edges") {
    TF_df <- interpreted_input %>% 
      filter(int == "tfac2gene") # Only include rows regarding transcription regulation
    
    unique_genes <- TF_df %>% dplyr::select(target) %>% unique()
    
    df <- getBM(attributes = c(
      "external_gene_name",
      "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
      "external_synonym"), # I.e. gene synonyms
      filters = "external_gene_name",
      values = unique_genes,
      mart = ensembl.con) 
    
    df <- df %>% 
      filter(uniprotswissprot != "") %>% # Remove rows which don't have a Uniprot ID
      dplyr::select(-external_synonym) %>% # Remove the synonyms
      distinct(.keep_all = T) %>% # Remove duplicated rows
      filter(uniprotswissprot %in% updated_df$UniprotKB) # Keep only those that have proteins in network
    
    another_df <- left_join(TF_df, df, by=c("target" = "external_gene_name")) %>% 
      filter(!is.na(uniprotswissprot)) %>% 
      left_join(., updated_df[, c("id", "UniprotKB")], by=c("source" = "UniprotKB")) %>% 
      left_join(., updated_df[, c("id", "UniprotKB")], by=c("uniprotswissprot" = "UniprotKB")) %>% 
      dplyr::select("source" = id.x, 
                    "target" = id.y, 
                    "interaction" = int, 
                    -source, 
                    -target, 
                    -uniprotswissprot) %>% 
      mutate(Curator = "Marius",
             Added_from = "BioGateway",
             target = str_c(target, "_gene"))
    
    
    return_df <- bind_rows(input_data_edges, another_df)
    return(return_df)
  }
}


newest_nodes <- finished_df %>% 
  add_genes_BioGateway("positive regulation of JNK cascade.csv", ., "nodes") %>% 
  add_genes_BioGateway("interleukin-1-mediated signaling pathway.csv", ., "nodes") %>% 
  add_genes_BioGateway("apoptotic process.csv", ., "nodes")


 
newest_edges <- add_genes_BioGateway("positive regulation of JNK cascade.csv", finished_df, "edges") %>% 
  bind_rows(., add_genes_BioGateway("interleukin-1-mediated signaling pathway.csv", finished_df, "edges")) %>% 
  bind_rows(., add_genes_BioGateway("apoptotic process.csv", finished_df, "edges")) %>% 
  distinct(.keep_all = T)



## Adding all edges for genes encoding proteins:
# Finding all genes in the spreadsheet:
genes <- newest_nodes %>% 
  filter(molecule_type == "gene") %>% 
  dplyr::select(UniprotKB) %>% 
  pull() 

# Creating a data frame from the nodes table, and reorganize it to be similar
# to the edge table:
df <- newest_nodes %>%
  filter(UniprotKB %in% genes) %>%
  dplyr::select(id, UniprotKB, molecule_type) %>%
  mutate(molecule_type = ifelse(molecule_type == "gene", "source", "target")) %>%
  group_by(molecule_type) %>% pivot_wider(names_from =  molecule_type, values_from = id) %>% 
  mutate(interaction = "encodes",
         Added_from = "BioGateway",
         Curator = "Marius") %>% 
  dplyr::select(-UniprotKB) %>% 
  relocate(source, target)

# Add the data frame just created to the edge data frame:
newest_edges <- bind_rows(newest_edges, df) %>% 
  distinct(.keep_all = T)


write.csv(newest_nodes, str_c(getwd(), "/new_nodes.csv"), row.names = F)
write.csv(newest_edges, str_c(getwd(), "/new_edges.csv"), row.names = F)




