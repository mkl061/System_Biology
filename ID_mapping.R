### Package handling ----

# Bioconductor (i.e. BiocManager):
if (!require("BiocManager", quietly = T)) {
  install.packages("BiocManager", quiet = T)
}
#BiocManager::install(update = T, ask = F, ... = c(quiet = T))
# OBS! The command above usually give a warning, but apparently this
# can be ignored. See: https://support.bioconductor.org/p/128825/ 



if (!require("biomaRt", quietly = T)) {
  #install.packages("biomaRt", quiet = T)
  BiocManager::install("biomaRt")
}



# Tidyverse:
if (!require("tidyverse", quietly = T)) {
  install.packages("tidyverse", quiet = T)
}

# Readxl (NB! Is part of tidyverse, but needs to be loaded separately)
if (!require("readxl", quietly = T)) {
  install.packages("readxl", quiet = T)
}



### Working directory ----
if (Sys.info()["nodename"] == "MARIUSPC") { # Marius 1
  setwd("C:/Users/Legion/Downloads")
} else if (Sys.info()["nodename"] == "penguin") { # Marius 2
  setwd("/mnt/chromeos/MyFiles/Downloads")
} else if (Sys.info()["nodename"] == "dhcp-10-22-28-162.wlan.ntnu.no") { # Jann
  setwd("/Users/Jann/Library/CloudStorage/OneDrive-SharedLibraries-NTNU/o365_System Bio - group 6 - Documents/General")
} else if (Sys.info()["nodename"] == "LAPTOP-3BEG4HPK") { # Vibeke
  setwd("C:/Users/Vibeke/Downloads")
}


### Custom function ----

# Function that looks through all character columns, and removes white 
# spaces at front and end, and convert multiple spaces to single spaces:
remove_whitespaces <- function(df) {
  df <- df %>% 
    mutate( # Change column values
      across( # All columns which ...
        where(is.character), # ... are of type character/string
        str_squish # Apply this function. Which removes whitespaces
      )
    )
  
  return(df) # Return the data frame
}



# Function to reverse all elements in a character column:
str_rev_in_vec <- function(vec) {
  vec <- vec %>% 
    strsplit(NULL) %>% 
    lapply(rev) %>% 
    lapply(paste, collapse = "")
  
  return(vec)
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
files <- dir(pattern = ".xlsx")
files <- files[str_starts(files, "Project_workbook")]
mtimes <- file.mtime(files) %>% order()
file <- files[mtimes][length(files)]

input_data_nodes <- read_excel(file, sheet = "nodes")

#input_data_nodes$Common_name <- input_data_nodes$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

# Tidy up the file by removing unwanted white spaces:
input_data_nodes <- remove_whitespaces(input_data_nodes)

# Fix capital letter in name:
input_data_nodes <- input_data_nodes %>% 
  mutate(Curator = str_to_title(Curator))


## Edges:

input_data_edges <- read_excel(file, sheet = "edges")

input_data_edges <- remove_whitespaces(input_data_edges)

# Fix capital letter in name:
input_data_edges <- input_data_edges %>% 
  mutate(Curator = str_to_title(Curator))


# Tidy up:
rm(files, file, mtimes)



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
# Moreover, you can't get all identifiers in a single step,
# so here we create two data frames, which we merge later on.
multiple_ids1 <- getBM(
  attributes = c(
    "uniprotswissprot", # Uniprot IDs - Swiss-Prot (i.e. Reviewed content)
    "entrezgene_id",  # Entrez (i.e. NCBI) gene ID/number
    "external_gene_name", # Common gene name
    "ensembl_gene_id" # Ensemble gene ID
    ), 
  filters = "uniprotswissprot", # Map using these IDs (here: Uniprot IDs)
  values = input_data_nodes$UniprotKB, # The column the IDs are found
  mart = ensembl.con) %>% 
  mutate(dup = duplicated(uniprotswissprot)) %>% # OBS! Because Ensemble have multiple IDs per protein, we choose to only keep one
  filter(dup != T) %>% 
  dplyr::select(-dup)


multiple_ids2 <- getBM(
  attributes = c(
    "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
    "hgnc_id", # HGNC ID
    "hgnc_symbol"), # I.e. Gene symbol as reported in HGNC
  filters = "uniprotswissprot",
  values = input_data_nodes$UniprotKB,
  mart = ensembl.con)


# Merge the two data frames into one data frame:
multiple_ids <- left_join(multiple_ids1, multiple_ids2, 
                          by=c("uniprotswissprot" = "uniprotswissprot"))

# Tidy up:
rm(multiple_ids1, multiple_ids2)


# Creates a data frame with all the synonyms:
# Each synonym has its own row, meaning that there are multiple rows
# with the same id/value in the "uniprotswissprot" column.
raw_gene_syn_df <- getBM(
  attributes = c(
    "uniprotswissprot", # Uniprot ID, from the Swiss-Prot data base
    "external_synonym"), # I.e. gene synonyms
  filters = "uniprotswissprot",
  values = input_data_nodes$UniprotKB,
  mart = ensembl.con)


# Use the raw_gene_syn_df to create another data frame with only
# one row per Uniprot ID, and all the synonyms combined in a single 
# string per Uniprot ID:
fin_gene_syn_df <- raw_gene_syn_df %>% 
  group_by(uniprotswissprot) %>% 
  summarise(gene_synonyms_BioMart = 
              paste0(external_synonym, 
                     collapse = "|"))


# Combine the data frame with different IDs with the the data frame
# with all the gene synonyms:
BioMart_df <- left_join(
  multiple_ids, fin_gene_syn_df, 
  by=c("uniprotswissprot" = "uniprotswissprot"))


## Tidying up:
rm(multiple_ids,
   raw_gene_syn_df, 
   fin_gene_syn_df
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
  # Download the file to the directory:
  write_tsv(swiss_prot, str_c(getwd(), "/Swiss_Prot.tsv"))
}


# All the Uniprot IDs from the spreadsheet:
Uniprot_IDs <- input_data_nodes %>%
  filter(UniprotKB != "") %>% 
  dplyr::select(UniprotKB) %>% 
  pull() # Results in character vector



# The strings may contain information about what the result of
# protein cleavage is. This function determine what should be done
# with that information:
cleavage_info <- function(string, remove_or_extract) {
  if (str_detect(string, "\\[Cleaved into: ")) {
    ind_start <- str_locate(string, "\\[Cleaved into: ")[1]
    ind_end <- str_locate(string, "\\]")[2]
    
    
    string_first <- str_sub(string, end = ind_start-1)
    string_end <- str_sub(string, start = ind_end+1)
    
    if (remove_or_extract == "remove") {
      string <- str_c(string_first, string_end) %>% str_squish()
      return(string)
    } else if (remove_or_extract == "extract") {
      string <- str_sub(string, ind_start, ind_end) %>% str_squish()
      return(string)
    }
  }
  
  if (remove_or_extract == "remove") {
    return(string)
  } else {
    return("")
  }
}



# Convert from "name1 (name2) (name3)" to "name1|name2|name3":
str_change_sep <- function(string) {
  # Checks if sting include a left parenthesis, which would
  # mean that the string consists of multiple terms:
  if (str_detect(string, "\\(")) {
    # Checks if string includes ") (", which would mean that
    # the string consists of at least 3 terms:
    if (str_detect(string, "\\) \\(")) {
      string <- string %>% 
        str_replace_all( # Replace ") (" with "|"
          pattern = "\\) \\(",
          replacement = "|") 
    }
    
    # If only 2 terms in string, or after the ") (" are replaced
    # in the previous IF-statement:
    string <- string %>% 
      str_replace( # Replace the " (" with "|"
        pattern = " \\(",
        replacement = "|") 
    
    # The only thing missing is the last parenthesis, if present:
    if (str_ends(string, "\\)")) {
      string <- str_sub(string, end = -2)
    }
  }
  
  # Return the string:
  # If the string only consisted of a single term, all the if-else statements
  # would have been skipped, and the string is unaltered.
  return(string)
}



# Function to remove duplicates in a string:
#str_rm_duplicates <- function(string, sep, collapse) {
str_rm_duplicates <- function(string, sep, collapse) {
  # Handle situations where there is only a single term, or 
  # the "string" is NA:
  if (is.na(string)) {
    return("")
  } else if (sep == "|" && str_detect(string, "\\|") == F) {
    return(string)
  } else if (sep != "|" && str_detect(string, sep) == F) {
    return(string)
  }
  
  
  if (sep == "|") {
    sep_bar <- T
  } else {
    sep_bar <- F
  }
  if (collapse == "|") {
    col_bar <- T
  } else {
    col_bar <- F
  }
  
  
  if (all(sep_bar, col_bar)) {
    string <- string %>% 
      str_split(pattern = "\\|") %>% 
      unlist() %>% 
      base::unique() %>% 
      str_c(collapse = "|")
  } else if (any(sep_bar, col_bar)) {
    if (sep_bar) {
      string <- string %>% 
        str_split(pattern = "\\|") %>% 
        unlist() %>% 
        base::unique() %>% 
        str_c(collapse = "|")
    } else {
      string <- string %>% 
        str_split(pattern = sep) %>% 
        unlist() %>% 
        base::unique() %>% 
        str_c(collapse = "|")
    }
  } else {
    string <- string %>% 
      str_split(pattern = sep) %>% 
      unlist() %>% 
      base::unique() %>% 
      str_c(collapse = collapse)
  }
  
  return(string)
}

# Vectorize the function str_rm_duplicates:
str_rm_duplicates_vec <- Vectorize(str_rm_duplicates, "string")



# Function to remove the EC-number which may be part of the string:
str_rm_ECnumber <- function(string) {
  string <- string %>% 
    str_split(pattern = "\\|") %>% # Splits the string at "|", into list
    unlist()
  
  # Does at least one element starts with "EC "?:
  if (any(str_detect(string, "EC "))) {
    index <- str_which(string, pattern = "EC ") # Locate the elements
    
    string <- string[-index] # Removes by index values
  }
  
  string <- str_c(string, collapse = "|") # Collapse vector to string
  return(string)
}


# Subset of swiss_prot df, with only IDs equal to those in our input data:
Swiss_df <- swiss_prot %>% 
  filter(Entry %in% Uniprot_IDs) %>% 
  group_by(Entry) %>% 
  mutate( # Change the content of columns:
    # Extract the cleavage information:
    Cleaved = cleavage_info(Protein.names, "extract"),
    # Remove the clevage information from the column Protein.names:
    Protein.names = cleavage_info(
      Protein.names,
      remove_or_extract = "remove"),
    # Change the separator of terms in the column Protein.names:
    Protein.names = str_change_sep(Protein.names),
    # Remove the EC-numbers in the column Protein.names:
    Protein.names = str_rm_ECnumber(Protein.names),
    # Construct a new column with only the first term from Protein.names:
    Protein_name =
      ifelse(
        str_detect(Protein.names, "\\|"),
        str_sub(
          Protein.names,
          end = str_locate(Protein.names, "\\|")[1]-1
        ),
        Protein.names
      ),
    Protein_synonyms =
      ifelse(
        str_detect(Protein.names, "\\|"), # IF statement
        str_sub(Protein.names, # True; cell = only protein synonyms (i.e. without protein name)
                start = str_locate(
                  Protein.names,
                  pattern = "\\|")[1]+1),
        ""), # False; cell = empty (i.e. no protein synonyms)
    # Add the Entry name as a synonym:
    Protein_synonyms = str_c(Entry.Name, Protein_synonyms, sep="|"),
    # Remove posible "|" at end of string:
    Protein_synonyms = ifelse(str_ends(Protein_synonyms,
                                       pattern = "\\|"),
                              str_sub(Protein_synonyms, end = -2),
                              Protein_synonyms),
    # Construct a new column of gene synonyms from the columns Gene.Names 
    # and Gene.Names..synonym.:
    Gene_synonyms = str_c(Gene.Names, Gene.Names..synonym.,
                          sep = " "),
    # Remove duplicated synonyms:
    Gene_synonyms = str_rm_duplicates_vec(Gene_synonyms, 
                                      sep = " ",
                                      collapse = "|"),
    # Remove posible "|" at end of string:
    Gene_synonyms = ifelse(str_ends(Gene_synonyms,
                                    pattern = "\\|"),
                           str_sub(Gene_synonyms, end = -2),
                           Gene_synonyms)
  ) %>% 
  dplyr::select(Entry, 
                "protein_name"=Protein_name, 
                "protein_synonyms_swiss"=Protein_synonyms, 
                Cleaved, 
                "gene_synonyms_swiss"=Gene_synonyms)






## Double merge, the BioMart_df and Swiss_df, to the input_data_nodes df:
merge_df <- input_data_nodes %>% 
  left_join(., BioMart_df, # First merging the input_data_nodes with BioMart_df
            by=c("UniprotKB" = "uniprotswissprot")) %>% 
  left_join(., Swiss_df, # Then merging the previously created df with Swiss_df
            by=c("UniprotKB" = "Entry"))



merge_df <- merge_df %>% 
  mutate(gene_synonyms = 
           str_c(gene_synonyms_BioMart, 
                 gene_synonyms_swiss,
                 sep = "|"),
         gene_synonyms = ifelse(
           str_starts(gene_synonyms, "\\|"),
           str_sub(gene_synonyms, 2),
           ifelse(
             str_ends(gene_synonyms, "\\|"),
             str_sub(gene_synonyms, -2),
             gene_synonyms
           )
         ),
         gene_synonyms = 
           str_rm_duplicates_vec(gene_synonyms,
                             sep = "|",
                             collapse = "|"))

# Tidy up:
rm(Swiss_df,
   BioMart_df)




#
finished_df <- merge_df %>% 
  mutate(
    Common_name = fill_cells(Common_name, protein_name),
    NCBI_gene = fill_cells(NCBI_gene, entrezgene_id),
    Ensembl = fill_cells(Ensembl, ensembl_gene_id),
    HGNC = fill_cells(HGNC, hgnc_id),
    Origin = ifelse(is.na(Origin), "added", ifelse(Origin == "", "added", Origin))
    ) %>% 
  dplyr::select(
    id,
    Common_name,
    UniprotKB,
    hgnc_symbol,
    molecule_type,
    gene_synonyms,
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
             Receptor_or_Ligand = "",
             Origin = "added",
             Reactome_ID = NA,
             Signor_ID = NA)
    
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

# Creating, or updating, new CSV files:
write.csv(newest_nodes, str_c(getwd(), "/new_nodes.csv"), row.names = F)
write.csv(newest_edges, str_c(getwd(), "/new_edges.csv"), row.names = F)
