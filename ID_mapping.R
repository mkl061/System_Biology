library(org.Hs.eg.db)

# API URL for ALL human proteins in the Swiss-Prot database:
# https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names&format=tsv&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29

swiss_prot <- read.delim("Swiss.tsv")

## Your spreadsheet:
#dat <- read.csv("original_set.csv")
dat <- read.csv("nodes.csv")
dat$Common_name <- dat$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

## Trim all whitespaces in all cells with strings:
for (col_ind in 1:ncol(dat)) {
  if (is.character(dat[,col_ind])) {
    dat[,col_ind] <- sapply(dat[,col_ind], str_trim)
  }
}

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

dat2 <- left_join(dat, merg_table)
