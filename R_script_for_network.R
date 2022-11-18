### Packages ----

# Bioconductor (i.e. BiocManager):
if (!require("BiocManager", quietly = T)) {
  install.packages("BiocManager", quiet = T)
}



# RCy3:
if (!require("RCy3", quietly = T)) {
  BiocManager::install("RCy3")
}


# Tidyverse:
if (!require("tidyverse", quietly = T)) {
  install.packages("tidyverse", quiet = T)
}



### Functions: ----
# Function that looks through all character columns, and removes
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



################ NB!! ##################
# Set working directory:
#setwd("C:/Users/Legion/Downloads") # Working directory. Change to your directory of choice (e.g. Downloads)
# NB! Short cut for setting working directory: Ctrl+Shift+H or Cmd+Shift+H
if (Sys.info()["nodename"] == "MARIUSPC") { # Marius 1
  setwd("C:/Users/Legion/Downloads")
} else if (Sys.info()["nodename"] == "penguin") { # Marius 2
  setwd("/mnt/chromeos/MyFiles/Downloads")
} else if (Sys.info()["nodename"] == "dhcp-10-22-28-162.wlan.ntnu.no") { # Jann
  setwd("~/Desktop/run everything")
} else if (Sys.info()["nodename"] == "LAPTOP-3BEG4HPK") { # Vibeke
  setwd("C:/Users/Vibeke/Downloads")
}


#### Omics data ----
# GSE14905:
GSE14905 <- read.csv("GSE14905.csv")[,1:3] # The file should only contain 3 columns
# COAD_CMS2:
COAD_CMS2 <- read.csv("COAD_CMS2.csv")[,1:3]
# TCGA_BRCA:
TCGA_BRCA <- read.csv("TCGA_BRCA.csv")[,1:3]
# GSE120245:
GSE120245 <- read.csv("GSE120245.csv")[,1:3]




### Inputs and first tidy-up ----
# Read the CSVs: OBS! Change "sep" if the CSV file use other separating character, i.e. ","

## Nodes:
nodes <- read.csv("new_nodes.csv", sep = ",", header = T, stringsAsFactors = F)

# General tidying up:
colnames(nodes)[1] <- "id"
nodes$Common_name <- nodes$Common_name %>% str_replace("\xa0", "") # Fixing a recurring problem

nodes <- remove_whitespaces(nodes) # Remove all whitespaces where there should be none


# Incorporate omics data:
inc_omics <- function(df, omics_df) {
  # Edit the column names of the omics df to identify the source of the
  # omics data:
  
  df.name <- deparse(substitute(omics_df))
  colnames(omics_df)[2] <- str_c("logFC_", df.name)
  colnames(omics_df)[3] <- str_c("FDR_", df.name)
  
  # Join the data frames
  df <- df %>% 
    left_join(., omics_df, by=c("hgnc_symbol" = "hgnc_symbol"))
  
  return(df)
}


# Use the function to merge the dfs
nodes <- nodes %>% 
  inc_omics(GSE14905) %>% 
  inc_omics(COAD_CMS2) %>% 
  inc_omics(TCGA_BRCA) %>% 
  inc_omics(GSE120245) 


# Obs! The RCy3 package and Cytoscape struggles with changing the column
# type (e.g. string -> numeric). The following lines of code will create a
# dummy row with values just to force the column type. The mentioned row
# will be removed after the network is created, by other lines of code.
test_row <- rep(NA, ncol(nodes))
test_row[1] <- 12345

col_index <- colnames(nodes)[
  str_detect(colnames(nodes), "logFC") | str_detect(colnames(nodes), "FDR")
  ] %>% 
  sapply(., grep, colnames(nodes)) %>% 
  as.vector()

test_row[col_index] <- 1.1

nodes <- rbind(test_row, nodes[1:nrow(nodes), ])



## Edges:
edges <- read.csv("new_edges.csv", sep = ",", header = T, stringsAsFactors = F)

# General tyding up: 
colnames(edges)[1] <- "source"
colnames(edges)[2] <- "target"

edges <- remove_whitespaces(edges)



### Adding missing nodes, and standardizations ----
# Check for uniqueness:
n_un <- unique(nodes$id)
e_un <- unique(c(edges$source, edges$target))
dif <- setdiff(e_un, n_un) # Contain elements of edges (source or target) which is not part of nodes
# Ads the missing nodes:
if (length(dif) > 0) {
  for (i in 1:length(dif)) {
    nodes <- rbind(nodes,
                   c(dif[i], rep("", times=ncol(nodes)-1))
    )
  }  
}



# Standardize activation/inhibition:
# unique(edges$Interaction_result)
edges$Interaction_result <- edges$Interaction_result %>%
  replace(., . == "deactivation", "inhibits") %>%
  replace(., . == "inactivation", "inhibits")



### Cytoscape ----
# Checks if Cytoscape is open, IF not open Cytoscape first!
cytoscapePing()


### Creates a network in Cytoscape:
today <- Sys.Date() %>% # Finds today's date
  format(format="%d.%m.%Y") # Date on format: dd.mm.yyyy

# The network creation code-line:
createNetworkFromDataFrames(nodes, edges, 
                            str_c("NETWORK GROUP 6, updated: ", today)
                            )

# Remove the dummy-row:
selectNodes("12345", "id")
deleteSelectedNodes()


# Generate a custom style:
my.style <- "Group6_style"
if (my.style %in% getVisualStyleNames() == F) {
  
  defaults <- list(NODE_SHAPE="ROUND_RECTANGLE",
                   NODE_HEIGHT = 35,
                   NODE_WIDTH = 75,
                   EDGE_TRANSPARENCY=255,
                   EDGE_WIDTH = 1.5)
  
  nodeLabels <- mapVisualProperty('node label',
                                  'id',
                                  'p')
  
  nodeShape <- mapVisualProperty("NODE_SHAPE",
                                 'molecule_type',
                                 'd',
                                 c("protein", "complex", "TF", "gene"), 
                                 c("ROUND_RECTANGLE", "HEXAGON", "ROUND_RECTANGLE", "ELLIPSE"))
  
  nodeFills <- mapVisualProperty('node fill color',
                                 'molecule_type',
                                 'd',
                                 c("protein", "complex", "TF", "gene"), 
                                 c("#a0b9f2","#789ef5", "#ff5a36", "#ffff99"))
  
  arrowShapes <- mapVisualProperty("EDGE_TARGET_ARROW_SHAPE",
                                   'Interaction_result',
                                   'd',
                                   c("activation","inhibits","interacts", ""),
                                   c("Arrow","T","None", "None"))
  
  arrowType <- mapVisualProperty("EDGE_LINE_TYPE",
                                 "interaction",
                                 "d",
                                 c("pp", "tfac2gene", "encodes"),
                                 c("SOLID", "VERTICAL_SLASH", "SEPARATE_ARROW"))
  
  arrowColor <- mapVisualProperty("EDGE_STROKE_UNSELECTED_PAINT",
                                  "Interaction_result",
                                  "d",
                                  c("activation", "inhibits", "interacts", ""),
                                  c("#1fb805", "#d60000", "#4e4e4e", "#4e4e4e"))
  
  arrowTargetColor <- mapVisualProperty("EDGE_TARGET_ARROW_UNSELECTED_PAINT",
                                        "Interaction_result",
                                        "d",
                                        c("activation","inhibits","interacts", ""),
                                        c("#1fb805", "#d60000", "#4e4e4e", "#4e4e4e"))
  
  edgeWidth <- mapVisualProperty('edge width', 
                                 "interaction",
                                 'd', 
                                 table.column.values = c("pp"), 
                                 visual.prop.values = c(2))
  
  
  createVisualStyle(my.style, defaults, list(nodeLabels,
                                             nodeShape,
                                             nodeFills,
                                             arrowShapes,
                                             arrowType,
                                             arrowColor,
                                             arrowTargetColor,
                                             edgeWidth))
  
}

# Apply the visual style:
setVisualStyle(my.style)
lockNodeDimensions(FALSE, my.style)


# Create a style for interpreting expression data:
my.style2 <- "Group6_omics_style"
copyVisualStyle(my.style, my.style2)
setNodeColorDefault("#bcbcbc", my.style2)
setNodeShapeMapping(table.column = "molecule_type",
                    table.column.values = c("protein", "TF", "gene"),
                    shapes = c("ROUND_RECTANGLE", "TRIANGLE", "ELLIPSE"),
                    style.name = my.style2)
setNodeColorMapping(table.column = "logFC_GSE14905",
                          style.name = my.style2)


# Saves the file, with overwrite previalages:
saveSession("Updated_network")
