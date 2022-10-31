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



################ NB!! ##################
# Set working directory:
#setwd("C:/Users/Legion/Downloads") # Working directory. Change to your directory of choice (e.g. Downloads)
# NB! Short cut for setting working directory: Ctrl+Shift+H or Cmd+Shift+H



### Inputs and first tidy-up ----
# Read the CSVs: OBS! Change "sep" if the CSV file use other separating character, i.e. ","

# Nodes:
nodes <- read.csv("new_nodes.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(nodes)[1] <- "id"
nodes$Common_name <- nodes$Common_name %>% str_replace("\xa0", "") # Fixing a recurring problem

nodes <- remove_whitespaces(nodes)

# Edges:
edges <- read.csv("new_edges.csv", sep = ",", header = T, stringsAsFactors = F)
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
unique(edges$Interaction_result)
edges$Interaction_result <- edges$Interaction_result %>% 
  replace(., . == "deactivation", "inhibits") %>% 
  replace(., . == "inactivation", "inhibits") %>% 
  replace(., . == "up-regulation", "activation")



### Cytoscape ----
# Checks if Cytoscape is open, IF not open Cytoscape first!
cytoscapePing()


### Creates a network in Cytoscape:
createNetworkFromDataFrames(nodes, edges)



my.style <- "group6"
if (my.style %in% getVisualStyleNames() == F) {
  ### Create a custom style: Use: getVisualPropertyNames()
  
  
  
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

setVisualStyle(my.style)
lockNodeDimensions(FALSE, my.style)
