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
nodes <- read.csv("nodes.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(nodes)[1] <- "id"
nodes$Common_name <- nodes$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

nodes <- remove_whitespaces(nodes)

# Edges:
edges <- read.csv("edges.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(edges)[1] <- "source"
colnames(edges)[2] <- "target"

edges <- remove_whitespaces(edges)



### Adding missing nodes, and standardizations ----
# Check for uniqueness:
n_un <- unique(nodes$id)
e_un <- unique(c(edges$source, edges$target))
dif <- setdiff(e_un, n_un) # Contain elements of edges (source or target) which is not part of nodes
# Ads the missing nodes:
for (i in 1:length(dif)) {
  nodes <- rbind(nodes,
                 c(dif[i], rep("", times=ncol(nodes)-1))
                 )
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


### Create a custom style:
my.style <- "group6"
copyVisualStyle("default", my.style)
lockNodeDimensions(FALSE, my.style)


# Mapping:
unique(nodes$molecule_type) # Checks what the categories are
setNodeColorMapping("molecule_type", 
                    style.name = my.style,
                    table.column.values = c("complex", "protein"), 
                    mapping.type = "d", 
                    colors = c("#789ef5", "#a0b9f2"))

setNodeShapeMapping("molecule_type", 
                    style.name = my.style, 
                    table.column.values = c("protein", "complex"),
                    shapes = c("ROUND_RECTANGLE", "HEXAGON"))

setEdgeTargetArrowMapping("Interaction_result", 
                          style.name = my.style, 
                          table.column.values = c("inhibits", "activation", ""), 
                          shapes = c("T", "ARROW", "NONE"))

setEdgeTargetArrowColorMapping("Interaction_result", 
                               style.name = my.style, 
                               table.column.values = c(
                                 "inhibits", 
                                 "activation", 
                                 ""), 
                               mapping.type = "d", 
                               colors = c("#a1453d", "#5ea63f", "#d6d6d6"))

setEdgeColorMapping("Interaction_result", 
                    style.name = my.style, 
                    table.column.values = c("inhibits", "activation", ""), 
                    mapping.type = "d", 
                    colors = c("#a1453d", "#5ea63f", "#d6d6d6"))

