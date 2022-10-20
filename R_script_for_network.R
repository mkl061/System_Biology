is_installed <- function(package_id) {
  pacs <- installed.packages()[ , "Package"]
  if (!(package_id %in% pacs)) {
    install.packages(package_id)
  } else {
    cat("The package is installed")
    return(T)
  }
}

# Checks if the packages are installed, and installs them if not:
is_installed("tidyverse")
is_installed("BiocManager")

if (is_installed("RCy3")) {
} else {
  BiocManager::install("RCy3")
}

# Load the packages:
library(tidyverse)
library(RCy3)

# Set working directory:
#setwd("C:/Users/Legion/Downloads") # Working directory. Change to your directory of choice (e.g. Downloads)
# NB! Short cut for setting working directory: Ctrl+Shift+H or Cmd+Shift+H


# Read the CSVs: OBS! Change "sep" if the CSV file use other separating character, i.e. ","
nodes <- read.csv("nodes.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(nodes)[1] <- "id"
nodes$Common_name <- nodes$Common_name %>% str_replace("\xa0", "") # Fixing a reocurring problem

edges <- read.csv("edges.csv", sep = ",", header = T, stringsAsFactors = F)
colnames(edges)[1] <- "source"
colnames(edges)[2] <- "target"


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

###
#nodes <- nodes %>% select(id)

# Standardize activation/inhibition:
unique(edges$Interaction_result)
edges$Interaction_result <- edges$Interaction_result %>% 
  replace(., . == "deactivation", "inhibits") %>% 
  replace(., . == "inactivation", "inhibits") %>% 
  replace(., . == "up-regulation", "activation")  


# Checks if Cytoscape is open, IF not open Cytoscape first!
cytoscapePing()


# Creates a network in Cytoscape:
createNetworkFromDataFrames(nodes, edges)


# Create a custom style:
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

