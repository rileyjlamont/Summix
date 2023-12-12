library(visNetwork)
library(igraph)

#visNetwork pulls from columns w/ specific names
#currently modifying the excel sheet bc lazy but can just append data frame

nodes <- read.csv(("~/Downloads/plot_dat.csv"), header=T, as.is=T)
edges <- read.csv(("~/Downloads/plot_link.csv"), header=T, as.is=T)
colnames(nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(edges) <- c("from", "to", "weight")


nodes$size = 90*nodes$props
edges$width = 3*(3-100*edges$weight)
visNetwork(nodes, edges, width = "100%") %>%
  visNodes(
    shape = "dot",
    shadow = TRUE,
    size = 100,
    font = list(size = 14),
    color = list(
      border = "#013848",
      highlight = "#FF8000"
    ),
  ) %>%
  visEdges(
    shadow = FALSE,
    smooth = FALSE, # remove if want curved edges
    color = list(color = "lightgray", highlight = "#C62F4B")
  )  %>%
  # visLegend(addEdges = data.frame(label = "link"), useGroups = FALSE)
  # %>%
  visLayout(randomSeed = 20)

cont_nodes <- read.csv(("~/Downloads/cont_node.csv"), header=T, as.is=T)
cont_edges <- read.csv(("~/Downloads/cont_link.csv"), header=T, as.is=T)
colnames(cont_nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(cont_edges) <- c("from", "to", "weight")


cont_nodes$size = 70*cont_nodes$props
cont_edges$width = 10*(1.1-6*cont_edges$weight)
visNetwork(cont_nodes, cont_edges, width = "100%") %>%
  visNodes(
    shape = "dot",
    shadow = TRUE,
    size = 100,
    font = list(size = 14),
    color = list(
      border = "#013848",
      highlight = "#FF8000"
    ),
  ) %>%
  visEdges(
    shadow = FALSE,
    smooth = FALSE, # remove if want curved edges
    color = list(color = "lightgray", highlight = "#C62F4B")
  )  %>%
  # visLegend(addEdges = data.frame(label = "link"), useGroups = FALSE)
  # %>%
  visLayout(randomSeed = 20)

unk_nodes <- read.csv(("~/Downloads/unk_dat.csv"), header=T, as.is=T)
unk_edges <- read.csv(("~/Downloads/unk_link.csv"), header=T, as.is=T)
colnames(unk_nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(unk_edges) <- c("from", "to", "weight")


unk_nodes$size = 70*unk_nodes$props
unk_edges$width = 10*(1.1-6*unk_edges$weight)
visNetwork(unk_nodes, unk_edges, width = "100%") %>%
  visNodes(
    shape = "dot",
    shadow = TRUE,
    size = 100,
    font = list(size = 14),
    color = list(
      border = "#013848",
      highlight = "#FF8000"
    ),
  ) %>%
  visEdges(
    shadow = FALSE,
    smooth = FALSE, # remove if want curved edges
    color = list(color = "lightgray", highlight = "#C62F4B")
  )  %>%
  # visLegend(addEdges = data.frame(label = "link"), useGroups = FALSE)
  # %>%
  visLayout(randomSeed = 20)

