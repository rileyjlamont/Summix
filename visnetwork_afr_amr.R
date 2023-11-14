library(visNetwork)
library(igraph)


my_directory = "~/Users/rileylamont/Downloads/"
### AFRICAN
afr_nodes <- read.csv(("~/Downloads/afr_fs_dat.csv"), header=T, as.is=T)
afr_edges <- read.csv(("~/Downloads/afr_link_dat.csv"), header=T, as.is=T)

#visNetwork pulls from columns w/ specific names
#currently modifying the excel sheet bc lazy but can just append data frame
colnames(afr_nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(afr_edges) <- c("from", "to", "weight")

#pal <- c('AFR' = "#FDE725FF", 'EAS' = "#5DC863FF", 'EUR' = "#21908CFF", 'IAM' = "#3B528BFF", 'SAS' = "#440154FF")
#56B4E9'
#scale
afr_nodes$size = 150*afr_nodes$props
afr_edges$width = 1250*afr_edges$weight

visNetwork(afr_nodes, afr_edges, width = "100%") %>%
  visNodes(
    shape = "dot",
    shadow = TRUE,
    size = 100,
    font = list(size = 22),
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
  visLayout(randomSeed = 97)

### ADMIXED AMERICAN
amr_nodes <- read.csv(("~/Downloads/amr_fs_dat.csv"), header=T, as.is=T)
amr_edges <- read.csv(("~/Downloads/amr_link_dat.csv"), header=T, as.is=T)

colnames(amr_nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(amr_edges) <- c("from", "to", "weight")

#pal <- c('AFR' = "#FDE725FF", 'EAS' = "#5DC863FF", 'EUR' = "#21908CFF", 'IAM' = "#3B528BFF", 'SAS' = "#440154FF")


amr_nodes$size = 200*amr_nodes$props
amr_edges$width = 1250*amr_edges$weight

visNetwork(amr_nodes, amr_edges, width = "100%") %>%
  visNodes(
    shape = "dot",
    shadow = TRUE,
    size = 100,
    font = list(size = 22),
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
  visLayout(randomSeed = 8)

# drag italian point out of cluster
# figure out legend

nodes <- read.csv(("~/Downloads/hayley_plot_dat.csv"), header=T, as.is=T)
edges <- read.csv(("~/Downloads/hayley_plot_link.csv"), header=T, as.is=T)
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

nodes <- read.csv(("~/Downloads/adelle_plot_dat.csv"), header=T, as.is=T)
edges <- read.csv(("~/Downloads/adelle_plot_link.csv"), header=T, as.is=T)
colnames(nodes) <- c("id", "props", "ancs", "group", "label", "color")
colnames(edges) <- c("from", "to", "weight")


nodes$size = 90*nodes$props
edges$width = 3*(3-300*edges$weight)
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

