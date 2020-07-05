# This script contains codes for generating Fig. 2

library(igraph)

# Load data matrix containing the hazard ratio of acquiring an ARO in association to prior exposure (past 30 days) to antibiotics or pre-existing ARO colonization. 

final_mat = readRDS("2020-07-04_org_abx_mat.RData")

tip_org_abx_graph_p0.05 = graph.adjacency(t(final_mat), weighted = TRUE, mode = "directed")
E(tip_org_abx_graph_p0.05)$width = E(tip_org_abx_graph_p0.05)$weight*0.5
V(tip_org_abx_graph_p0.05)$name = gsub("_", " ", V(tip_org_abx_graph_p0.05)$name)

tip_org_abx_graph_p0.05 = graph.adjacency(t(org_abx30_graph_mat_p0.05), weighted = TRUE, mode = "directed")
E(tip_org_abx_graph_p0.05)$width = E(tip_org_abx_graph_p0.05)$weight*0.5
V(tip_org_abx_graph_p0.05)$name = gsub("_", " ", V(tip_org_abx_graph_p0.05)$name)


V(tip_org_abx_graph_p0.05)$color <- ifelse(V(tip_org_abx_graph_p0.05)$name %in% names(abx_groups), "gray80", "limegreen") # abx_groups defined in Fig2_data.R

# Plot interaction network
# First plot with tkplot to set up coordinates for future use
org_abx_graph = tkplot(tip_org_abx_graph_p0.05, 
                       vertex.label.font = 2, 
                       vertex.label.color = "black", 
                       vertex.frame.color =NA,
                       edge.label = round(E(tip_org_abx_graph_p0.05)$weight,1), 
                       edge.label.cex = 1.25, 
                       edge.label.color = 'red', 
                       edge.label.font = 2, 
                       vertex.label.cex = 1.25, 
                       vertex.size = 30, 
                       edge.color = "black")

org_abx_graph_coords = tk_coords(org_abx_graph, norm = FALSE)

# write.table(org_abx_graph_coords, paste0("org_abx_network_coords.txt'))

# org_abx_graph_coords = read.table('_org_abx_network_coords.txt')

# To save, create a filename
pdf(filename)

plot(tip_org_abx_graph_p0.05,
vertex.label.font = 2,
vertex.label.color = "black",
vertex.frame.color =NA,
edge.label = round(E(tip_org_abx_graph_p0.05)$weight,1),
edge.label.cex = 1.25,
edge.label.color = 'red',
edge.label.font = 2,
vertex.label.cex = 1.25,
vertex.size = 15,
edge.color = "black",
layout = as.matrix(org_abx_graph_coords))

dev.off()