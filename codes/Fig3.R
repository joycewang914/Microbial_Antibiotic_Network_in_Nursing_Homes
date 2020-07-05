# This script contains codes for generating the "backbone" of Fig 3 - HR of polymicrobial interaction in CAUTI was added using Adobe Illustrator. 

library(igraph)

# Load data matrix containing the hazard ratio of mono- and co-colonization of AROs associated with developing a CAUTI

outcome_network_matrix = readRDS("2020-07-04_org_cauti_HR_mat.RData")

outcome_network_matrix_graph = graph.adjacency(t(outcome_network_matrix), weighted = TRUE, mode = "directed")
E(outcome_network_matrix_graph)$width = E(outcome_network_matrix_graph)$weight*1

V(outcome_network_matrix_graph)$name = gsub("_", " ", V(outcome_network_matrix_graph)$name)

V(outcome_network_matrix_graph)$color <- ifelse(grepl("CAUTI", V(outcome_network_matrix_graph)$name),"brown1", "blue")


plot.igraph(outcome_network_matrix_graph , vertex.label.font = 2, vertex.label.color = "black", edge.label = round(E(outcome_network_matrix_graph)$weight,2), edge.label.cex = 1, edge.label.color = 'red', edge.label.font = 2, vertex.label.cex = 1, edge.color = "black")

