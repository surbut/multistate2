library(visNetwork)

# Nodes data with reduced label text
nodes <- data.frame(
  id = 1:10, 
  label = c("Health", "Diabetes", "Hypertension", "HyperLipid", 
            "D + Hypertension", "Hypertension + HyperLipid", "HyperLipid + D", 
            "Triple", "CAD", "Death"),
  level = c(1, 2, 2, 2, 3, 3, 3, 4, 5, 5),
  color = c("lightblue", "lightcoral", "lightcoral", "lightcoral", 
            "lightgreen", "lightgreen", "lightgreen", "lightyellow", "grey", "black"),
  font.size = 25,  # adjusted font size
  size = 25  # adjusted node size
)

# Edges data
edges <- data.frame(
  from = c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9),
  to = c(2, 3, 4, 9, 10, 5, 9, 10, 6, 9, 10, 7, 9, 10, 8, 9, 10, 8, 9, 10, 8, 9, 10, 9, 10, 10)
)

# Plot
visNetwork(nodes, edges, width = "100%", height = "1000px") %>%
  visLayout(hierarchical = list(enabled = TRUE, direction = "UD", levelSeparation = 300)) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(
    hierarchicalRepulsion = list(centralGravity = 0.0, springLength = 300, springConstant = 0.05, nodeDistance = 350)
  )


