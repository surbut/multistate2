# Nodes data with adjusted levels and colors
nodes <- data.frame(
  id = 1:8, 
  label = c("Health", "Diabetes", "Hypertension", "Hyperlipidemia", "Double risk factor", 
            "Triple risk factor", "CAD", "Death"),
  level = c(1, 2, 2, 2, 3, 3, 4, 4),
  color = c("lightblue", "lightcoral", "lightcoral", "lightcoral", "lightgreen", "lightyellow", "grey", "black")
)

# Edges data
edges <- data.frame(
  from = c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7),
  to = c(2, 3, 4, 7, 8, 5, 7, 8, 5, 7, 8, 5, 7, 8, 6, 7, 8, 7, 8, 8)
)

# Plot
visNetwork(nodes, edges, width = "100%") %>%
  visLayout(hierarchical = list(enabled = TRUE, direction = "UD", levelSeparation = 150)) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
