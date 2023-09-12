# Nodes (representing states)
nodes <- data.frame(id = c("Health", "HyperLip", "Hypertension", "Diabetes", 
                           "HyperLip & Hypertension", "Diabetes & Hypertension", 
                           "Diabetes & HyperLip", "Triple Risk", "CAD", "Death"),
                    level = c(1, 2, 2, 2, 3, 3, 3, 4, 5, 6),
                    group = c(1, 2, 2, 2, 3, 3, 3, 4, 5, 5))

# Edges (representing transitions)
edges <- data.frame(from = c("Health", "Health", "Health", "Health", "Health",
                             "HyperLip", "HyperLip", "HyperLip", "HyperLip", "HyperLip",
                             "Hypertension", "Hypertension", "Hypertension", "Hypertension", "Hypertension",
                             "Diabetes", "Diabetes", "Diabetes", "Diabetes", "Diabetes",
                             "HyperLip & Hypertension", "HyperLip & Hypertension", "HyperLip & Hypertension",
                             "Diabetes & Hypertension", "Diabetes & Hypertension", "Diabetes & Hypertension",
                             "Diabetes & HyperLip", "Diabetes & HyperLip", "Diabetes & HyperLip",
                             "Triple Risk", "Triple Risk", "Triple Risk",
                             "CAD"),
                      to = c("HyperLip", "Hypertension", "Diabetes", "CAD", "Death",
                             "HyperLip & Hypertension", "Diabetes & HyperLip", "HyperLip & Hypertension", "CAD", "Death",
                             "HyperLip & Hypertension", "Diabetes & Hypertension", "Diabetes & Hypertension", "CAD", "Death",
                             "Diabetes & HyperLip", "Diabetes & Hypertension", "Diabetes & HyperLip", "CAD", "Death",
                             "Triple Risk", "CAD", "Death",
                             "Triple Risk", "CAD", "Death",
                             "Triple Risk", "CAD", "Death",
                             "CAD", "Death", "Death",
                             "Death"),
                      arrows = 'to')

# Create the visualization
visNetwork(nodes, edges, height = "600px", width = "100%") %>%
  visLayout(randomSeed = 123, 
            hierarchical = list(enabled = TRUE, direction = "UD", 
                                levelSeparation = 200, nodeSpacing = 150)) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), nodesIdSelection = TRUE) %>%
  visLegend()
