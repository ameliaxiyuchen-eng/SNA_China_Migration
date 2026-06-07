rm(list = ls())
library(readxl)
library(igraph)

setwd("/Users/xiyu/Desktop/社会网络")
data <- read_excel("第七次.xls", col_names = TRUE)


##weighted network

data <- as.data.frame(data)
# Set row names and remove the first column
rownames(data) <- data[[1]]
data <- data[, -1]
# Convert to numeric matrix
migration_matrix <- as.matrix(data)

# Read the data
data <- read_excel("第七次.xls", col_names = TRUE)

# Convert the data frame to a matrix
data <- as.data.frame(data)
rownames(data) <- data[[1]]  # Set row names using the first column
data <- data[, -1]           # Remove the first column
migration_matrix <- as.matrix(data)

# Function to calculate modularity for directed, weighted networks
calculate_modularity <- function(graph, community) {
  A <- as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)
  V <- sum(A)  # Total edge weight (volume)
  NS_out <- rowSums(A)  # Out-strength of nodes
  NS_in <- colSums(A)   # In-strength of nodes
  
  Q <- 0  # Initialize modularity
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if (community[i] == community[j]) {
        expected_weight <- (NS_out[i] * NS_in[j]) / (2 * V)
        Q <- Q + (A[i, j] - expected_weight)
      }
    }
  }
  return(Q / V)
}
# Step 2: Tabu Search Function
tabu_search <- function(graph, max_iterations = 20, tabu_list_size = 5, max_no_improve = 10) {
  n <- vcount(graph)
  nodes <- V(graph)
  current_community <- 1:n
  best_community <- current_community
  best_modularity <- calculate_modularity(graph, current_community)
  tabu_list <- list()
  no_improve <- 0
  iteration <- 1
  
  while (iteration <= max_iterations && no_improve <= max_no_improve) {
    cat("Iteration:", iteration, " | Best Modularity:", best_modularity, "\n")
    
    # Step 1: Calculate edge betweenness
    betweenness_weights <- edge_betweenness(graph, directed = TRUE, weights = E(graph)$weight)
    
    # Step 2: Find the highest-betweenness edge and remove it
    high_betweenness_edges <- order(betweenness_weights, decreasing = TRUE)
    graph <- delete_edges(graph, E(graph)[high_betweenness_edges[1]])
    
    # Recalculate betweenness for the new graph structure
    # (this happens automatically in the next iteration due to recalculating `edge_betweenness`)
    
    # Step 3: Recalculate community structure
    community_candidate <- cluster_edge_betweenness(graph, directed = TRUE, weights = E(graph)$weight)
    candidate_membership <- membership(community_candidate)
    candidate_modularity <- calculate_modularity(graph, candidate_membership)
    
    # Step 4: Update the best solution if modularity improves
    if (candidate_modularity > best_modularity && !(list(candidate_membership) %in% tabu_list)) {
      best_modularity <- candidate_modularity
      best_community <- candidate_membership
      no_improve <- 0
      
      # Update Tabu List
      tabu_list <- append(tabu_list, list(candidate_membership))
      if (length(tabu_list) > tabu_list_size) {
        tabu_list <- tabu_list[-1]
      }
    } else {
      no_improve <- no_improve + 1
    }
    
    iteration <- iteration + 1
  }
  
  return(list(best_community = best_community, best_modularity = best_modularity))
}

# Application
g <- graph_from_adjacency_matrix(migration_matrix, mode = "directed", weighted = TRUE)


# Step 4: Run Tabu Search
result <- tabu_search(g, max_iterations = 20, tabu_list_size = 2, max_no_improve = 10)

# Step 5: Display Results
cat("Best Modularity Achieved:", result$best_modularity, "\n")
cat("Best Community Structure:\n")
print(result$best_community)

# Visualize the graph with communities
plot(g, vertex.color = result$best_community, vertex.size = 10, vertex.label = NA, edge.arrow.size = 0.5)


