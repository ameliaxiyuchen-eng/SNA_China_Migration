rm(list = ls())
library(readxl)
library(igraph)

setwd("~/Desktop/network")

#Step 1:  Function to calculate modularity for directed, weighted networks
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
tabu_search <- function(graph, max_iterations = 50, tabu_list_size = 5, max_no_improve = 10) {
  n <- vcount(graph)  # Number of nodes
  nodes <- V(graph)
  
  # Initial random assignment of nodes to communities
  current_community <- sample(1:n, n, replace = TRUE)
  best_community <- current_community
  best_modularity <- calculate_modularity(graph, best_community)
  
  # Initialize Tabu List
  tabu_list <- list()
  
  # Initialize counters
  no_improve <- 0
  iteration <- 1
  
  while (iteration <= max_iterations && no_improve <= max_no_improve) {
    cat("Iteration:", iteration, " | Best Modularity:", best_modularity, "\n")
    neighborhood <- list()
    modularity_scores <- c()
    
    # Generate neighbors by moving nodes to different communities
    for (node in 1:n) {
      current_assignment <- current_community
      for (new_community in setdiff(1:n, current_assignment[node])) {
        current_assignment[node] <- new_community
        neighborhood <- append(neighborhood, list(current_assignment))
        modularity_scores <- c(modularity_scores, calculate_modularity(graph, current_assignment))
      }
    }
    
    # Select the best neighbor not in the Tabu List
    best_neighbor_idx <- which.max(modularity_scores)
    best_neighbor <- neighborhood[[best_neighbor_idx]]
    best_neighbor_modularity <- modularity_scores[best_neighbor_idx]
    
    if (!(list(best_neighbor) %in% tabu_list)) {
      current_community <- best_neighbor
      current_modularity <- best_neighbor_modularity
      
      # Update Tabu List
      tabu_list <- append(tabu_list, list(best_neighbor))
      if (length(tabu_list) > tabu_list_size) {
        tabu_list <- tabu_list[-1]  # Remove oldest element
      }
      
      # Update best solution if improved
      if (current_modularity > best_modularity) {
        best_community <- current_community
        best_modularity <- current_modularity
        no_improve <- 0  # Reset no improvement counter
      } else {
        no_improve <- no_improve + 1
      }
    }
    
    iteration <- iteration + 1
  }
  
  return(list(best_community = best_community, best_modularity = best_modularity))
}

# Step 3: Prepare Migration Matrix and Graph


migration7 <- read_excel("migration_matrix_7.xlsx")


migration7 <- migration7[, -1]  
rownames(migration7) <- colnames(migration7)
migration_matrix <- as.matrix(migration7)

# Create directed graph from migration matrix
g <- graph_from_adjacency_matrix(migration_matrix, mode = "directed", weighted = TRUE)
result <- tabu_search(g, max_iterations = 100, tabu_list_size = 10, max_no_improve = 40)
cat("Best Modularity:", result$best_modularity, "\n")
print("Best Community Structure:")
print(result$best_community)

community_assignment <- data.frame(Node = V(g)$name, Community = result$best_community)
print(community_assignment)


