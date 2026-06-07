rm(list = ls())
library(readxl)
library(openxlsx)
library(dplyr)
library(geosphere) # For distance calculation using latitude and longitude
library(poweRlaw) 
library(igraph)
library(sna)
library(intergraph)
library(ggplot2)

setwd("/Users/xiyu/Desktop/network")


################################################ descriptive data ###########################################

##weighted network


# Read the data
data <- read_excel("第七次.xls", col_names = TRUE)

# Convert the data frame to a matrix
data <- as.data.frame(data)
rownames(data) <- data[[1]]  # Set row names using the first column
data <- data[, -1]           # Remove the first column
migration_matrix <- as.matrix(data)

# Transpose the matrix to switch rows and columns
migration_matrix <- t(migration_matrix)

# describe the matrix
print(migration_matrix)
total_migrants <- sum(migration_matrix)
max_link_weight <- max(migration_matrix)
avg_link_weight <- mean(migration_matrix)
median_link_weight <- median(as.vector(migration_matrix))
min_link_weight <- min(migration_matrix[migration_matrix > 0], na.rm = TRUE)
num_links <- sum(migration_matrix > 0)

# average weighted index

##Compute the transpose of the matrix
A_T <- t(migration_matrix)

## Calculate Frobenius norms
frobenius_A <- sqrt(sum(migration_matrix^2))               # ||A||_F
frobenius_diff <- sqrt(sum((migration_matrix - A_T)^2))    # ||A - A^T||_F

##Compute the asymmetry index
S <- (frobenius_diff^2) / (2 * frobenius_A^2)

## Print the result
cat("Weighted Asymmetry Index (S):", S, "\n")

# Link_weight powerlaw
hist(migration_matrix[migration_matrix > 0], breaks = 50, main = "Link Weight Distribution", xlab = "Weight")
log_weights <- log(migration_matrix[migration_matrix > 0])
plot(log10(1:length(log_weights)), sort(log_weights, decreasing = TRUE), type = "l", 
     main = "Log-Log Plot of Link Weights", xlab = "Rank (Log Scale)", ylab = "Weight (Log Scale)")
pl_model <- displ$new(migration_matrix[migration_matrix > 0])
estimate_pars(pl_model)

link_weights <- as.vector(migration_matrix)
link_weights <- link_weights[link_weights > 0]
pl_model <- displ$new(link_weights)

pl_model$setXmin(1689) 
pl_model$setPars(estimate_pars(pl_model))

cat("Estimated power-law exponent:", pl_model$pars, "\n")
cat("Estimated xmin (cutoff):", pl_model$getXmin(), "\n")
plot(pl_model)
lines(pl_model, col = "red")


# Load necessary library
library(poweRlaw)

# Example data: replace with your migration data
link_weights <- as.vector(migration_matrix)  # Convert matrix to vector
link_weights <- link_weights[link_weights > 0]  # Filter positive values

# Initialize the power-law model
pl_model <- displ$new(link_weights)

# Estimate the power-law exponent and other parameters
pl_model$setXmin(1689)  # Set xmin (minimum value for the power law)
pl_model$setPars(estimate_pars(pl_model))  # Estimate the parameters

# Print the estimated power-law exponent and xmin (cutoff)
cat("Estimated power-law exponent:", pl_model$pars, "\n")
cat("Estimated xmin (cutoff):", pl_model$getXmin(), "\n")

# Plot the power-law distribution
# The plot function already includes the power-law fit, but we add more customizations
plot(pl_model, main = "Power-Law Distribution of Migration Link Weights",
     xlab = "Migration Link Weight (log scale)", ylab = "Probability Density (log scale)",
     col = "blue", pch = 16, cex = 0.7)

# Add the fitted power-law line to the plot
lines(pl_model, col = "red", lwd = 2)

# Add a legend for clarity
legend("bottomleft", legend = c("Observed Data", "Fitted Power-Law"), 
       col = c("blue", "red"), lty = c(NA, 1), pch = c(16, NA), lwd = c(NA, 2))

# Add a grid for better visualization
grid()




##binary network
binary_matrix <- ifelse(migration_matrix > avg_link_weight  , 1, 0)

n <- nrow(binary_matrix)
bi_links <- sum(binary_matrix > 0)
density <- bi_links / (n * (n - 1))

bilateral_links <- sum(binary_matrix * t(binary_matrix)) / 2
bilateral_density <- bilateral_links /bi_links

###SCC of binary network
g_bi <- graph_from_adjacency_matrix(binary_matrix, mode = "directed", diag = FALSE)
SCC <- clusters(g_bi, mode="strong")  
SCC$membership
SCC$csize
SCC$no


# V(g)$color <- rainbow(SCC$no)[SCC$membership]
# plot(g, mark.groups = split(1:vcount(g), SCC$membership))

###shortest path of binary network
g_bi_non <- graph_from_adjacency_matrix(binary_matrix, mode = "undirected", diag = FALSE)
graph_diameter <- diameter(g_bi_non, directed = FALSE)
avg_path_length <- mean_distance(g_bi_non, directed = FALSE)

####################################### Assortative ###########################################

##weighted network NS ANNS
g_wei <- graph_from_adjacency_matrix(migration_matrix, mode = "directed", weighted = TRUE, diag = FALSE)

in_strength <- strength(g_wei, mode = "in", weights = E(g_wei)$weight)
out_strength <- strength(g_wei, mode = "out", weights = E(g_wei)$weight)
NS <- strength(g_wei, mode = "all")

filtered_matrix <- migration_matrix
filtered_matrix[filtered_matrix < 500] <- 0
g_filtered <- graph_from_adjacency_matrix(filtered_matrix, mode = "directed", weighted = TRUE)
knn_result_wei <- knn(g_filtered, mode = "all", weights = E(g_filtered)$weight)
ANNS <- knn_result_wei$knn
average_ANNS <- mean(ANNS, na.rm = TRUE)
cor(NS, ANNS)
cor_coefficient <- cor(NS, ANNS, use = "complete.obs")
cat("Pearson Correlation Coefficient:", cor_coefficient, "\n")

# plot(NS, ANNS, 
#xlab = "Node Strength (NS)", 
#ylab = "Average Nearest Neighbor Strength (ANNS)", 
#main = "Scatterplot of NS vs. ANNS",
#pch = 19, col = "blue") 
# Add a regression line for better visualization
#abline(lm(ANNS ~ NS), col = "red", lwd = 2)

### log-transform way1
log_migration_matrix <- log1p(migration_matrix)
g_log <- graph_from_adjacency_matrix(log_migration_matrix, mode = "directed", weighted = TRUE, diag = FALSE)

in_strength_log <- strength(g_log, mode = "in", weights = E(g_log)$weight)
out_strength_log <- strength(g_log, mode = "out", weights = E(g_log)$weight)
NS_log <- strength(g_log, mode = "all")

filtered_matrix_log <- log_migration_matrix
filtered_matrix_log[filtered_matrix_log < 7.31] <- 0
g_filtered_log <- graph_from_adjacency_matrix(filtered_matrix_log, mode = "directed", weighted = TRUE)
knn_result_wei_log <- knn(g_filtered_log, mode = "all", weights = E(g_filtered_log)$weight)
ANNS_log <- knn_result_wei_log$knn
cor(NS, ANNS_log)
cor_coefficient_log <- cor(NS, ANNS_log, use = "complete.obs")
cat("Pearson Correlation Coefficient:", cor_coefficient_log, "\n")

### log-transform way2
log_migration_matrix <- log1p(migration_matrix)
filtered_matrix_log <- log_migration_matrix
avg_link_weight_log <- mean(log_migration_matrix)
filtered_matrix_log[filtered_matrix_log < 7.31] <- 0


g_filtered <- graph_from_adjacency_matrix(filtered_matrix_log, mode = "directed", weighted = TRUE, diag = FALSE)

in_strength_fil <- strength(g_filtered, mode = "in", weights = E(g_filtered)$weight)
out_strength_fil <- strength(g_filtered, mode = "out", weights = E(g_filtered)$weight)
NS_fil <- strength(g_filtered, mode = "all")

knn_result_wei_fil <- knn(g_filtered, mode = "all", weights = E(g_filtered)$weight)
ANNS_fil <- knn_result_wei_fil$knn
cor(NS_fil, ANNS_fil)
cor_coefficient_fil <- cor(NS_fil, ANNS_fil, use = "complete.obs")
cat("Pearson Correlation Coefficient:", cor_coefficient_fil, "\n")



##binary network ND ANND
g_bi <- graph_from_adjacency_matrix(binary_matrix, mode = "directed", weighted = TRUE, diag = FALSE)

in_degree <- degree(g_bi, mode = "in")
out_degree <- degree(g_bi, mode = "out")
ND <- degree(g_bi, mode = "all")
knn_result_bi <- knn(graph = g_bi, mode = "all")
ANND <- knn_result_bi$knn
average_ANND <- mean(ANND, na.rm = TRUE)
cor(ND, ANND, use = "complete.obs")


###sum ND and NS
node_metrics <- data.frame(
  Node = V(g_wei)$name,
  In_Degree = in_degree,
  Out_Degree = out_degree,
  In_Strength = in_strength,
  Out_Strength = out_strength
)

#write.xlsx(node_metrics, file = "/Users/xiyu/Desktop/社会网络/node_metrics.xlsx")

#write.xlsx(migration_matrix,file = "/Users/xiyu/Desktop/社会网络/migration_matrix_7.xlsx" )



####################################### Clustering ###########################################

#### binary node clustering coefficient BCC

# Step 1: Adjacency matrix and transpose
# Binary matrix
A <- binary_matrix  # Binary adjacency matrix
A_T <- t(A)         # Transpose of the adjacency matrix

# Total degree (sum of in-degree and out-degree)
degree_tot <- rowSums(A) + colSums(A)

# Reciprocated links
reciprocal_links <- A * A_T
reciprocal_count <- rowSums(reciprocal_links)

# Triangles (diagonal of (A + A_T)^3)
triangle_matrix <- (A + A_T) %*% (A + A_T) %*% (A + A_T)
triangles <- diag(triangle_matrix)

# Calculate BCC for each node
BCC <- numeric(nrow(A))
for (i in 1:nrow(A)) {
  if (degree_tot[i] > 1) {  # Avoid division by zero
    potential_triangles <- 2 * (degree_tot[i] * (degree_tot[i] - 1) - 2 * reciprocal_count[i])
    if (potential_triangles > 0) {
      BCC[i] <- triangles[i] / potential_triangles
    } else {
      BCC[i] <- NA
    }
  } else {
    BCC[i] <- NA
  }
}

# Display the BCC values
print(BCC)
average_BCC <- mean(BCC, na.rm = TRUE)

#ND NS BCC
cor_ND_BCC <- cor(ND, BCC, use = "complete.obs")
cat("Pearson Correlation between ND and BCC:", cor_ND_BCC, "\n")
cor_NS_BCC <- cor(NS, BCC, use = "complete.obs")
cat("Pearson Correlation between NS and BCC:", cor_NS_BCC, "\n")

#### weighted node clustering coefficient WCC
# Weighted matrix
migration_matrix_normalized <- migration_matrix / max(migration_matrix)
W <- migration_matrix_normalized # Weighted adjacency matrix
W_T <- t(W)            # Transpose of the adjacency matrix

# Normalize weights (cube root)
W_normalized <- W^(1/3)
W_T_normalized <- W_T^(1/3)

# Total degree (binary degree from weighted graph)
degree_tot <- rowSums(W > 0) + colSums(W > 0)

# Reciprocated links
reciprocal_links <- (W_normalized > 0) * (W_T_normalized > 0)
reciprocal_count <- rowSums(reciprocal_links)

# Weighted triangles
weighted_triangle_matrix <- (W_normalized + W_T_normalized) %*% 
  (W_normalized + W_T_normalized) %*% 
  (W_normalized + W_T_normalized)
weighted_triangles <- diag(weighted_triangle_matrix)

# Calculate WCC for each node
WCC <- numeric(nrow(W))
for (i in 1:nrow(W)) {
  if (degree_tot[i] > 1) {  # Avoid division by zero
    potential_triangles <- 2 * (degree_tot[i] * (degree_tot[i] - 1) - 2 * reciprocal_count[i])
    if (potential_triangles > 0) {
      WCC[i] <- weighted_triangles[i] / potential_triangles
    } else {
      WCC[i] <- NA
    }
  } else {
    WCC[i] <- NA
  }
}

# Display the WCC values
print(WCC)
mean(WCC)

#ND NS WCC
cor_ND_WCC <- cor(ND, WCC, use = "complete.obs")
cat("Pearson Correlation between ND and WCC:", cor_ND_WCC, "\n")
cor_NS_WCC <- cor(NS, WCC, use = "complete.obs")
cat("Pearson Correlation between NS and WCC:", cor_NS_WCC, "\n")


####################################### Community ###########################################

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
        expected_weight <- (NS_out[i] * NS_in[j]) / V
        Q <- Q + (A[i, j] - expected_weight)
      }
    }
  }
  return(Q / V )
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

# Create directed graph from migration matrix

result <- tabu_search(g_7, max_iterations = 100, tabu_list_size = 10, max_no_improve = 40)
cat("Best Modularity:", result$best_modularity, "\n")
print("Best Community Structure:")
print(result$best_community)

community_assignment <- data.frame(Node = V(g_7)$name, Community = result$best_community)
print(community_assignment)




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
    
    # Check if the number of unique communities is 4, stop if true
    if (length(unique(current_community)) == 4) {
      cat("Stopping early: Found 4 communities\n")
      break
    }
    
    iteration <- iteration + 1
  }
  
  return(list(best_community = best_community, best_modularity = best_modularity))
}

tabu_search(g_7, max_iterations = 50, tabu_list_size = 5, max_no_improve = 10)

################################## testing


cluster_edge_betweenness (g_bi, weights = NULL,
                          directed = TRUE,
                          edge.betweenness = TRUE)

####################################### Centrality ###########################################

#####Path Length

##Betweenness

betweenness_binary <- betweenness(g_bi, directed = TRUE)
betweenness_weighted <- betweenness(g_wei, directed = TRUE, weights = E(g_wei)$weight)
g_filtered
betweenness_weighted_filtered <- betweenness(g_filtered, directed = TRUE, weights = E(g_filtered)$weight)
##Closeness
closeness_binary_in <- closeness(g_bi, mode = "in", normalized = TRUE)
closeness_binary_out <- closeness(g_bi, mode = "out", normalized = TRUE)

closeness_weighted_in <- closeness(g_wei, mode = "in", weights = E(g_wei)$weight, normalized = TRUE)
closeness_weighted_out <- closeness(g_wei, mode = "out", weights = E(g_wei)$weight, normalized = TRUE)

##Harmonic
 
harmonic_centrality( g_bi, mode = "out")
harmonic_centrality( g_bi, mode = "in")
harmonic_centrality( g_bi, mode = "all")

harmonic_centrality( g_wei, mode = "out")
harmonic_centrality( g_wei, mode = "in")
harmonic_centrality( g_wei, mode = "all")

#####Eigenvector-based 

##Bonaich
bonaich_binary <- eigen_centrality(g_bi,directed = TRUE, weights = NULL)
bonaich_weighted <- eigen_centrality(g_wei,directed = TRUE, weights = E(g_wei)$weight)

alpha_centrality(g_bi, nodes = V(g_bi), alpha = 1, loops = FALSE, exo = 1, weights = NA, tol = 1e-07,sparse = TRUE)
alpha_centrality(g_wei, nodes = V(g_wei), alpha = 1, loops = FALSE, exo = 1, weights = E(g_wei)$weight, tol = 1e-07,sparse = TRUE)

##Page Rank
page_rank(g_bi, vids = V(g_bi), directed = TRUE, damping = 0.85)
page_rank(g_wei, vids = V(g_wei), directed = TRUE, damping = 0.85)

##Hits:authority & hub
hits_scores(g_bi,scale = TRUE)
hits_scores(g_wei,scale = TRUE,weights = NULL)

####################################### Null Model###########################################

###binary
binary_matrix
g_bi
in_degree <- degree(g_bi, mode = "in")
out_degree <- degree(g_bi, mode = "out")



####################################### Gravity Model ###########################################