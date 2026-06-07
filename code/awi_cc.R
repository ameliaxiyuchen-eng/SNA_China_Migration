rm(list = ls())

library(readxl)
library(openxlsx)
library(dplyr)
library(poweRlaw) 
library(igraph)
library(sna)
library(intergraph)
library(ggplot2)
library(showtext)
library(wdnet)

setwd("/Users/xiyu/Desktop/network")


migration7 <- read_excel("migration_matrix_7.xlsx")
migration7 <- migration7[, -1]  
rownames(migration7) <- colnames(migration7)  # Set row names using the first column
matrix7 <- as.matrix(migration7)

g_7 <- graph_from_adjacency_matrix(matrix7, mode = "directed", weighted = TRUE, diag = FALSE)



ins_7 <- strength(g_7, mode = "in", weights = E(g_7)$weight)
out_7 <- strength(g_7, mode = "out", weights = E(g_7)$weight)
ns_7 <- strength(g_7, mode = "all")

########################################## Observed WAI ##############################
# Calculate (A - A^T)
diff_matrix_7 <- matrix7 - t(matrix7)

# Frobenius norm of (A - A^T)
numerator_7 <- sum(diff_matrix^2)

# Frobenius norm of A
denominator_7 <- 2 * sum(matrix7^2)

# Weighted Asymmetry Index
wai_7 <- numerator_7 / denominator_7

####### 6
# Calculate (A - A^T)
diff_matrix_6 <- matrix6 - t(matrix6)

# Frobenius norm of (A - A^T)
numerator_6 <- sum(diff_matrix^2)

# Frobenius norm of A
denominator_6 <- 2 * sum(matrix6^2)

# Weighted Asymmetry Index
wai_6 <- numerator_6 / denominator_6


####### 5
# Calculate (A - A^T)
diff_matrix_5 <- matrix5 - t(matrix5)

# Frobenius norm of (A - A^T)
numerator_5 <- sum(diff_matrix^2)

# Frobenius norm of A
denominator_5 <- 2 * sum(matrix5^2)

# Weighted Asymmetry Index
wai_5 <- numerator / denominator

########################################## Observed CC ##############################
netwk_7 <- igraph_to_wdnet(g_7)
netwk_6 <- igraph_to_wdnet(g_6)
netwk_5 <- igraph_to_wdnet(g_5)

cc_7 <- clustcoef(netwk_7, method = "Fagiolo")
cc_6 <- clustcoef(netwk_6, method = "Fagiolo")
cc_5 <- clustcoef(netwk_5, method = "Fagiolo")

cc <- clustcoef(netwk_7, method = "Fagiolo")
locallcc <- cc$total$localcc
# Save local clustering coefficients to Excel
write.xlsx(data.frame(LocalCC = locallcc), "local_clustering_coefficients.xlsx", rowNames = TRUE)


############################## NUll MODEL ##############################
library(igraph)

# Step 1: Create the original graph
matrix7 <- as.matrix(migration7)  # Replace with your adjacency matrix
g_7 <- graph_from_adjacency_matrix(matrix7, mode = "directed", weighted = TRUE, diag = FALSE)

# Step 2: Get in-strength and out-strength
ins_7 <- strength(g_7, mode = "in", weights = E(g_7)$weight)
out_7 <- strength(g_7, mode = "out", weights = E(g_7)$weight)

# Step 3: Asymmetry calculation function & Clustering Coefficient function
calculate_asymmetry <- function(mat) {
  diff_matrix <- mat - t(mat)
  numerator <- sum(diff_matrix^2)
  denominator <- 2 * sum(mat^2)
  wai <- numerator / denominator
  return(wai)
}

calculate_clustering <- function(mat) {
  g <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE, diag = FALSE)
  netwk <- igraph_to_wdnet(g)
  cc <- clustcoef(netwk, method = "Fagiolo")
  globalcc <- cc$total$globalcc
  return(globalcc)
}

# Step 4: Generate Fixed Strength Sequence Model (FSSM)
generate_fssm <- function(ins_7, out_7) {
  N <- length(ins_7)
  new_matrix <- matrix(0, N, N)
  
  for (i in 1:N) {
    out_weights <- runif(N, 0, 1)
    out_weights <- out_weights / sum(out_weights) * out_7[i]
    new_matrix[i, ] <- out_weights
  }
  
  for (j in 1:N) {
    correction_factor <- ins_7[j] / sum(new_matrix[, j])
    new_matrix[, j] <- new_matrix[, j] * correction_factor
  }
  
  return(new_matrix)
}

# Step 5: Calculate observed asymmetry & clustering coefficient
observed_asymmetry <- calculate_asymmetry(matrix7)
cat("Observed Asymmetry:", observed_asymmetry, "\n")

observed_globalcc <- calculate_clustering(matrix7)
cat("Observed Clustering Coefficient:", observed_globalcc, "\n")

# Step 6: Generate 1000 random graphs and calculate indices once
set.seed(42)  # For reproducibility
random_matrices <- replicate(1000, generate_fssm(ins_7, out_7), simplify = FALSE)

random_asymmetries <- sapply(random_matrices, calculate_asymmetry)
random_clustering <- sapply(random_matrices, calculate_clustering)

# Step 7: Compare observed vs random
cat("Mean Random Asymmetry:", mean(random_asymmetries), "\n")
cat("SD of Random Asymmetry:", sd(random_asymmetries), "\n")

cat("Mean Random Cluster:", mean(random_clustering), "\n")
cat("SD of Random Cluster:", sd(random_clustering), "\n")

# Visualization 

showtext_auto()
# 绘图
hist(random_asymmetries, breaks = 50, col = "lightblue", main = "",
     xlab = "不对称性指数", ylab = "频数", xlim = c(predicted_asymmetry - 0.05, max(random_asymmetries)))
abline(v = observed_asymmetry, col = "red", lwd = 2, lty = 2)
abline(v = predicted_asymmetry, col = "blue", lwd = 2, lty = 2)
legend("topright", 
       legend = c("观察值", "预测值"), 
       col = c("red", "blue"), 
       lty = 2, 
       lwd = 2)

hist(random_clustering, breaks = 50, col = "lightblue", main = "",
     xlab = "聚类系数", ylab = "频数")
abline(v = observed_globalcc, col = "red", lwd = 2, lty = 2)
abline(v = predicted_globalcc, col = "blue", lwd = 2, lty = 2)
legend("topright", 
       legend = c("观察值", "预测值"), 
       col = c("red", "blue"), 
       lty = 2, 
       lwd = 2)

# Step 8: Calculate p-value
p_value_asymmetry <- mean(random_asymmetries >= observed_asymmetry)
cat("Asymmetry Index P-value:", p_value_asymmetry, "\n")

p_value_clustering <- mean(random_clustering >= observed_globalcc)
cat("Clustering Coefficient P-value:", p_value_clustering, "\n")
