rm(list = ls())

library(readxl)
library(readr)
library(openxlsx)
library(dplyr)
library(igraph)
library(poweRlaw)
library(showtext)
library(ggplot2)

setwd("~/Desktop/network")
############################################## MIGRATION MATRIX ########################################

migration5 <- read_excel("migration_matrix_5.xlsx")
migration5 <- migration5[, -1] 
rownames(migration5) <- colnames(migration5)  # Set row names using the first column
matrix5 <- as.matrix(migration5)


migration6 <- read_excel("migration_matrix_6.xlsx")
migration6 <- migration6[, -1]  
rownames(migration6) <- colnames(migration6)  # Set row names using the first column
matrix6 <- as.matrix(migration6)


migration7 <- read_excel("migration_matrix_7.xlsx")
migration7 <- migration7[, -1]  
rownames(migration7) <- colnames(migration7)  # Set row names using the first column
matrix7 <- as.matrix(migration7)

######################################################  Rescale ########################################################

# Create graphs from adjacency matrices
g_5 <- graph_from_adjacency_matrix(matrix5, mode = "directed", weighted = TRUE, diag = FALSE)
g_6 <- graph_from_adjacency_matrix(matrix6, mode = "directed", weighted = TRUE, diag = FALSE)
g_7 <- graph_from_adjacency_matrix(matrix7, mode = "directed", weighted = TRUE, diag = FALSE)

# Rescale by matrix volume
scale_by_volume <- function(matrix) {
  total_volume <- sum(matrix)
  matrix / total_volume
}

matrix5_rescaled <- scale_by_volume(matrix5)
matrix6_rescaled <- scale_by_volume(matrix6)
matrix7_rescaled <- scale_by_volume(matrix7)

######################################################  Calculate Node Strength and Link Weights ########################################################
# Calculate node strength (in + out) for each matrix
calc_node_strength <- function(graph) {
  strength(graph, mode = "all", weights = E(graph)$weight)
}

ns_5 <- calc_node_strength(graph_from_adjacency_matrix(matrix5_rescaled, mode = "directed", weighted = TRUE, diag = FALSE))
ns_6 <- calc_node_strength(graph_from_adjacency_matrix(matrix6_rescaled, mode = "directed", weighted = TRUE, diag = FALSE))
ns_7 <- calc_node_strength(graph_from_adjacency_matrix(matrix7_rescaled, mode = "directed", weighted = TRUE, diag = FALSE))

# Extract non-zero values for power-law fitting
ns_5_data <- ns_5[ns_5 > 0]
ns_6_data <- ns_6[ns_6 > 0]
ns_7_data <- ns_7[ns_7 > 0]

# Extract link weights
link_weights_5 <- matrix5_rescaled[matrix5_rescaled > 0]
link_weights_6 <- matrix6_rescaled[matrix6_rescaled > 0]
link_weights_7 <- matrix7_rescaled[matrix7_rescaled > 0]


# Combine node strengths (filter non-zero)
combined_ns <- c(ns_5[ns_5 > 0], ns_6[ns_6 > 0], ns_7[ns_7 > 0])

# Extract link weights and combine
link_weights_5 <- matrix5_rescaled[matrix5_rescaled > 0]
link_weights_6 <- matrix6_rescaled[matrix6_rescaled > 0]
link_weights_7 <- matrix7_rescaled[matrix7_rescaled > 0]
combined_lw <- c(link_weights_5, link_weights_6, link_weights_7)

######################################################  Fit Power-Law and Calculate Parameters ######################################################  
# Function to fit power law and get parameters
fit_power_law <- function(data) {
  pl <- conpl$new(data)  # Continuous power-law model
  est <- estimate_xmin(pl, xmax = max(data))  # Estimate xmin
  pl$setXmin(est)  # Set xmin
  list(xmin = est$xmin, alpha = est$pars, gof = est$gof)
}

# Calculate parameters for combined node strength and link weights
ns_fit <- fit_power_law(combined_ns)
lw_fit <- fit_power_law(combined_lw)

# Print Results
ns_fit
lw_fit

######################################################  plot  ######################################################  
# Function to calculate CCDF while preserving matrix information
calculate_ccdf_with_matrix <- function(data_list) {
  # Combine data while preserving matrix information
  combined_data <- data.frame(
    value = unlist(data_list),
    matrix = factor(rep(c("1995-2000", "2005-2010", "2015-2020"), 
                       times = sapply(data_list, length)))
  )
  
  # Sort by value in descending order
  combined_data <- combined_data[order(combined_data$value, decreasing = TRUE), ]
  
  # Calculate CCDF
  n <- nrow(combined_data)
  combined_data$ccdf <- (1:n) / n
  
  return(combined_data)
}

# Calculate CCDF for combined data while preserving matrix information
ccdf_ns <- calculate_ccdf_with_matrix(list(ns_5_data, ns_6_data, ns_7_data))
ccdf_lw <- calculate_ccdf_with_matrix(list(link_weights_5, link_weights_6, link_weights_7))

# Function to generate power law line
generate_power_law_line <- function(xmin, alpha, x_range, ccdf_data) {
  # Find the CCDF value at xmin
  ccdf_at_xmin <- ccdf_data$ccdf[which.min(abs(ccdf_data$value - xmin))]
  
  # For power law: P(X > x) ∝ x^(-α+1)
  # We normalize by the value at xmin to match the empirical CCDF
  y <- ccdf_at_xmin * (x_range/xmin)^(-alpha + 1)
  data.frame(x = x_range, y = y)
}

# Generate power law lines starting from xmin
x_range_ns <- seq(ns_fit$xmin, max(ccdf_ns$value), length.out = 100)
x_range_lw <- seq(lw_fit$xmin, max(ccdf_lw$value), length.out = 100)

power_law_ns <- generate_power_law_line(ns_fit$xmin, ns_fit$alpha, x_range_ns, ccdf_ns)
power_law_lw <- generate_power_law_line(lw_fit$xmin, lw_fit$alpha, x_range_lw, ccdf_lw)

# Create node strength distribution plot
showtext_auto()

p1 <- ggplot(ccdf_ns, aes(x = value, y = ccdf, color = matrix)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_line(data = power_law_ns, aes(x = x, y = y), 
            color = "black", linetype = "dashed", linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "人口迁移总量 (Node Strength)",
       y = "累积概率 (Cumulative Probability)",
       color = "") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.9))

# Create link weight distribution plot
p2 <- ggplot(ccdf_lw, aes(x = value, y = ccdf, color = matrix)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_line(data = power_law_lw, aes(x = x, y = y), 
            color = "black", linetype = "dashed", linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "省际人口迁移量 (Link Weight)",
       y = "累积概率 (Cumulative Probability)",
       color = "") +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.9))

# Save plots
ggsave("node_strength_ccdf.pdf", p1, width = 8, height = 6)
ggsave("link_weight_ccdf.pdf", p2, width = 8, height = 6)

# Print power law parameters
cat("\nNode Strength Power Law Parameters:\n")
cat("xmin:", ns_fit$xmin, "\n")
cat("alpha:", ns_fit$alpha, "\n")
cat("goodness of fit:", ns_fit$gof, "\n\n")

cat("Link Weight Power Law Parameters:\n")
cat("xmin:", lw_fit$xmin, "\n")
cat("alpha:", lw_fit$alpha, "\n")
cat("goodness of fit:", lw_fit$gof, "\n")
