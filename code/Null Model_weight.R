rm(list = ls())
library(igraph)
library(ggplot2)

# Read the data
data <- read_excel("第七次.xls", col_names = TRUE)

# Convert the data frame to a matrix
data <- as.data.frame(data)
rownames(data) <- data[[1]]  # Set row names using the first column
data <- data[, -1]           # Remove the first column
migration_matrix <- as.matrix(data)

# Transpose the matrix to switch rows and columns
migration_matrix <- t(migration_matrix)

log_migration_matrix <- log1p(migration_matrix)

total_sum <- sum(log_migration_matrix)
log_migration_matrix_norm <- log_migration_matrix / total_sum

my_matrix <- log_migration_matrix_norm

g_my <- graph_from_adjacency_matrix(my_matrix,weighted = TRUE, mode = "directed", diag = FALSE)
in_strength <- strength(g_my, mode = "in", weights = E(g_my)$weight)
out_strength <- strength(g_my, mode = "out", weights = E(g_my)$weight)
NS <- in_strength + out_strength
########################################################### MLE ###########################################################
# Objective function for node strength (in and out)
objective_function_strengths <- function(params, my_matrix, in_strength, out_strength) {
  # Extract parameters y_in and y_out for each node
  y_in <- params[1:length(in_strength)]  # y_i^in for each node
  y_out <- params[(length(in_strength) + 1):(2 * length(in_strength))]  # y_i^out for each node
  
  # Initialize the objective value
  objective_value <- 0
  
  # Loop through each node for in-strength and out-strength
  for (i in 1:length(in_strength)) {
    
    # Compute the observed in-strength and out-strength
    observed_in_strength <- sum(my_matrix[, i])  # Sum of weights coming to node i
    observed_out_strength <- sum(my_matrix[i, ])  # Sum of weights going out of node i
    
    # Compute the expected in-strength for node i
    expected_in_strength <- sum(sapply(1:length(in_strength), function(j) {
      if (i != j) {
        return(y_in[i] * y_out[j] / (1 - y_in[i] * y_out[j]))
      }
      return(0)
    }))
    
    # Compute the expected out-strength for node i
    expected_out_strength <- sum(sapply(1:length(in_strength), function(j) {
      if (i != j) {
        return(y_out[i] * y_in[j] / (1 - y_out[i] * y_in[j]))
      }
      return(0)
    }))
    
    # Calculate the squared differences between observed and expected strengths
    objective_value <- objective_value + (observed_in_strength - expected_in_strength)^2 + (observed_out_strength - expected_out_strength)^2
  }
  
  # Return the total objective value (we want to minimize this)
  return(objective_value)
}

# Optimizing the parameters using the adjusted objective function for the weighted case
optimize_parameters_strengths <- function(my_matrix, in_strength, out_strength) {
  # Initial guess for the parameters (e.g., uniform distribution based on degree)
  #initial_params <- c(in_strength / mean(in_strength), out_strength / mean(out_strength))
  #initial_params <- c(rep(1, length(in_strength)), rep(1, length(in_strength)))
  # Log-transform the initial guess
  initial_y_in <- log(1 + in_strength / mean(in_strength))  # Avoid log(0) by adding 1
  initial_y_out <- log(1 + out_strength / mean(out_strength))
  
  # Combine the parameters
  initial_params <- c(initial_y_in, initial_y_out)
  # Perform the optimization
  result <- optim(initial_params, fn = objective_function_strengths, my_matrix = my_matrix,
                  in_strength = in_strength, out_strength = out_strength, method = "BFGS",
                  control = list(maxit = 1000, trace = 2))
  
  # Return the optimized parameters
  return(result$par)
}

# Call the function to optimize
optimized_parameters_strengths <- optimize_parameters_strengths(my_matrix, in_strength, out_strength)

# Print the optimized parameters
print(optimized_parameters_strengths)

########################################################### TESTING ###########################################################
# Function to check degree errors
check_strength_error <- function(params, my_matrix, in_strength, out_strength) {
  y_in <- params[1:length(in_strength)]  # x_i^in for each node
  y_out <- params[(length(in_strength) + 1):length(params)]  # x_i^out for each node
  
  total_error_in <- 0
  total_error_out <- 0
  
  # Loop through each node
  for (i in 1:length(in_strength)) {
    # Calculate expected in-degree
    expected_in_strength <- sum(sapply(1:length(in_strength), function(j) {
      if (i != j) {
        return(y_in[i] *y_out[j] / (1 - y_in[i] * y_out[j]))
      }
      return(0)
    }))
    
    # Calculate expected out-degree
    expected_out_strength <- sum(sapply(1:length(out_strength), function(j) {
      if (i != j) {
        return(y_out[i] * y_in[j] / (1 - y_out[i] * y_in[j]))
      }
      return(0)
    }))
    
    # Accumulate error for in-degree and out-degree
    total_error_in <- total_error_in + abs(in_strength[i] - expected_in_strength)
    total_error_out <- total_error_out + abs(out_strength[i] - expected_out_strength)
  }
  
  # Return total error for in and out degrees
  return(c(total_error_in, total_error_out))
}

# Example of checking errors
errors <- check_strength_error(optimized_parameters_strengths, my_matrix, in_strength, out_strength)
print(errors)  # Print the error values



# Plotting degree distributions: Empirical vs Model
par(mfrow=c(1, 2))  # Two side-by-side plots

# Empirical in-degree distribution
hist(in_strength, main="Empirical In-Strength", xlab="In-Strength", col=rgb(0.2, 0.5, 0.2, 0.5))

# Model in-degree distribution (from optimized parameters)
expected_in_strength <- sapply(1:length(in_strength), function(i) {
  sum(sapply(1:length(in_strength), function(j) {
    if (i != j) {
      return(optimized_parameters_strengths[i] * optimized_parameters_strengths[j + length(in_strength)] / (1 - optimized_parameters_strengths[i] * optimized_parameters_strengths[j + length(in_strength)]))
    }
    return(0)
  }))
})

hist(expected_in_strength, main="Model In-Strength", xlab="In-Strength", col=rgb(0.5, 0.2, 0.2, 0.5))

# Similarly for out-degrees
# Empirical out-degree distribution
hist(out_strength, main="Empirical Out-Strength", xlab="Out-Strength", col=rgb(0.2, 0.5, 0.2, 0.5))

# Model out-degree distribution
expected_out_strength <- sapply(1:length(out_strength), function(i) {
  sum(sapply(1:length(out_strength), function(j) {
    if (i != j) {
      return(optimized_parameters_strengths[i + length(in_strength)] * optimized_parameters_strengths[j] / (1 + optimized_parameters_strengths[i + length(in_strength)] * optimized_parameters_strengths[j]))
    }
    return(0)
  }))
})

hist(expected_out_strength, main="Model Out-Strength", xlab="Out-Strength", col=rgb(0.5, 0.2, 0.2, 0.5))


########################################################### EXPECTED ###########################################################
# P_matrix
B <- my_matrix  
y_in <- optimized_parameters_strengths[1:31]
y_out <- optimized_parameters_strengths[32:62]

# Step: Total Strength  (the same as NS)
strength_tot <- rowSums(B) + colSums(B)

# Step: Calculate p_ij and p_ji based on model parameters (y_in, y_out)
p_matrix <- matrix(0, nrow = nrow(B), ncol = ncol(B))  # Initialize probability matrix

# Compute p_ij based on model parameters
for (i in 1:nrow(B)) {
  for (j in 1:ncol(B)) {
    if (i != j) {
      p_matrix[i, j] <- (y_out[i] * y_in[j]) / (1 - y_out[i] * y_in[j])  # p_ij
    }
  }
}
################################ ANNS ################################
knn_observe_wei <- knn(g_my, mode = "all", weights = E(g_my)$weight)
ANNS_observed <- knn_result_wei$knn

cor(NS,ANNS_observed, use = "complete.obs")

g_p <- graph_from_adjacency_matrix(p_matrix,weighted = TRUE, mode = "directed", diag = FALSE)
knn_expect_wei <- knn(g_p, mode = "all", weights = E(g_p)$weight)
ANNS_expected <- knn_result_wei$knn
cor(NS,ANNS_expected, use = "complete.obs")

################################## observed
k_tot <- rowSums(my_matrix) + colSums(my_matrix)
# Initialize a vector to store the observed ANNS for each node
anns_observed <- numeric(length(k_tot))

# Loop through each node to calculate its observed ANNS
for (i in 1:length(k_tot)) {
  # Get the neighbors of node i (i.e., non-zero elements in row i and column i of the matrix)
  neighbors <- which(my_matrix[i, ] > 0 | my_matrix[, i] > 0)
  
  # Sum (a_ij + a_ji) * s_j_tot for each neighbor
  sum_k_j_tot <- sum((my_matrix[i, neighbors] + my_matrix[neighbors, i]) * k_tot[neighbors])
  
  # Compute the observed ANND for node i
  anns_observed[i] <- sum_k_j_tot / k_tot[i]
}


# Print observed ANNS
print(anns_observed)


################################### expected

# Initialize a vector to store the observed ANNS for each node
anns_expected<- numeric(length(strength_tot))

# Loop through each node to calculate its observed ANNS
for (i in 1:length(strength_tot)) {
  # Get the neighbors of node i (i.e., non-zero elements in row i and column i of the matrix)
  neighbors <- which(p_matrix[i, ] > 0 | p_matrix[, i] > 0)
  
  # Sum (p_ij + p_ji) * s_j_tot for each neighbor
  sum_k_j_tot <- sum((p_matrix[i, neighbors] + p_matrix[neighbors, i]) * strength_tot[neighbors])
  
  # Compute the observed ANND for node i
  anns_expected[i] <- sum_k_j_tot / strength_tot[i]
}


###################################### comparison

# Compute the difference between observed and expected ANND
diff <- anns_observed - anns_expected

# Plot the difference
plot(diff, main = "Difference between Observed and Expected ANNS", ylab = "Difference", xlab = "Node")

# Calculate standard deviation for observed and expected ANND
observed_sd <- sd(anns_observed,na.rm = TRUE)
expected_sd <- sd(anns_expected,na.rm = TRUE)

# Calculate the 95% confidence intervals
N <- length(NS)  # Number of nodes

observed_CI_lower <- anns_observed - 1.96 * observed_sd / sqrt(N)
observed_CI_upper <- anns_observed + 1.96 * observed_sd / sqrt(N)

expected_CI_lower <- anns_expected - 1.96 * expected_sd / sqrt(N)
expected_CI_upper <- anns_expected + 1.96 * expected_sd / sqrt(N)


library(ggplot2)

# Create a data frame for ggplot
df <- data.frame(
  Node = 1:N,
  Observed_ANNS = anns_observed,
  Expected_ANNS = anns_expected,
  Observed_CI_lower = observed_CI_lower,
  Observed_CI_upper = observed_CI_upper,
  Expected_CI_lower = expected_CI_lower,
  Expected_CI_upper = expected_CI_upper
)

# Plot the graph with observed and expected ANNS with confidence bands
ggplot(df, aes(x = Node)) +
  # Plot observed ANND with confidence interval
  geom_line(aes(y = Observed_ANNS), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Observed_CI_lower, ymax = Observed_CI_upper), fill = "blue", alpha = 0.2) +
  
  # Plot expected ANND with confidence interval
  geom_line(aes(y = Expected_ANNS), color = "red", size = 1) +
  geom_ribbon(aes(ymin = Expected_CI_lower, ymax = Expected_CI_upper), fill = "red", alpha = 0.2) +
  
  # Customize labels and titles
  labs(title = "Observed and Expected ANNS with 95% Confidence Intervals",
       x = "Node", y = "ANNS") +
  theme_minimal()



cor_observed <- cor(NS, anns_observed, use = "complete.obs")
cor_expected <- cor(NS, anns_expected, use = "complete.obs")

# Function to calculate 95% confidence interval for Pearson correlation
ci_pearson <- function(cor_value, n) {
  z <- 0.5 * log((1 + cor_value) / (1 - cor_value))  # Fisher's Z-transformation
  se <- 1 / sqrt(n - 3)  # Standard error
  z_lower <- z - 1.96 * se
  z_upper <- z + 1.96 * se
  cor_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)  # Inverse Fisher Z
  cor_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)
  return(c(cor_lower, cor_upper))
}

# Number of data points
n <- length(NS)

# Calculate confidence intervals for both correlations
ci_observed <- ci_pearson(cor_observed, n)
ci_expected <- ci_pearson(cor_expected, n)

# Print the Pearson correlations and their confidence intervals
cat("Pearson Correlation (Observed ANNS vs NS):", cor_observed, "\n")
cat("95% Confidence Interval (Observed ANNS vs NS):", ci_observed, "\n")
cat("Pearson Correlation (Expected ANNS vs NS):", cor_expected, "\n")
cat("95% Confidence Interval (Expected ANNS vs NS):", ci_expected, "\n")

# Plotting with ggplot2
library(ggplot2)

# Create a data frame for plotting
df <- data.frame(NS, anns_observed, anns_expected)

# Create the plot
ggplot(df, aes(x = NS)) +
  geom_point(aes(y = anns_observed), color = "blue", alpha = 0.6) +
  geom_smooth(aes(y = anns_observed), method = "lm", color = "blue", se = TRUE) +  # Linear fit with confidence bands
  geom_point(aes(y = anns_expected), color = "red", alpha = 0.6) +
  geom_smooth(aes(y = anns_expected), method = "lm", color = "red", se = TRUE) +  # Linear fit with confidence bands
  labs(x = "Node Strength (NS)", y = "Average Nearest Neighbor strength (ANNS)", 
       title = "Observed vs Expected ANNS vs Node Strength") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "none")


########################################### WCC ######################################################################################

######################################### observed

### TRY 1
# Function to calculate weighted clustering coefficient (WCC) based on triangles
calculate_WCC <- function(my_matrix) {
  # Step 1: Cube-root the weights (adjust the matrix element-wise)
  M_cube_root <- my_matrix^(1/3)  # Element-wise cube root
  
  # Step 2: Compute total degree (k_i^{tot}) for each node
  k_tot <- rowSums(M_cube_root) + colSums(M_cube_root)
  
  # Step 3: Compute reciprocated edges (k_i^{↔}) for each node
  reciprocal_links <- (M_cube_root > 0) * (t(M_cube_root) > 0)
  reciprocal_count <- rowSums(reciprocal_links)
  
  # Step 4: Initialize the vector to store WCC values for each node
  WCC <- numeric(length(k_tot))
  
  # Step 5: Loop through each node to compute WCC
  for (i in 1:length(k_tot)) {
    # Get neighbors for node i (nodes j and k where w_ij > 0 or w_ji > 0)
    neighbors <- which(M_cube_root[i, ] > 0 | M_cube_root[, i] > 0)
    
    # Initialize the sum for the numerator of WCC
    numerator <- 0
    
    # Loop through pairs of neighbors (j, k) to calculate the weighted interactions
    for (j in neighbors) {
      for (k in neighbors) {
        if (i != j && i != k && j != k) {
          # Sum the weighted interactions as described in the formula
          numerator <- numerator + (M_cube_root[i, j] + M_cube_root[j, i]) * 
            (M_cube_root[j, k] + M_cube_root[k, j]) * 
            (M_cube_root[k, i] + M_cube_root[i, k])
          
          # Only print the numerator if it's negative
          if (numerator < 0) {
            print(numerator)   # This will print when numerator is negative
          }
        }
      }
    }
    
    # Step 6: Compute the denominator for WCC (based on node degree and reciprocated edges)
    denominator <- 2 * (k_tot[i] * (k_tot[i] - 1) - 2 * reciprocal_count[i])
    if (denominator < 0) {
      print(denominator)   # This will print when numerator is negative
    }
    # Step 7: Compute the WCC for node i (if denominator is non-zero)
    if (denominator != 0) {
      WCC[i] <- numerator / denominator
    } else {
      WCC[i] <- NA  # Assign NA if denominator is 0 (isolated or non-connected nodes)
    }
  }
  
  # Return the WCC values
  return(WCC)
}

#calculate
observed_WCC <- calculate_WCC(my_matrix)
print(observed_WCC)

observed_WCC <- calculate_WCC_with_weights(my_matrix)



#### TRY 2

# Function to calculate WCC with weight of reciprocated links
calculate_WCC_with_weights <- function(my_matrix) {
  # Step 1: Compute M^(1/3)
  M_cube_root <- my_matrix^(1/3)  # Element-wise cube root
  
  # Step 2: Compute total degree (sum of row or column)
  k_tot <- rowSums(my_matrix) + colSums(my_matrix)
  # Step 3: Calculate the weight of reciprocated links
  reciprocated_weights <- (M_cube_root > 0) * (t(M_cube_root) > 0)  # Matrix of reciprocated link indicators
  weight_reciprocal_links <- sum(M_cube_root * reciprocated_weights)  # Sum of reciprocated link weights
  
  # Initialize WCC vector
  WCC <- numeric(length(k_tot))
  
  # Step 4: Loop through each node to calculate WCC
  for (i in 1:length(k_tot)) {
    # Get the neighbors of node i (non-zero values in row i and column i)
    neighbors <- which(my_matrix[i, ] > 0 | my_matrix[, i] > 0)
    
    numerator <- 0
    for (j in neighbors) {
      for (k in neighbors) {
        if (i != j && i != k && j != k) {
          # Sum the weighted interactions as described in the formula
          numerator <- numerator + (M_cube_root[i, j] + M_cube_root[j, i]) * 
            (M_cube_root[j, k] + M_cube_root[k, j]) * 
            (M_cube_root[k, i] + M_cube_root[i, k])
        }
      }
    }
    
    # Denominator: k_i_tot * (k_i_tot - 1) - 2 * weight of reciprocated links for node i
    denominator <- k_tot[i] * (k_tot[i] - 1) - 2 * weight_reciprocal_links
    if (denominator < 0) {
      print(denominator)   # This will print when numerator is negative
    }
    if (denominator != 0) {
      WCC[i] <- numerator / (2 * denominator)
    } else {
      WCC[i] <- NA  # Handle case where denominator is 0 or negative
    }
  }
  
  return(WCC)
}


######### TRY 3

# Function to calculate WCC with weight of reciprocated links for each node

calculate_WCC_with_weights <- function(my_matrix) {
  # Step 1: Compute M^(1/3)
  M_cube_root <- my_matrix^(1/3)  # Element-wise cube root
  
  # Step 2: Compute total degree (sum of row or column)
  k_tot <- rowSums(my_matrix) + colSums(my_matrix)
  
  # Initialize WCC vector
  WCC <- numeric(length(k_tot))
  
  # Step 3: Loop through each node to calculate WCC
  for (i in 1:length(k_tot)) {
    # Get the neighbors of node i (non-zero values in row i and column i)
    neighbors <- which(my_matrix[i, ] > 0 | my_matrix[, i] > 0)
    
    # Step 4: Calculate weight of reciprocated links for node i
    reciprocated_weight_i <- 0
    for (j in neighbors) {
      if (my_matrix[i, j] > 0 && my_matrix[j, i] > 0) {
        # Sum the weight of the reciprocated link between i and j
        reciprocated_weight_i <- reciprocated_weight_i + M_cube_root[i, j] + M_cube_root[j, i]
      }
    }
    
    numerator <- 0
    # Step 5: Calculate numerator (sum of weighted interactions in triangles involving node i)
    for (j in neighbors) {
      for (k in neighbors) {
        if (i != j && i != k && j != k) {
          # Sum the weighted interactions as described in the formula
          numerator <- numerator + (M_cube_root[i, j] + M_cube_root[j, i]) * 
            (M_cube_root[j, k] + M_cube_root[k, j]) * 
            (M_cube_root[k, i] + M_cube_root[i, k])
        }
      }
    }
    
    # Denominator: k_i_tot * (k_i_tot - 1) - 2 * weight of reciprocated links for node i
    denominator <- k_tot[i] * (k_tot[i] - 1) - 2 * reciprocated_weight_i
    if (denominator < 0) {
      print(denominator)   # This will print when numerator is negative
    }
    if (denominator != 0) {
      WCC[i] <- numerator / (2 * denominator)
    } else {
      WCC[i] <- NA  # Handle case where denominator is 0 or negative
    }
  }
  
  return(WCC)
}

observed_WCC <- calculate_WCC_with_weights(my_matrix)

cor(NS,observed_WCC)
############ TRY 4
# Function to calculate WCC with a safer adjustment to the denominator
calculate_WCC_with_adjusted_denominator <- function(my_matrix) {
  # Step 1: Compute M^(1/3)
  M_cube_root <- my_matrix^(1/3)  # Element-wise cube root
  
  # Step 2: Compute total degree (sum of row or column)
  k_tot <- rowSums(my_matrix) + colSums(my_matrix)
  
  # Initialize WCC vector
  WCC <- numeric(length(k_tot))
  
  # Step 3: Loop through each node to calculate WCC
  for (i in 1:length(k_tot)) {
    # Get the neighbors of node i (non-zero values in row i and column i)
    neighbors <- which(my_matrix[i, ] > 0 | my_matrix[, i] > 0)
    
    # Step 4: Calculate weight of reciprocated links for node i
    reciprocated_weight_i <- 0
    for (j in neighbors) {
      if (my_matrix[i, j] > 0 && my_matrix[j, i] > 0) {
        # Sum the weight of the reciprocated link between i and j
        reciprocated_weight_i <- reciprocated_weight_i + M_cube_root[i, j] + M_cube_root[j, i]
      }
    }
    
    numerator <- 0
    # Step 5: Calculate numerator (sum of weighted interactions in triangles involving node i)
    for (j in neighbors) {
      for (k in neighbors) {
        if (i != j && i != k && j != k) {
          # Sum the weighted interactions as described in the formula
          numerator <- numerator + (M_cube_root[i, j] + M_cube_root[j, i]) * 
            (M_cube_root[j, k] + M_cube_root[k, j]) * 
            (M_cube_root[k, i] + M_cube_root[i, k])
        }
      }
    }
    
    # Step 6: Adjust the denominator to prevent negative values
    denominator <- k_tot[i] * (k_tot[i] - 1) - 2 * reciprocated_weight_i
    # Add a small regularization factor if the denominator is too small or negative
    if (denominator <= 0) {
      denominator <- denominator + 15  # how to deal with denominator < 0 ??
      print(denominator)
    }
    # Step 7: Compute WCC
    WCC[i] <- numerator / (2 * denominator)
  }
  
  return(WCC)
}

wcc_observed <- calculate_WCC_with_adjusted_denominator(my_matrix)

cor(NS,wcc_observed)

