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
# switch to binary
avg_link_weight <- mean(migration_matrix)
binary_matrix <- ifelse(migration_matrix > avg_link_weight  , 1, 0)

my_matrix <- binary_matrix

g_my <- graph_from_adjacency_matrix(my_matrix,weighted = TRUE, mode = "directed", diag = FALSE)
in_degree <- degree(g_my, mode = "in")
out_degree <- degree(g_my, mode = "out")
ND <- in_degree + out_degree

########################################################### MLE ###########################################################
# Objective function for node strength (in and out)
objective_function_degrees <- function(params, my_matrix, in_degree, out_degree) {
  # Extract parameters x_in and x_out for each node
  x_in <- params[1:length(in_degree)]  
  x_out <- params[(length(in_degree) + 1):(2 * length(in_degree))]  
  
  # Initialize the objective value
  objective_value <- 0
  
  # Loop through each node for in-strength and out-strength
  for (i in 1:length(in_degree)) {
    
    # Compute the observed in-strength and out-strength
    observed_in_degree <- sum(my_matrix[, i])  # Sum of weights coming to node i
    observed_out_degree <- sum(my_matrix[i, ])  # Sum of weights going out of node i
    
    # Compute the expected in-strength for node i
    expected_in_degree <- sum(sapply(1:length(in_degree), function(j) {
      if (i != j) {
        return(x_in[i] * x_out[j] / (1 + x_in[i] * x_out[j]))
      }
      return(0)
    }))
    
    # Compute the expected out-strength for node i
    expected_out_degree <- sum(sapply(1:length(in_degree), function(j) {
      if (i != j) {
        return(x_out[i] * x_in[j] / (1 + x_out[i] * x_in[j]))
      }
      return(0)
    }))
    
    # Calculate the squared differences between observed and expected strengths
    objective_value <- objective_value + (observed_in_degree - expected_in_degree)^2 + (observed_out_degree - expected_out_degree)^2
  }
  
  # Return the total objective value (we want to minimize this)
  return(objective_value)
}

# Optimizing the parameters using the adjusted objective function for the weighted case
optimize_parameters <- function(my_matrix, in_degree, out_degree) {
  # Initial guess for the parameters (e.g., uniform distribution based on degree)
  initial_params <- c(in_degree / mean(in_degree), out_degree / mean(out_degree))
  #initial_params <- c(rep(1, length(in_strength)), rep(1, length(in_strength)))

  # Perform the optimization
  result <- optim(initial_params, fn = objective_function_degrees, my_matrix = my_matrix,
                  in_degree = in_degree, out_degree = out_degree, method = "BFGS",
                  control = list(maxit = 1000, trace = 2))
  
  # Return the optimized parameters
  return(result$par)
}

# Call the function to optimize
optimized_parameters <- optimize_parameters(my_matrix, in_degree, out_degree)

# Print the optimized parameters
print(optimized_parameters)

###################################### Exaiming the result ######################################
################## method 1 degree error
# Function to check degree errors
check_degree_error <- function(params, my_matrix, in_degree, out_degree) {
  x_in <- params[1:length(in_degree)]  # x_i^in for each node
  x_out <- params[(length(in_degree) + 1):length(params)]  # x_i^out for each node
  
  total_error_in <- 0
  total_error_out <- 0
  
  # Loop through each node
  for (i in 1:length(in_degree)) {
    # Calculate expected in-degree
    expected_in_degree <- sum(sapply(1:length(in_degree), function(j) {
      if (i != j) {
        return(x_in[i] * x_out[j] / (1 + x_in[i] * x_out[j]))
      }
      return(0)
    }))
    
    # Calculate expected out-degree
    expected_out_degree <- sum(sapply(1:length(out_degree), function(j) {
      if (i != j) {
        return(x_out[i] * x_in[j] / (1 + x_out[i] * x_in[j]))
      }
      return(0)
    }))
    
    # Accumulate error for in-degree and out-degree
    total_error_in <- total_error_in + abs(in_degree[i] - expected_in_degree)
    total_error_out <- total_error_out + abs(out_degree[i] - expected_out_degree)
  }
  
  # Return total error for in and out degrees
  return(c(total_error_in, total_error_out))
}

# checking errors
errors <- check_degree_error(optimized_parameters, my_matrix, in_degree, out_degree)
print(errors)  # Print the error values

################## method 2 degree distribution

# Plotting degree distributions: Empirical vs Model
par(mfrow=c(1, 2))  # Two side-by-side plots

# Empirical in-degree distribution
hist(in_degree, main="Empirical In-Degree", xlab="In-Degree", col=rgb(0.2, 0.5, 0.2, 0.5))

# Model in-degree distribution (from optimized parameters)
expected_in_degrees <- sapply(1:length(in_degree), function(i) {
  sum(sapply(1:length(in_degree), function(j) {
    if (i != j) {
      return(optimized_parameters[i] * optimized_parameters[j + length(in_degree)] / (1 + optimized_parameters[i] * optimized_parameters[j + length(in_degree)]))
    }
    return(0)
  }))
})

hist(expected_in_degrees, main="Model In-Degree", xlab="In-Degree", col=rgb(0.5, 0.2, 0.2, 0.5))

# Similarly for out-degrees
# Empirical out-degree distribution
hist(out_degree, main="Empirical Out-Degree", xlab="Out-Degree", col=rgb(0.2, 0.5, 0.2, 0.5))

# Model out-degree distribution
expected_out_degrees <- sapply(1:length(out_degree), function(i) {
  sum(sapply(1:length(out_degree), function(j) {
    if (i != j) {
      return(optimized_parameters[i + length(in_degree)] * optimized_parameters[j] / (1 + optimized_parameters[i + length(in_degree)] * optimized_parameters[j]))
    }
    return(0)
  }))
})

hist(expected_out_degrees, main="Model Out-Degree", xlab="Out-Degree", col=rgb(0.5, 0.2, 0.2, 0.5))

################## method 3 Compare predicted vs observed edges

# Calculate expected edge probability matrix P_ij
P_matrix <- matrix(0, nrow = length(in_degree), ncol = length(out_degree))

for (i in 1:length(in_degree)) {
  for (j in 1:length(out_degree)) {
    P_matrix[i, j] <- optimized_parameters[i] * optimized_parameters[j + length(in_degree)] / (1 + optimized_parameters[i] * optimized_parameters[j + length(in_degree)])
  }
}

# Compare with actual adjacency matrix (real network)
# Note: a_ij is the observed adjacency matrix
diff <- abs(P_matrix - my_matrix)  # Compare predicted vs observed edges

# Sum of differences (lower value means better fit)
total_diff <- sum(diff)
print(total_diff)

###################################### Expected properties ############################################################################
x_in <- optimized_parameters[1:31]
x_out <- optimized_parameters[32:62]

A <- binary_matrix  # Binary adjacency matrix
#Total degree (sum of in-degree and out-degree)
degree_tot <- rowSums(A) + colSums(A)

#Calculate p_ij and p_ji based on model parameters (x_in, x_out)
p_matrix <- matrix(0, nrow = nrow(A), ncol = ncol(A))  # Initialize probability matrix

# Compute p_ij based on model parameters
for (i in 1:nrow(A)) {
  for (j in 1:ncol(A)) {
    if (i != j) {
      p_matrix[i, j] <- (x_out[i] * x_in[j]) / (1 + x_out[i] * x_in[j])  # p_ij
    }
  }
}


########################################### ANND ######################################################################################
k_tot <- rowSums(my_matrix) + colSums(my_matrix)

######################################### observed
# Initialize a vector to store the observed ANND for each node
k_nn_tot_observed <- numeric(length(k_tot))

# Loop through each node to calculate its observed ANND
for (i in 1:length(k_tot)) {
  # Get the neighbors of node i (i.e., non-zero elements in row i and column i of the matrix)
  neighbors <- which(my_matrix[i, ] == 1 | my_matrix[, i] == 1)
  
  # Sum (a_ij + a_ji) * k_j_tot for each neighbor
  sum_k_j_tot <- sum((my_matrix[i, neighbors] + my_matrix[neighbors, i]) * k_tot[neighbors])
  
  # Compute the observed ANND for node i
  k_nn_tot_observed[i] <- sum_k_j_tot / k_tot[i]
}

# Print observed ANND
print(k_nn_tot_observed)


###################################### expected
# Assuming x_in and x_out are vectors of the model parameters for in-degrees and out-degrees
# Initialize a vector to store the expected ANND for each node
k_nn_tot_expected <- numeric(length(k_tot))

# Loop through each node to calculate its observed ANND
for (i in 1:length(k_tot)) {
  # Get the neighbors of node i (i.e., non-zero elements in row i and column i of the matrix)
  neighbors <- which(p_matrix[i, ] > 0 | p_matrix[, i] > 0)
  
  # Sum (p_ij + p_ji) * s_j_tot for each neighbor
  sum_k_j_tot <- sum((p_matrix[i, neighbors] + p_matrix[neighbors, i]) * k_tot[neighbors])
  
  # Compute the expected ANND for node i, handling division by zero
  if (k_tot[i] != 0) {
    k_nn_tot_expected[i] <- sum_k_j_tot / k_tot[i]
  } else {
    k_nn_tot_expected[i] <- NA  
  }
}

###################################### comparison

# Compute the difference between observed and expected ANND
diff <- k_nn_tot_observed - k_nn_tot_expected

# Plot the difference
plot(diff, main = "Difference between Observed and Expected ANND", ylab = "Difference", xlab = "Node")


# Calculate standard deviation for observed and expected ANND
observed_sd <- sd(k_nn_tot_observed,na.rm = TRUE)
expected_sd <- sd(k_nn_tot_expected,na.rm = TRUE)

# Calculate the 95% confidence intervals
N <- length(k_tot)  # Number of nodes

observed_CI_lower <- k_nn_tot_observed - 1.96 * observed_sd / sqrt(N)
observed_CI_upper <- k_nn_tot_observed + 1.96 * observed_sd / sqrt(N)

expected_CI_lower <- k_nn_tot_expected - 1.96 * expected_sd / sqrt(N)
expected_CI_upper <- k_nn_tot_expected + 1.96 * expected_sd / sqrt(N)


library(ggplot2)

# Create a data frame for ggplot
df <- data.frame(
  Node = 1:N,
  Observed_ANND = k_nn_tot_observed,
  Expected_ANND = k_nn_tot_expected,
  Observed_CI_lower = observed_CI_lower,
  Observed_CI_upper = observed_CI_upper,
  Expected_CI_lower = expected_CI_lower,
  Expected_CI_upper = expected_CI_upper
)

# Plot the graph with observed and expected ANND with confidence bands
ggplot(df, aes(x = Node)) +
  # Plot observed ANND with confidence interval
  geom_line(aes(y = Observed_ANND), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Observed_CI_lower, ymax = Observed_CI_upper), fill = "blue", alpha = 0.2) +
  
  # Plot expected ANND with confidence interval
  geom_line(aes(y = Expected_ANND), color = "red", size = 1) +
  geom_ribbon(aes(ymin = Expected_CI_lower, ymax = Expected_CI_upper), fill = "red", alpha = 0.2) +
  
  # Customize labels and titles
  labs(title = "Observed and Expected ANND with 95% Confidence Intervals",
       x = "Node", y = "ANND") +
  theme_minimal()


Observed_ANND = k_nn_tot_observed
Expected_ANND = k_nn_tot_expected

cor_observed <- cor(ND, Observed_ANND, use = "complete.obs")
cor_expected <- cor(ND, Expected_ANND, use = "complete.obs")

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
n <- length(ND)

# Calculate confidence intervals for both correlations
ci_observed <- ci_pearson(cor_observed, n)
ci_expected <- ci_pearson(cor_expected, n)

# Print the Pearson correlations and their confidence intervals
cat("Pearson Correlation (Observed ANND vs ND):", cor_observed, "\n")
cat("95% Confidence Interval (Observed ANND vs ND):", ci_observed, "\n")
cat("Pearson Correlation (Expected ANND vs ND):", cor_expected, "\n")
cat("95% Confidence Interval (Expected ANND vs ND):", ci_expected, "\n")

# Plotting with ggplot2
library(ggplot2)

# Create a data frame for plotting
df <- data.frame(ND, Observed_ANND, Expected_ANND)

# Create the plot
ggplot(df, aes(x = ND)) +
  geom_point(aes(y = Observed_ANND), color = "blue", alpha = 0.6) +
  geom_smooth(aes(y = Observed_ANND), method = "lm", color = "blue", se = TRUE) +  # Linear fit with confidence bands
  geom_point(aes(y = Expected_ANND), color = "red", alpha = 0.6) +
  geom_smooth(aes(y = Expected_ANND), method = "lm", color = "red", se = TRUE) +  # Linear fit with confidence bands
  labs(x = "Node Degree (ND)", y = "Average Node Neighbor Degree (ANND)", 
       title = "Observed vs Expected ANND vs Node Degree") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "none")

########################################### BCC ######################################################################################

######################################### observed
# Function to calculate observed BCC
observed_BCC <- function(binary_matrix) {
  # Step 1: Adjacency matrix and transpose
  A <- binary_matrix  # Binary adjacency matrix
  A_T <- t(A)         # Transpose of the adjacency matrix
  
  # Step 2: Total degree (sum of in-degree and out-degree)
  degree_tot <- rowSums(A) + colSums(A)
  
  # Step 3: Reciprocated links
  reciprocal_links <- A * A_T
  reciprocal_count <- rowSums(reciprocal_links)
  
  # Step 4: Calculate (A + A_T)^3 (this gives the number of triangles)
  sum_matrix <- A + A_T
  triangle_matrix <- sum_matrix %*% sum_matrix %*% sum_matrix  # Cube of (A + A_T)
  triangles <- diag(triangle_matrix)  # Diagonal elements are the triangles for each node
  
  # Step 5: Calculate BCC for each node using the second formula
  BCC <- numeric(nrow(A))
  for (i in 1:nrow(A)) {
    if (degree_tot[i] > 1) {  # Avoid division by zero for nodes with degree <= 1
      potential_triangles <- 2 * (degree_tot[i] * (degree_tot[i] - 1) - 2 * reciprocal_count[i])
      if (potential_triangles > 0) {
        BCC[i] <- (triangles[i]) / potential_triangles  # Cube the triangle count for each node
      } else {
        BCC[i] <- NA  # If no potential triangles, set as NA
      }
    } else {
      BCC[i] <- NA  # Set as NA for nodes with degree 1 or less
    }
  }
  
  return(BCC)
}

######################################### expected 

# Function to calculate expected BCC using model parameters (x_in and x_out)
expected_BCC <- function(binary_matrix, x_in, x_out) {
  # Step 1: Adjacency matrix and transpose (not used for expected BCC)
  A <- binary_matrix  # Binary adjacency matrix (this is only for observed BCC)
  A_T <- t(A)         # Transpose of the adjacency matrix
  
  # Step 2: Total degree (sum of in-degree and out-degree)
  degree_tot <- rowSums(A) + colSums(A)
  
  # Step 3: Calculate p_ij and p_ji based on model parameters (x_in, x_out)
  p_matrix <- matrix(0, nrow = nrow(A), ncol = ncol(A))  # Initialize probability matrix
  
  # Compute p_ij based on model parameters
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if (i != j) {
        p_matrix[i, j] <- (x_out[i] * x_in[j]) / (1 + x_out[i] * x_in[j])  # p_ij
      }
    }
  }
  
  # Step 4: Calculate (p_ij + p_ji) * k_j_tot for each pair of neighbors (j, k)
  sum_matrix <- p_matrix + t(p_matrix)  # Sum of p_ij and p_ji for each node pair
  triangle_matrix <- sum_matrix %*% sum_matrix %*% sum_matrix  # Cube of (p_ij + p_ji)
  triangles <- diag(triangle_matrix)  # Diagonal elements are the expected triangles for each node
  
  # Step 5: Calculate expected BCC for each node using the expected formula
  expected_BCC_values <- numeric(nrow(A))
  for (i in 1:nrow(A)) {
    if (degree_tot[i] > 1) {  # Avoid division by zero for nodes with degree <= 1
      # Calculate the expected number of triangles based on p_ij, p_ji
      expected_triangles <- 0
      for (j in 1:nrow(A)) {
        for (k in 1:nrow(A)) {
          if (j != k) {
            p_ij <- p_matrix[i, j]  # p_ij
            p_ji <- p_matrix[j, i]  # p_ji
            expected_triangles <- expected_triangles + (p_ij + p_ji) * (sum(A[j, ]) + sum(A[k, ])) 
          }
        }
      }
      
      potential_triangles_expected <- 2 * (degree_tot[i] * (degree_tot[i] - 1) - 2 * sum(A[i, ] * A[, i]))
      if (potential_triangles_expected > 0) {
        expected_BCC_values[i] <- expected_triangles / potential_triangles_expected
      } else {
        expected_BCC_values[i] <- NA
      }
    } else {
      expected_BCC_values[i] <- NA
    }
  }
  
  return(expected_BCC_values)
}

###Apply
observed_BCC_values <- observed_BCC(binary_matrix)
expected_BCC_values <- expected_BCC(binary_matrix, x_in, x_out)


# Compute the difference between observed and expected BCC
diff_bcc <- observed_BCC_values - expected_BCC_values

# Plot the difference
plot(diff_bcc, main = "Difference between Observed and Expected BCC", ylab = "Difference", xlab = "Node")



# Combine the observed and expected BCC into a data frame
bcc_data <- data.frame(
  Observed_BCC = observed_BCC_values,
  Expected_BCC = expected_BCC_values
)

# Plot
ggplot(bcc_data, aes(x = Observed_BCC, y = Expected_BCC)) +
  geom_point(color = "blue", alpha = 0.6) +               # Scatter plot of observed vs predicted BCC
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of perfect correlation
  labs(title = "Observed vs Predicted BCC",
       x = "Observed BCC",
       y = "Predicted (Expected) BCC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        panel.grid.major = element_line(color = "gray90"),  # Lighter grid for better readability
        panel.grid.minor = element_line(color = "gray95"))  # Minor grid lines for clarity


#### visualiztion
# Calculate standard deviations for confidence intervals
observed_sd <- sd(observed_BCC_values, na.rm = TRUE)
predicted_sd <- sd(expected_BCC_values, na.rm = TRUE)

# Create a confidence interval (assuming 95% confidence)
ci_multiplier <- 1.96  # for 95% confidence
observed_upper <- observed_BCC_values + ci_multiplier * observed_sd
observed_lower <- observed_BCC_values - ci_multiplier * observed_sd

predicted_upper <- expected_BCC_values + ci_multiplier * predicted_sd
predicted_lower <- expected_BCC_values - ci_multiplier * predicted_sd


# Calculate the 95% confidence intervals
N <- length(k_tot)  # Number of nodes

library(ggplot2)

# Create a data frame for ggplot
df <- data.frame(
  Node = 1:N,
  Observed_BCC = observed_BCC_values,
  Expected_BCC = expected_BCC_values,
  Observed_CI_lower = observed_lower,
  Observed_CI_upper = observed_upper,
  Expected_CI_lower = predicted_lower,
  Expected_CI_upper = predicted_upper
)

# Plot the graph with observed and expected ANND with confidence bands
ggplot(df, aes(x = Node)) +
  # Plot observed ANND with confidence interval
  geom_line(aes(y = Observed_BCC), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = Observed_CI_lower, ymax = Observed_CI_upper), fill = "blue", alpha = 0.2) +
  
  # Plot expected ANND with confidence interval
  geom_line(aes(y = Expected_BCC), color = "red", size = 1) +
  geom_ribbon(aes(ymin = Expected_CI_lower, ymax = Expected_CI_upper), fill = "red", alpha = 0.2) +
  
  # Customize labels and titles
  labs(title = "Observed and Expected ANND with 95% Confidence Intervals",
       x = "Node", y = "ANND") +
  theme_minimal()

###Pearson Correlation

Observed_BCC = observed_BCC_values
Expected_BCC = expected_BCC_values

cor_observed <- cor(ND, Observed_BCC, use = "complete.obs")
cor_expected <- cor(ND, Expected_BCC, use = "complete.obs")

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
n <- length(ND)

# Calculate confidence intervals for both correlations
ci_observed <- ci_pearson(cor_observed, n)
ci_expected <- ci_pearson(cor_expected, n)

# Print the Pearson correlations and their confidence intervals
cat("Pearson Correlation (Observed BCC vs ND):", cor_observed, "\n")
cat("95% Confidence Interval (Observed BCC vs ND):", ci_observed, "\n")
cat("Pearson Correlation (Expected BCC vs ND):", cor_expected, "\n")
cat("95% Confidence Interval (Expected BCC vs ND):", ci_expected, "\n")

# Plotting with ggplot2
library(ggplot2)

# Create a data frame for plotting
df <- data.frame(ND, Observed_BCC, Expected_BCC)

# Create the plot
ggplot(df, aes(x = ND)) +
  geom_point(aes(y = Observed_BCC), color = "blue", alpha = 0.6) +
  geom_smooth(aes(y = Observed_BCC), method = "lm", color = "blue", se = TRUE) +  # Linear fit with confidence bands
  geom_point(aes(y = Expected_BCC), color = "red", alpha = 0.6) +
  geom_smooth(aes(y = Expected_BCC), method = "lm", color = "red", se = TRUE) +  # Linear fit with confidence bands
  labs(x = "Node Degree (ND)", y = "Binary Clustering Coefficient (BCC)", 
       title = "Observed vs Expected BCC vs Node Degree") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "none")

