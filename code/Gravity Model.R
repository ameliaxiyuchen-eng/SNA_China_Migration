# Clear workspace and load required libraries
rm(list = ls())

# Install required packages if not already installed
if (!require("car")) {
  # Install dependencies first
  install.packages(c("carData", "MatrixModels", "quantreg", "lme4", "pbkrtest", "Matrix", "MatrixModels", "SparseM"))
  # Then install car
  install.packages("car", dependencies = TRUE)
}

if (!require("reshape2")) install.packages("reshape2")

# Load required libraries
library(readxl)
library(readr)
library(openxlsx)
library(dplyr)
library(geosphere)  # For distance calculation using latitude and longitude
library(igraph)
library(wdnet)
library(ggplot2)
library(car)  # For VIF calculation
library(reshape2)  # For correlation heatmap


# Set working directory
setwd("~/Desktop/network")

#=============================================================================
# 1. Data Loading and Initial Setup
#=============================================================================

# Read migration matrix data
migration7 <- read_excel("migration_matrix_7.xlsx")
migration7 <- migration7[, -1]  
rownames(migration7) <- colnames(migration7)
migration_matrix <- as.matrix(migration7)

# Read nodes data and perform log transformations
df <- read_csv("nodes_sum.csv")
df <- df %>%
  mutate(
    ln_pop = log(pop_2020),
    ln_pgdp = log(pgdp_2020)
  )

#=============================================================================
# 2. Matrix Calculations
#=============================================================================

# Calculate distance matrix using Vincenty formula
d_matrix <- matrix(0, nrow = nrow(df), ncol = nrow(df))
for (i in 1:nrow(df)) {
  for (j in 1:nrow(df)) {
    if (i != j) {
      d_matrix[i, j] <- distVincentySphere(
        c(df$lng[i], df$lat[i]), 
        c(df$lng[j], df$lat[j])
      )
    }
  }
}

# Calculate relative per capita GDP matrix
rY_matrix <- outer(df$pgdp_2020, df$pgdp_2020, FUN = "/")
rY_matrix <- 1 / rY_matrix

# Calculate common language and climate matrices
lang_matrix <- outer(df$lang, df$lang, FUN = function(x, y) as.integer(x == y))
climate_matrix <- outer(df$climate, df$climate, FUN = function(x, y) as.integer(x == y))

#=============================================================================
# 3. Data Preparation for OLS
#=============================================================================

# Create migration data frame
n <- nrow(migration_matrix)
migration_data <- expand.grid(source_id = 1:n, destination_id = 1:n)
migration_data$migration_flow <- as.vector(migration_matrix)
migration_data <- migration_data[migration_data$source_id != migration_data$destination_id, ]

# Prepare vectors for OLS
pop_2020 <- df$pop_2020
gdp_2020 <- df$pgdp_2020

pop_dest <- rep(pop_2020, each = length(pop_2020))
pop_source <- rep(pop_2020, times = length(pop_2020))
gdp_source <- rep(gdp_2020, times = length(gdp_2020))
gdp_dest <- rep(gdp_2020, each = length(gdp_2020))

# Create OLS data frame
ols_2020 <- data.frame(
  migration = as.vector(migration_matrix),
  dist = as.vector(d_matrix),
  pop_source = pop_source,
  pop_dest = pop_dest,
  gdp_source = gdp_source,
  gdp_dest = gdp_dest,
  rY = as.vector(rY_matrix),
  comlang = as.vector(lang_matrix),
  comclimate = as.vector(climate_matrix)
)

# Filter and transform data
ols_2020_filtered <- ols_2020[ols_2020$migration > 0 & ols_2020$dist > 0, ]

# Log transformations for model variables
ln_migration <- log(ols_2020_filtered$migration)
ln_dist <- log(ols_2020_filtered$dist)
ln_pop_source <- log(ols_2020_filtered$pop_source)
ln_pop_dest <- log(ols_2020_filtered$pop_dest)
ln_gdp_source <- log(ols_2020_filtered$gdp_source)
ln_gdp_dest <- log(ols_2020_filtered$gdp_dest)
ln_rY <- log(ols_2020_filtered$rY)
COMLANG <- ols_2020_filtered$comlang
COMCLIMATE <- ols_2020_filtered$comclimate
  
#=============================================================================
# 4. Model Fitting
#=============================================================================

# Fit different model specifications
model_1 <- lm(ln_migration ~ ln_pop_source + ln_pop_dest + ln_dist + ln_rY)
model_2 <- lm(ln_migration ~ ln_pop_source + ln_pop_dest + ln_dist + ln_rY + COMLANG + COMCLIMATE) 
model_3 <- lm(ln_migration ~ ln_pop_source + ln_pop_dest + ln_dist + ln_rY + 
              COMLANG + COMCLIMATE + 
              factor(migration_data$source_id) + factor(migration_data$destination_id))

# Print model summaries
print("Model 1 Summary:")
print(summary(model_1))
print("\nModel 2 Summary:")
print(summary(model_2))
print("\nModel 3 Summary:")
print(summary(model_3))

#=============================================================================
# 4.1 Multicollinearity Analysis
#=============================================================================

# Create correlation matrix for independent variables
independent_vars <- data.frame(
  ln_pop_source = ln_pop_source,
  ln_pop_dest = ln_pop_dest,
  ln_dist = ln_dist,
  ln_rY = ln_rY,
  COMLANG = COMLANG,
  COMCLIMATE = COMCLIMATE
)

# Calculate correlation matrix
cor_matrix <- cor(independent_vars)

# Print correlation matrix
cat("\nCorrelation Matrix:\n")
print(cor_matrix)

# Calculate VIF manually
# Function to calculate VIF
calculate_vif <- function(model) {
  # Get the design matrix
  X <- model.matrix(model)[,-1]  # Remove intercept
  # Calculate VIF for each variable
  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)
  
  for(i in 1:ncol(X)) {
    # Fit model with current variable as response
    vif_model <- lm(X[,i] ~ X[,-i])
    # Calculate VIF
    vif_values[i] <- 1/(1 - summary(vif_model)$r.squared)
  }
  return(vif_values)
}

# Calculate VIF for model_2
vif_results <- calculate_vif(model_2)

# Print VIF results
cat("\nVariance Inflation Factors (VIF):\n")
print(vif_results)

# Create correlation heatmap
library(reshape2)

# Convert correlation matrix to long format for ggplot
cor_matrix_long <- melt(cor_matrix)

# Create heatmap
ggplot(cor_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                      midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print interpretation of VIF values
cat("\nVIF Interpretation:\n")
cat("VIF > 5 indicates moderate multicollinearity\n")
cat("VIF > 10 indicates severe multicollinearity\n")
cat("VIF > 1 indicates some correlation with other predictors\n")

#=============================================================================
# 5. Generate New Migration Matrix
#=============================================================================

# Extract coefficients from model_2
coefs <- coef(model_2)
intercept <- coefs[1]
alpha_pop_source <- coefs[2]
alpha_pop_dest <- coefs[3]
alpha_dist <- coefs[4]
alpha_rY <- coefs[5]
alpha_lang <- coefs[6]
alpha_climate <- coefs[7]

# Create new migration matrix
migration_matrix_new <- matrix(NA, nrow = length(df$pop_2020), ncol = length(df$pop_2020))
for (i in 1:length(df$pop_2020)) {
  for (j in 1:length(df$pop_2020)) {
    if (i != j) {
      migration_matrix_new[i, j] <- exp(
        intercept + 
        alpha_pop_source * log(df$pop_2020[i]) +
        alpha_pop_dest * log(df$pop_2020[j]) +
        alpha_dist * log(d_matrix[i, j]) +
        alpha_rY * log(rY_matrix[i, j]) +
        alpha_lang * lang_matrix[i, j] +
        alpha_climate * climate_matrix[i, j]
      )
    } else {
      migration_matrix_new[i, j] <- 0
    }
  }
}

# Add row and column names from the original migration matrix
rownames(migration_matrix_new) <- rownames(migration_matrix)
colnames(migration_matrix_new) <- colnames(migration_matrix)

# Save the new migration matrix to Excel
write.xlsx(as.data.frame(migration_matrix_new), "predicted_migration_matrix.xlsx", rowNames = TRUE)


#=============================================================================
# 6. Network Analysis and Comparison
#=============================================================================

# Create graphs
g_wei <- graph_from_adjacency_matrix(migration_matrix, mode = "directed", weighted = TRUE, diag = FALSE)
g_new <- graph_from_adjacency_matrix(migration_matrix_new, mode = "directed", weighted = TRUE, diag = FALSE)

# Plot the network



# Calculate network statistics
in_strength <- strength(g_wei, mode = "in")
out_strength <- strength(g_wei, mode = "out")
NS <- in_strength + out_strength

in_strength_new <- strength(g_new, mode = "in")
out_strength_new <- strength(g_new, mode = "out")
NS_new <- in_strength_new + out_strength_new

# Plot link weight
# Create density plot for both matrices
showtext_auto()
plot(density(log10(migration_matrix[migration_matrix > 0])), 
     main = "",
     xlab = "log10(人口迁移量)", 
     ylab = "频率",
     col = "blue", 
     lwd = 2,
     xlim = c(0.5, max(log10(migration_matrix)) + 0.5),
     ylim = c(0, 0.8)) 

# Add density line for predicted values
lines(density(log10(migration_matrix_new[migration_matrix_new > 0])),
      col = "red", 
      lwd = 2)

# Add legend
legend(x = max(log10(migration_matrix)) * 0.7, y = 0.7,
       legend = c("观察值", "预测值"), 
       col = c("blue", "red"),
       lwd = 2,
       cex = 0.7,
       bty = "n")  # bty = "n" removes the box/frame around the legend



# Create dot plot for node strength comparison
# First, create a data frame for plotting
node_strength_df <- data.frame(
  Strength = c(NS, NS_new),
  Type = factor(rep(c("观察值", "预测值"), each = length(NS)))
)

# Add box plot on top of dot plot for additional visualization
ggplot(node_strength_df, aes(x = Type, y = Strength)) +
  geom_boxplot(alpha = 0.5) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  theme_minimal() +
  labs(title = "",
       x = "类型",
       y = "迁移总量") +
  theme(plot.title = element_text(hjust = 0.5))

# gravity model cant show the fat tail of the observed data


#=============================================================================
# 6. Network Analysis and Comparison
#=============================================================================
predicted_asymmetry <- calculate_asymmetry(migration_matrix_new)
predicted_globalcc <- calculate_clustering(migration_matrix_new)
