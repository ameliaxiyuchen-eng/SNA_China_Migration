############################################ PROPERTY CHOOSING ##################################
rm(list = ls())

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


####################################################### WEIGHTED
############################### node strength distribution
g_5 <- graph_from_adjacency_matrix(matrix5, mode = "directed", weighted = TRUE, diag = FALSE)
g_6 <- graph_from_adjacency_matrix(matrix6, mode = "directed", weighted = TRUE, diag = FALSE)
g_7 <- graph_from_adjacency_matrix(matrix7, mode = "directed", weighted = TRUE, diag = FALSE)



ins_7 <- strength(g_7, mode = "in", weights = E(g_7)$weight)
out_7 <- strength(g_7, mode = "out", weights = E(g_7)$weight)
ns_7 <- strength(g_7, mode = "all")
rank(ns_7)

ins_6 <- strength(g_6, mode = "in", weights = E(g_6)$weight)
out_6 <- strength(g_6, mode = "out", weights = E(g_6)$weight)
ns_6 <- strength(g_6, mode = "all")
rank(ns_6)
ins_5 <- strength(g_5, mode = "in", weights = E(g_5)$weight)
out_5 <- strength(g_5, mode = "out", weights = E(g_5)$weight)
ns_5 <- strength(g_5, mode = "all")
rank(ns_5)

# Libraries
library(ggplot2)
library(ggrepel)
library(dplyr)


# Function to compute density for a given vector
compute_density <- function(strength, network, type) {
  dens <- density(strength) # Compute density
  data.frame(
    Strength = dens$x,
    Density = dens$y,
    Network = network,
    Type = type  )
}

# Compute densities for all node strength types
ins_density <- compute_density(ins_5, "1995-2000", "In-Strength") %>%
  bind_rows(compute_density(ins_6, "2005-2010", "In-Strength")) %>%
  bind_rows(compute_density(ins_7, "2015-2020", "In-Strength"))

out_density <- compute_density(out_5, "1995-2000", "Out-Strength") %>%
  bind_rows(compute_density(out_6, "2005-2010", "Out-Strength")) %>%
  bind_rows(compute_density(out_7, "2015-2020", "Out-Strength"))

ns_density <- compute_density(ns_5, "1995-2000", "Total-Strength") %>%
  bind_rows(compute_density(ns_6, "2005-2010", "Total-Strength")) %>%
  bind_rows(compute_density(ns_7, "2015-2020", "Total-Strength"))

# Combine all densities into a single data frame
strength_density_df <- bind_rows(ins_density, out_density, ns_density)

# Plot density curves for node strengths
ggplot(strength_density_df, aes(x = Strength, y = Density, color = Network)) +
  geom_line(size = 1) + # Density curves
  facet_wrap(~Type, scales = "free_x") + # Separate plots for In, Out, and Total
  scale_x_log10() + # Use log scale for node strength
  labs(
    title = "Node Strength Density Distributions by Network",
    x = "Node Strength (log scale)",
    y = "Density",
    color = "Time Period"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )


# Install poweRlaw package if not already installed
if (!requireNamespace("poweRlaw", quietly = TRUE)) {
  install.packages("poweRlaw")
}
library(poweRlaw)

# Example with in-strength (repeat for out-strength and total-strength)

# Define the data (replace ins_7 with your actual in-strength vector)
in_strength <- ns_7[ns_7 > 0] # Ensure positive values only

# Fit a power-law distribution
pl_in <- displ$new(in_strength)

# Estimate xmin (minimum value where power-law holds)
pl_in$setXmin(5000)

# Estimate the alpha (power-law exponent)
alpha_est <- estimate_pars(pl_in)

# Results
cat("In-Strength Power-Law Exponent (alpha):", alpha_est$pars, "\n")

# Summary of the fit
summary(pl_in)



################################## link wieght & power law 

hist(migration_matrix[migration_matrix > 0], breaks = 50, main = "Link Weight Distribution", xlab = "Weight")

###### in and out strength are correlated
cor(in_strength, out_strength)

# Link_weight powerlaw Estimated power-law exponent: 1.81947 
log_weights <- log(migration_matrix[migration_matrix > 0])
plot(log10(1:length(log_weights)), sort(log_weights, decreasing = TRUE), type = "l", 
     main = "Log-Log Plot of Link Weights", xlab = "Rank (Log Scale)", ylab = "Weight (Log Scale)")
pl_model <- displ$new(migration_matrix[migration_matrix > 0])
estimate_pars(pl_model)

link_weights <- as.vector(migration_matrix)
link_weights <- link_weights[link_weights > 0]
pl_model <- displ$new(link_weights)

pl_model$setXmin(1130) 
pl_model$setPars(estimate_pars(pl_model))

cat("Estimated power-law exponent:", pl_model$pars, "\n")
cat("Estimated xmin (cutoff):", pl_model$getXmin(), "\n")
plot(pl_model)
lines(pl_model, col = "red")

# Load necessary library
library(poweRlaw)

# Example data: replace with your migration data
link_weights <- as.vector(matrix7)  # Convert matrix to vector
link_weights <- link_weights[link_weights > 0]  # Filter positive values
mean(link_weights)
# Initialize the power-law model
pl_model <- displ$new(link_weights)

# Estimate the power-law exponent and other parameters
pl_model$setXmin(5000)  # Set xmin (minimum value for the power law)
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




#################################### link weight powerlaw distribution
# Load necessary libraries
library(ggplot2)

# Flatten the matrices into vectors
weights_1995_2000 <- unlist(migration5, use.names = FALSE)
weights_2005_2010 <- unlist(migration6, use.names = FALSE)
weights_2015_2020 <- unlist(migration7, use.names = FALSE)

# Combine all weights into a single data frame
weights_df <- data.frame(
  weight = c(weights_1995_2000, weights_2005_2010, weights_2015_2020),
  period = rep(c("1995-2000", "2005-2010", "2015-2020"),
               times = c(length(weights_1995_2000),
                         length(weights_2005_2010),
                         length(weights_2015_2020)))
)


weights_df <- weights_df[weights_df$weight > 0, ]

# Rank the weights within each time period
weights_df <- weights_df[order(weights_df$weight, decreasing = TRUE), ] # Sort by weight
weights_df$rank <- ave(weights_df$weight, weights_df$period, FUN = seq_along)

weights_df <- weights_df[order(weights_df$weight, decreasing = TRUE), ] # Sort by weight
weights_df$rank <- ave(weights_df$weight, weights_df$period, FUN = seq_along)


# Log-log plot
library(ggplot2)

ggplot(weights_df, aes(x = weight, y = rank, color = period)) +
  geom_point(size = 1) + # Use points instead of lines
  scale_x_log10() + # Log scale for weight
  scale_y_log10() + # Log scale for rank
  labs(
    title = "Log-Log Plot of Link Weight vs Rank",
    x = "Migrants",
    y = "Rank",
    color = "Time Period"
  ) +
  theme_minimal()


max_weight <- max(E(g_7)$weight)
largest_edge <- which(E(g_7)$weight == max_weight)
source_node <- ends(g_7, largest_edge)[1, 1]
target_node <- ends(g_7, largest_edge)[1, 2]


####################################################### Homophyly & Assortative
df <- read_csv("nodes_sum.csv")
df <- na.omit(df)
df <- df %>%
  mutate(lang_man = ifelse(lang_num <= 7, 1, lang_num))

gdp_2020 <- df$pgdp_2020
pop_2020 <- df$pop_2020
comlan <- df$COMLANG
lang_num <- df$lang_num
lang_man <- df$lang_man
climate <- factor(df$climate, 
                       levels = c("中温带", "暖温带", "北亚热带", "中亚热带", "南亚热带", "高原温带"),
                       labels = 1:6)
typeof(climate)
numeric_climate <- as.numeric(climate)
################# Assortativity coefficient is used do decide 
#whether nodes with the same features are tended to connect with each other more
### it ranges from [-1,1] and higher value represents higher homophyly

netwk <- igraph_to_wdnet(g_7)
################ disordered characteristic

########## 方言区 0.09848036 0.3253989
assor_comlan <- assortcoef(netwk, f1 = comlan, f2 = comlan )$'f1-f1'
assor_lang <- assortcoef(netwk, f1 = lang, f2 = lang )$'f1-f1'

assor_lang_man <- assortcoef(netwk, f1 = lang_man, f2 = lang_man )$'f1-f1'

########## 气候区 0.3862887
assor_comcli <- assortcoef(netwk, f1 = numeric_climate, f2 = numeric_climate )$'f1-f1'

################ ordered characteristic
########## gdp 相似的省之间关联更少 -0.1074815
gdp_2020
assor_gdp <- assortcoef(netwk, f1 = gdp_2020, f2 = gdp_2020 )$'f1-f1'

df$gdp_level <- cut(df$pgdp_2020,
                    breaks = 3,            # Divide the data into 4 intervals
                    labels = c("Low", "Medium", "High"), # Assign labels
                    include.lowest = TRUE)  # Include the lowest value in the first level
attr(df$gdp_level, "breaks")


######### 人口总量相似的省之间关联更少  -0.05803819
assor_pop <- assortcoef(netwk, f1 = pop_2020, f2 = pop_2020 )$'f1-f1'
df$pop_level <- cut(df$pop_2020,
                    breaks = 3,            # Divide the data into 4 intervals
                    labels = c("Low", "Medium", "High"), # Assign labels
                    include.lowest = TRUE)  # Include the lowest value in the first level


################ node strength
#### 人口流出量大的省份（起始地）与人口流出量大（目的地）的省份有更密切的联系 ：0.04446579
assor_outs <- assortcoef(netwk)$outout
#### 人口流出量大 （起始地）的省份与人口流入量大 （目的地）的省份有更密切的联系 ：0.1435518
assor_outins <- assortcoef(netwk)$outin
#### 人口流入量大 （起始地）的省份与人口流出量大 （目的地）的省份有更密切的联系 ：0.1069547
assor_inouts <- assortcoef(netwk)$inout
#### 人口流入量大的省份（起始地）与人口流入量大（目的地）的省份关联更少 ：-0.1823657
assor_ins <- assortcoef(netwk)$inin
#### 人口迁移量大的省份之间关联更少
rank(NS)
assort_ns <- assortcoef(netwk, f1 = ns_7, f2 = ns_7)

########## This code don't consider node_strength
assortativity(g_wei, values = gdp_2020,  directed = TRUE)
assortativity(g_wei, values = pop_2020,  directed = TRUE)



###################################### Centrality DONE
###################################### Comunity DONE

##################################### Clustering Coefficient
library(wdnet)
netwk_7 <- igraph_to_wdnet(g_7)
netwk_6 <- igraph_to_wdnet(g_6)
netwk_5 <- igraph_to_wdnet(g_5)

cc_7 <- clustcoef(netwk_7, method = "Fagiolo")
cc_6 <- clustcoef(netwk_6, method = "Fagiolo")
cc_5 <- clustcoef(netwk_5, method = "Fagiolo")
