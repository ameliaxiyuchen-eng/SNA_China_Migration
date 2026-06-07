rm(list = ls())
library(readxl)
library(openxlsx)
library(dplyr)
library(geosphere) # For distance calculation using latitude and longitude
library(poweRlaw) 
library(igraph)
library(sna)
library(intergraph)
library(stats)
library(ggplot2)
library(wdnet)

setwd("/Users/xiyu/Desktop/network")


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



binary_matrix <- ifelse(migration_matrix > avg_link_weight  , 1, 0)
g_bi <- graph_from_adjacency_matrix(binary_matrix, mode = "directed", diag = FALSE)


in_degree <- degree(g_bi, mode = "in")
out_degree <- degree(g_bi, mode = "out")
ND = in_degree + out_degree

########################### Microcanonical graph ensembles ###########################


###################################### BINARY #####################################################

# Define the degree sequences for in-degrees and out-degrees
in_degrees <- in_degree
out_degrees <- out_degree
edge_density(g_bi) # actual edges/ largest possible edges

# Properties for the original graph
observed_BCC <- transitivity(g_bi, type = "local",isolates = c("zero"))
observed_ANND <- knn(g_bi)$knn  # Average nearest neighbor degree
observed_reci <- reciprocity(g_bi)

cor(ND,observed_BCC)
cor(ND, observed_ANND, use = "complete" )

# Create a directed graph with specific in-degree and out-degree sequences

random_graph<- sample_degseq(
  out.deg = out_degrees,
  in.deg = in_degrees,
  method = c("configuration")
)

plot(random_graph,
     vertex.size = 5,
     vertex.label = "",
     edge.width = 2,
     edge.arrow.size = 0)

edge_density(random_graph)


random.graphs.500.prob.list <- lapply(1:500, 
                                     function(x) 
                                       random_graph<- sample_degseq(
                                         out.deg = out_degrees,
                                         in.deg = in_degrees,
                                         method = c("configuration")
                                       ))

densities <- lapply(random.graphs.500.prob.list, 
                    function(x) edge_density(x)) # returns a list
unlist(densities)

calculate_properties <- function(graph) {
  BCC <- transitivity(graph, type = "local",isolates = c("zero"))  # Local clustering coefficient for each node
  ANND <- knn(graph)$knn      # Average nearest neighbor degree
  reci <- reciprocity(graph)
  return(list(BCC = BCC, ANND = ANND, reci = reci, centr = centr))
}


random_graph_properties <- lapply(random.graphs.500.prob.list, calculate_properties)

# Extract average properties for the random graphs
random_BCCs <- sapply(random_graph_properties, function(x) mean(x$BCC))
random_ANNDs <- sapply(random_graph_properties, function(x) x$ANND)
reci <-sapply(random_graph_properties, function(x) x$reci)


# Now, let's compare with the observed (original) properties
average_random_BCC <- mean(random_BCCs)
average_random_ANND <- mean(random_ANNDs)
average_random_reci <- mean(reci)

# Calculate the deviation (standard deviation or error)
sd_random_BCC <- sd(random_BCCs)
sd_random_degree <- sd(random_degrees)
sd_random_ANND <- sd(random_ANNDs)

# Plotting the results
library(ggplot2)

# Data for plotting
comparison_data <- data.frame(
  Metric = c("BCC", "ANND"),
  Observed = c(mean(observed_BCC), observed_ANND),
  Random_Avg = c(average_random_BCC, average_random_ANND),
  Random_SD = c(sd_random_BCC,sd_random_ANND)
)


# Plot the comparison with error bars for standard deviation
ggplot(comparison_data, aes(x = Metric)) +
  geom_bar(aes(y = Observed), stat = "identity", fill = "red", alpha = 0.7, width = 0.4, position = position_dodge(width = 0.8)) +
  geom_bar(aes(y = Random_Avg), stat = "identity", fill = "blue", alpha = 0.7, width = 0.4, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Random_Avg - Random_SD, ymax = Random_Avg + Random_SD), width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "Comparison of Original Graph with Random Graphs", 
       y = "Metric Value",
       x = "Metrics") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),  # Adjust axis labels
    axis.title = element_text(size = 14),  # Adjust axis title size
    plot.title = element_text(size = 16, hjust = 0.5)  # Adjust title size
    )
    

plot(density(random_BCCs), 
     xlim=c(0.2,0.7))
abline(v=mean(observed_BCC), lty=3) 
  
plot(density(random_ANNDs), 
     xlim=c(9,13))
abline(v= observed_ANND, lty=3)


###################################### WEIGHTED  #####################################################
migration_matrix
g_wei <- graph_from_adjacency_matrix(migration_matrix, mode = "directed", weighted = TRUE, diag = FALSE)

in_strength <- strength(g_wei, mode = "in")
out_strength <- strength(g_wei, mode = "out")

observed_WCC <- transitivity(g_wei, type = "local",isolates = c("zero"))
observed_ANNS <- knn(g_wei)$knn  # Average nearest neighbor degree
closeness_weighted_in <- closeness(g_wei, mode = "all", weights = E(g_wei)$weight, normalized = TRUE)

#wdnet testing
netwk <- igraph_to_wdnet(g_wei)
assortcoef(netwk)
WCC <- clustcoef(netwk, directed = TRUE,
          method = "Fagiolo",
          isolates = 0)$total$localcc

cor(NS,WCC)
cor(in_strength, out_strength)
##################### how to create weighted random graph? #############
etwk <- igraph_to_wdnet(random_graph)

assortcoef(etwk)
clustcoef(etwk, directed = TRUE,
          method = "Fagiolo")

plot(random_graph,
     vertex.size = 5,
     vertex.label = "",
     edge.width = 2,
     edge.arrow.size = 0)

random.graphs.500.list <- lapply(1:500, 
                                      function(x) 
                                        random_graph<- sample_degseq(
                                          out.deg = out_strength,
                                          in.deg = in_strength,
                                          method = c("configuration")
                                        ))

calculate_properties <- function(graph) {
  netwk <- igraph_to_wdnet(graph)
  assort <- assortcoef(netwk)
  cc <- clustcoef(netwk, directed = TRUE,
            method = "Fagiolo",
            isolates = 0)$outin
  return(list(assort = assort, cc = cc))
}

random_graph_properties <- lapply(random.graphs.500.prob.list, calculate_properties)

# Extract average properties for the random graphs
random_assorts <- sapply(random_graph_properties, function(x) x$assort)
random_ccs <- sapply(random_graph_properties, function(x) x$cc)

