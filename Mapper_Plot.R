# Install pacman if not already installed
if (!require(pacman)) install.packages("pacman")

# Load all required packages
pacman::p_load(
  TDA, ggplot2, plotly, FNN, cluster, matrixStats, dbscan, igraph,
  rgl, mappeR, grid, ks, tidyr, devtools, fastcluster,
  DescTools, pROC, MASS, fclust, umap, mclust, NbClust,
  proxy, boot, pls, dplyr, infotheo, sigclust,
  randomForest, irr, accSDA, brant, RColorBrewer,
  factoextra, nnet, ordinalForest, survival, parallelDist,Rtsne
)



#########################################
#Load all the functions
#########################################

max_cluster_no<-5

filter_olda_beta <- function(X, Y) {
  X <- as.matrix(X)
  Y <- ordered(Y)
  
  # Fit Ordinal Discriminant Analysis (OLDA) via ordASDA
  fit <- ordASDA(
    Xt = X,
    Yt = as.numeric(Y) 
  )
  
  return(fit$beta[1:ncol(X), , drop=FALSE])
}



optimal_mapper_auto <- function(n, alpha = .5) {
  
  # Optimal number of intervals
  l_star <- ((8 * (1-alpha) * n) / alpha)^(1/5)
  l_star <- max(2, round(l_star))  # ensure at least 2 intervals
  
  # Compute S*, q*, lambda*
  S_star <- (2 * n) / (l_star + 1)
  q_star <- n / (l_star + 1)
  lambda_star <- q_star / S_star
  
  return(list(
    n = n,
    l_star = l_star,
    S_star = S_star,
    q_star = q_star,
    lambda_star = lambda_star
  ))
}



# Function to build intervals from sorted filter values
create_intervals <- function(sorted_values, alpha = 0.5) {
  
  n <- length(sorted_values)
  
  opt <- optimal_mapper_auto(n, alpha)
  l <- opt$l_star
  S <- floor(opt$S_star)

  # Adaptive stride
  # adjust to ensure full coverage (balanced cover)
  q_adapt <- (n - S) / (l - 1)
  
  ord <- seq_len(n)
  intervals_idx <- vector("list", l)
  
  for (j in 1:l) {
    start_idx <- 1 + floor((j - 1) * q_adapt)
    end_idx   <- min(start_idx + S - 1, n)
    
    intervals_idx[[j]] <- ord[start_idx:end_idx]
  }
  
  # Create display intervals
  intervals <- data.frame(
    start = sapply(intervals_idx, function(idx) min(sorted_values[idx])),
    end   = sapply(intervals_idx, function(idx) max(sorted_values[idx]))
  )
  
  intervals
}





# Define the Mapper algorithm (manual implementation)
generate_mapper <- function(X, filter_values, intervals, cluster_function) {
  n_intervals <- nrow(intervals)
  nodes <- list()
  edges <- list()
  
  # Loop over intervals
  for (i in seq_len(n_intervals)) {
    # Get points in the current interval
    in_interval <- (filter_values >= intervals$start[i]) & (filter_values <= intervals$end[i])
    subset_data <- X[in_interval, ]
    
    if (nrow(subset_data) > 0) {
      # Apply clustering function to the subset
      clusters <- cluster_function(subset_data)
      
      
      # Create nodes for each cluster
      cluster_ids <- unique(clusters)
      for (cluster_id in cluster_ids) {
        nodes[[length(nodes) + 1]] <- which(in_interval)[clusters == cluster_id]
      }
    }
  }
  
  # Build edges between overlapping clusters
  for (i in seq_along(nodes)) {
    for (j in seq_along(nodes)) {
      if (i < j && length(intersect(nodes[[i]], nodes[[j]])) > 0) {
        edges <- append(edges, list(c(i, j)))
      }
    }
  }
  
  # Return the graph as a list
  return(list(nodes = nodes, edges = edges))
}


# Clustering function(hiararchical clustering)
cluster_function <- function(data) {

  d_matrix <- parallelDist::parDist(as.matrix(data), method = "euclidean") 
  hc <- fastcluster::hclust(d_matrix, method = "ward.D2")
  
  k_vals <- 2:max_cluster_no
  sil_vals <- sapply(k_vals, function(k) {
    cl <- cutree(hc, k)
    mean(cluster::silhouette(cl, d_matrix)[, 3])
  })
  
  optimal_k <- if (all(is.na(sil_vals))) 1 else k_vals[which.max(sil_vals)]
  
  # Assign clusters
  cluster_assignments <- cutree(hc, k = optimal_k)
  return(cluster_assignments)
}


#######################################


#########################
# UCSC Xena Data
#########################

# 1. Load the final data with final selected variables
# 2. Load all the previous functions

file_path <- "C:/Users/ahsanm8/Desktop/Dr. Nitai_Final_code/Github_Material_Mapper/Real Data Analysis/UCSC Xena"

data<-load(file.path(file_path, "UCSC_Xena_data.RData"))

Y<-data[,1]
X<-data[,-1]



filter_values <- as.matrix(X) %*%  filter_olda_beta(X, Y)
sorted_values <- sort(filter_values)

intervals <- create_intervals(sorted_values,alpha=.1)



subsets <- vector("list", nrow(intervals))
row_indices_sub <- vector("list", nrow(intervals))

for (j in seq_len(nrow(intervals))) {
  
  idx_j <- which(filter_values >= intervals$start[j] &
                   filter_values <= intervals$end[j])
  
  row_indices_sub[[j]] <- idx_j
  subsets[[j]] <- X[idx_j, , drop = FALSE]
}




# Generate the Mapper graph
mapper_graph <- generate_mapper(X, filter_values, intervals, cluster_function)

# Scale vertex sizes by cluster sizes
node_sizes <- sapply(mapper_graph$nodes, length)


# Convert nodes and edges to a graph object
graph <- make_empty_graph(n = length(mapper_graph$nodes), directed = FALSE)
for (edge in mapper_graph$edges) {
  graph <- add_edges(graph, edge)
}

# Compute class proportions for each node for all 5 classes (1–3)
node_class_props <- t(sapply(mapper_graph$nodes, function(node_indices) {
  if (length(node_indices) == 0) {
    return(rep(0, 5))
  }
  
  prop <- table(factor(Y[node_indices], levels = 1:3))
  prop / sum(prop)
}))

colnames(node_class_props) <- paste0("Class", 1:3)

# 5-color palette for 5 classes
class_colors <- c(
  "#1A9850",  # Class 1
  "#FEE10B",  # Class 2
  "#A50026"   # Class 3
)


# Attach pie chart data to vertices
V(graph)$pie <- split(node_class_props, row(node_class_props))

# Node sizes
V(graph)$size <- sqrt(node_sizes) * 2   # adjust for visibility

# Colors applied per slice
V(graph)$pie.color <- split(matrix(class_colors, nrow = nrow(node_class_props),
                                   ncol = 3, byrow = TRUE),
                            row(node_class_props))

set.seed(1234)
png(
  filename = "Final_UCSC_Xena.png",
  width = 6000,
  height = 6000,
  res = 600
)
lay <- layout_with_fr(graph)

plot(
  graph,
  layout = lay,
  vertex.label = NA,
  vertex.shape = "pie",
  margin = 0
)

dev.off()





#########################
# PPMI
#########################

# 1. Load the final data with final selected variables
# 2. Load all the previous functions

file_path <- "C:/Users/ahsanm8/Desktop/Dr. Nitai_Final_code/Github_Material_Mapper/Real Data Analysis"

data<-load(file.path(file_path, "PPMI_data.RData"))

Y<-data[,1]
X<-data[,-1]




filter_values <- as.matrix(X) %*%  filter_olda_beta(X, Y)
sorted_values <- sort(filter_values)

intervals <- create_intervals(sorted_values,alpha=.1)



subsets <- vector("list", nrow(intervals))
row_indices_sub <- vector("list", nrow(intervals))

for (j in seq_len(nrow(intervals))) {
  
  idx_j <- which(filter_values >= intervals$start[j] &
                   filter_values <= intervals$end[j])
  
  row_indices_sub[[j]] <- idx_j
  subsets[[j]] <- X[idx_j, , drop = FALSE]
}




# Generate the Mapper graph
mapper_graph <- generate_mapper(X, filter_values, intervals, cluster_function)

# Scale vertex sizes by cluster sizes
node_sizes <- sapply(mapper_graph$nodes, length)


# Convert nodes and edges to a graph object
graph <- make_empty_graph(n = length(mapper_graph$nodes), directed = FALSE)
for (edge in mapper_graph$edges) {
  graph <- add_edges(graph, edge)
}

# Compute class proportions for each node for all 5 classes (1–3)
node_class_props <- t(sapply(mapper_graph$nodes, function(node_indices) {
  
  prop <- table(factor(Y[node_indices], levels = 0:2))
  prop / sum(prop)
}))

colnames(node_class_props) <- paste0("Class", 0:2)

# 5-color palette for 5 classes
class_colors <- c(
  "#1A9850",  # Class 1
  "#FEE10B",  # Class 2
  "#A50026"   # Class 3
)

# Build igraph object (your "graph" object should already exist)
# graph <- graph_from_edgelist(...)

# Attach pie chart data to vertices
V(graph)$pie <- split(node_class_props, row(node_class_props))

# Node sizes
V(graph)$size <- sqrt(node_sizes) * 2   # adjust for visibility

# Colors applied per slice
V(graph)$pie.color <- split(matrix(class_colors, nrow = nrow(node_class_props),
                                   ncol = 3, byrow = TRUE),
                            row(node_class_props))



set.seed(123)
png(
  filename = "PPMI_cutoff5.png",
  width = 6000,
  height = 6000,
  res = 600
)

plot(
  graph,
  vertex.label = NA,
  vertex.shape = "pie"
)

dev.off()





###############################
# t-SNE plot
###############################

df <- data

# t-SNE requires only numeric features → remove Y column
X <- df[, -1]

# Convert to matrix
X_mat <- as.matrix(X)

# Run t-SNE (perplexity rule: < N/3)
set.seed(123)
tsne_out <- Rtsne(X_mat, perplexity = 30, dims = 2, verbose = TRUE, check_duplicates = FALSE)

# Prepare plotting dataframe
plot_df <- data.frame(
  TSNE1 = tsne_out$Y[,1],
  TSNE2 = tsne_out$Y[,2],
  Y = factor(df$Y, ordered = TRUE)
)




###########
# UCSC Xena
###########

# Relabel Y to G2/G3/G4
plot_df$Y <- factor(plot_df$Y,
                    levels = c(1, 2, 3),
                    labels = c("G2", "G3", "G4"))

grade_colors <- c(
  "G2" = "#1A9850",
  "G3" = "#FEE10B",
  "G4" = "#A50026"
)

ggplot(plot_df, aes(x = TSNE1, y = TSNE2, color = Y)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = stage_colors) +
  theme_minimal(base_size = 14) +
  theme(    panel.grid = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            axis.text = element_blank(),        # remove axis numbers
            axis.ticks = element_blank())       # remove tick marks)





###########
# PPMI
###########
# Relabel Y to Stage 0 / Stage 1 / Stage 2+
plot_df$Y <- factor(plot_df$Y,
                    levels = c(0, 1, 2),
                    labels = c("Stage 0", "Stage 1", "Stage 2+"))

# Correct color mapping
stage_colors <- c(
  "Stage 0"  = "#1A9850",
  "Stage 1"  = "#FEE10B",
  "Stage 2+" = "#A50026"
)

ggplot(plot_df, aes(x = TSNE1, y = TSNE2, color = Y)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = stage_colors) +
  theme_minimal(base_size = 14) +
  theme(    panel.grid = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            axis.text = element_blank(),        # remove axis numbers
            axis.ticks = element_blank())       # remove tick marks)
