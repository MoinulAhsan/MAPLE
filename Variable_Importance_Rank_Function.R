MAPLE_Importance_score <- function(data,secondary = FALSE,max_cv = 10) {
  
  ##############################
  ## Initialize objects
  ##############################
  
  mapper_Prediction_CV <- list()
  mapper_probability_prediction <- list()
  True_Y_CV <- list()
  
  mapper_qwk <- c()
  P_mapper_qwk_cv <- list()
  
  opt_intervals_CV <- c()
  
  accuracy_rate_CV <- list()
  
  max_cluster_no <- 5
  
  ##############################
  ## Cross-validation
  ##############################
  
    for(cv in 1:max_cv){
      
      print(paste("CV =", cv))
      
      set.seed(12345+cv)
      kk <- sample(1:nrow(data), size = round(0.2 * nrow(data))) #20% for test sample
      
      #################
      
      # Split the data into training and test sets
      train_data <- data[-kk, ]
      test_data <- data[kk, ]
      
      
      X<-as.matrix(train_data[,-1])
      Y<-train_data[, 1]
      
      n<-length(Y)
      
      
      
      
      
      
      
      ############Ordinal LDA ################
      # LDA will treat the class nominal not ordinal
      
      filter_olda_beta <- function(X, Y) {
        # Ensure predictors are matrix
        X <- as.matrix(X)
        # Convert response to ordered factor
        Y <- ordered(Y)
        # Fit Ordinal Discriminant Analysis (OLDA) via ordASDA
        fit <- ordASDA(
          Xt = X,
          Yt = as.numeric(Y),   # ordASDA needs numeric levels
        )
        
        return(fit$beta[1:ncol(X), , drop=FALSE])
      }
      
      
      filter_beta <- filter_olda_beta(X, Y)
      filter_values <-drop( X %*%  filter_beta)
      
      
      
      ###########################
      ## For balanced cover ##
      ###########################
      
      sorted_values <- sort(filter_values)
      
      create_intervals <- function(sorted_values, alpha = 0.5) {
        
        n <- length(sorted_values)
        
        opt <- optimal_mapper_auto(n, alpha)
        l <- opt$l_star
        S <- floor(opt$S_star)
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
      
      
      
      
      
      intervals_NUM<-c()
      mapper_Prediction_alp<-list()
      mapper_probability_prediction_alp<-list()
      SE_summary_CV_alp<-list()
      accuracy_rate<-c()
      
      
      saved_intervals<-list()
      saved_subsets<-list()
      saved_clusters<-list()
      saved_cross_interval_edges<-list()
      saved_row_indices_sub <- list()
      
      
      alpha<- seq(.9,.1,-.1)
      
      for(alp in 1:length(alpha)){
        intervals <- create_intervals(sorted_values,alpha= alpha[alp])
        
        intervals_NUM[alp]<-nrow(intervals)
        
        if (nrow(intervals) %in% intervals_NUM[-alp] ) next
        
        ######################
        
        
        
        
        # Initialize an empty list to store subsets
        subsets <- vector("list", nrow(intervals))
        row_indices_sub <- vector("list", nrow(intervals))
        
        for (j in seq_len(nrow(intervals))) {
          
          idx_j <- which(filter_values >= intervals$start[j] &
                           filter_values <= intervals$end[j])
          
          row_indices_sub[[j]] <- idx_j
          subsets[[j]] <- X[idx_j, , drop = FALSE]
        }
        
        
        
        #################
        # Use Hierarchical Clustering
        ################
        
        
        
        graphs <- list()
        clusters <- list()
        plot_list <- list()
        
        
        
        # Why hierarchical clustering (Ward.D2) ≈ spherical assumption
        #
        # Ward’s method minimizes within-cluster sum of squared Euclidean distances to the cluster centroid.
        # This implicitly means:
        #
        #   “Each cluster is best represented by its mean point, and all directions are equally important.
        
        # Store cluster assignments for each j for later cross-interval adjacency
        all_cluster_assignments <- list()
        
        for (j in 1:length(subsets)) {
          
          
          d_matrix <- parallelDist::parDist(subsets[[j]], method = "euclidean") #faster
          hc <- fastcluster::hclust(d_matrix, method = "ward.D2") #faster
          
          
          
          k_vals <- 2:max_cluster_no
          sil_vals <- sapply(k_vals, function(k) {
            cl <- cutree(hc, k)
            mean(cluster::silhouette(cl, d_matrix)[, 3])
          })
          
          optimal_k <- if (all(is.na(sil_vals))) 1 else k_vals[which.max(sil_vals)]
          
          
          
          
          cluster_assignments <- cutree(hc, k = optimal_k)
          
          
          
          
          all_cluster_assignments[[j]] <- data.frame(
            row_id = row_indices_sub[[j]],
            cluster = cluster_assignments,
            interval = j
          )
          
          adj_matrix <- outer(cluster_assignments, cluster_assignments, FUN = "==")
          diag(adj_matrix) <- 0
          g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
          clusters[[j]] <- components(g)
        }
        
        # ---- Build adjacency across intervals (both directions) ----
        cross_interval_edges <- data.frame(
          interval_from = integer(),
          cluster_from  = integer(),
          interval_to   = integer(),
          cluster_to    = integer(),
          stringsAsFactors = FALSE
        )
        
        for (j in 1:length(all_cluster_assignments)) {
          df_j <- all_cluster_assignments[[j]]
          
          # Look backward (j-1) if not the first interval
          if (j > 1) {
            df_prev <- all_cluster_assignments[[j - 1]]
            for (c1 in unique(df_j$cluster)) {
              rows_c1 <- df_j$row_id[df_j$cluster == c1]
              for (c2 in unique(df_prev$cluster)) {
                rows_c2 <- df_prev$row_id[df_prev$cluster == c2]
                if (length(intersect(rows_c1, rows_c2)) > 0) {
                  cross_interval_edges <- rbind(cross_interval_edges,
                                                data.frame(interval_from = j,
                                                           cluster_from  = c1,
                                                           interval_to   = j - 1,
                                                           cluster_to    = c2))
                }
              }
            }
          }
          
          # Look forward (j+1) if not the last interval
          if (j < length(all_cluster_assignments)) {
            df_next <- all_cluster_assignments[[j + 1]]
            for (c1 in unique(df_j$cluster)) {
              rows_c1 <- df_j$row_id[df_j$cluster == c1]
              for (c2 in unique(df_next$cluster)) {
                rows_c2 <- df_next$row_id[df_next$cluster == c2]
                if (length(intersect(rows_c1, rows_c2)) > 0) {
                  cross_interval_edges <- rbind(cross_interval_edges,
                                                data.frame(interval_from = j,
                                                           cluster_from  = c1,
                                                           interval_to   = j + 1,
                                                           cluster_to    = c2))
                }
              }
            }
          }
        }
        
        
        cross_interval_edges$node_from <- paste0("I", cross_interval_edges$interval_from, "_C", cross_interval_edges$cluster_from)
        cross_interval_edges$node_to   <- paste0("I", cross_interval_edges$interval_to, "_C", cross_interval_edges$cluster_to)
        
        nodes <- unique(c(cross_interval_edges$node_from, cross_interval_edges$node_to))
        adj_matrix_global <- matrix(0, nrow = length(nodes), ncol = length(nodes),
                                    dimnames = list(nodes, nodes))
        
        for (i in 1:nrow(cross_interval_edges)) {
          adj_matrix_global[cross_interval_edges$node_from[i],
                            cross_interval_edges$node_to[i]] <- 1
          adj_matrix_global[cross_interval_edges$node_to[i],
                            cross_interval_edges$node_from[i]] <- 1  # undirected
        }
        
        
        
        
        
        
        
        ############################
        #prediction
        ############################
        
        
        new_point<-as.matrix(test_data[,-1])
        true_Y<-test_data[, 1]
        
        
        filter_value <-drop( new_point %*%  filter_beta)
        
        
        
        
        interval_index <- list()
        for (i in 1:length(filter_value)) {
          if (filter_value[i] < min(intervals$start)) {
            # Assign to the first interval if the filter value is smaller than the starting point
            interval_index[[i]] <- 1
          } else if (filter_value[i] > max(intervals$end)) {
            # Assign to the last interval if the filter value is larger than the maximum endpoint
            interval_index[[i]] <- nrow(intervals)
          } else {
            # Otherwise, find the appropriate interval
            interval_index[[i]] <- which(
              filter_value[i] >= intervals$start &
                filter_value[i] <= intervals$end
            )
          }
        }
        
        
        
        # red=1, blue=0
        
        weighted_average_prediction<-c()
        weighted_probability_prediction<-c()
        weighted_average_prediction_prob<-list()
        closest_cluster_for_j_and_j1<-c()
        dist_closest_cluster_idx<-c()
        
        SE_store<-list()
        
        
        pred_class_QWK<-c()
        
        
        for(i in 1:length(interval_index)){
          
          if(length(interval_index[[i]])==1){
            
            interval_index_pos<-interval_index[[i]]
            # Step 3: Use the subset corresponding to the interval
            subset_new_point <- subsets[[interval_index_pos]]
            
            new_point_pos<-new_point[i,]
            
            
            
            
            
            # Find the representative points (centroids) for each cluster in the current interval
            cluster_centroids <- sapply(1:length(clusters[[interval_index_pos]]$csize), function(i) {
              # Extract the points that belong to the i-th cluster
              cluster_points <- subsets[[interval_index_pos]][clusters[[interval_index_pos]]$membership == i, ]
              
              # Calculate the centroid of the i-th cluster (mean of the points)
              if(length(cluster_points)==ncol(X)){    
                centroid <- cluster_points
                return(centroid)
              }else{
                centroid <- colMeans(cluster_points)
                return(centroid)
              }
              
            })
            
            
            
            # Calculate the distance from the new point to each cluster centroid
            distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
              # Calculate the Euclidean distance between the new point and the centroid
              dist(rbind(centroid, new_point_pos))
            })
            
            # Find the index of the closest cluster
            distances_to_centroids[which(clusters[[interval_index_pos]]$csize==1)]<-99999 # avoid taking the centroid that has only one cluster point. puting large number so that minimum does not take that
            closest_cluster_idx <- which.min(distances_to_centroids)
            
            # Assign the new point to the closest cluster
            closest_cluster <- closest_cluster_idx
            
            
            
            
            cluster_points <- subset_new_point[clusters[[interval_index_pos]]$membership == closest_cluster, ]
            
            # Find the row indices in X that correspond to each row in cluster_points
            row_indices <- row_indices_sub[[interval_index_pos]][clusters[[interval_index_pos]]$membership == closest_cluster]
            
            
            # Get the corresponding response values for the selected cluster points
            response_values <-Y[row_indices]
            
            
            
            if( secondary ==T  ){
              
              
              
              
              #################
              # cross_interval_edges
              
              closest_node<-paste0("I", interval_index_pos, "_C", closest_cluster)
              
              
              connected_rows <- subset(cross_interval_edges, node_from   == closest_node)
              
              node_to_list <- connected_rows$node_to
              
              all_nodes <- lapply(node_to_list, function(node) {
                # Extract interval and cluster from the node string
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster <- as.numeric(parts[2])
                
                # Subset the nodes automatically
                subsets[[interval]][clusters[[interval]]$membership == cluster, ]
              })
              
              # Combine everything into one data frame
              all_nodes_combined <- unique(do.call(rbind, all_nodes))
              
              
              
              
              # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
              cluster_points <- subset_new_point[clusters[[interval_index_pos]]$membership == closest_cluster, ]
              
              
              
              # Convert matrices to data frames
              all_nodes_df <- as.data.frame(all_nodes_combined)
              cluster_points_df <- as.data.frame(cluster_points)
              
              
              cluster_points<- as.matrix(unique(rbind(all_nodes_df, cluster_points_df)))
              
              
              
              #################
              
              
              
              # Find the row indices in X that correspond to each row in cluster_points
              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })
              
              
              # Get the corresponding response values for the selected cluster points
              response_values <-Y[row_indices]
              
              
              neighbor_2nd_count<-neighbor_2nd_count+1
              
              
            }
            
            
            
            combined_df <- data.frame(response = response_values, cluster_points)
            
            # Convert your vector into a one-row data frame
            new_point_df <- as.data.frame(as.list(new_point_pos))
            
            
            
            ###################################
            
            # Compute Euclidean distances from the new point
            distances <- apply(combined_df[, -1], 1, function(x) sum((x - unlist(new_point_df))^2))  # need to write why squared distance

            
            # Extract neighbor responses and distances
            neighbor_responses <- combined_df$response
            
            # Compute weights (inverse distance)
            weights <- 1 / (distances + 1e-6)  # avoid divide by zero
            
            # Normalize so weights sum to 1
            weights <- weights / sum(weights)
            
            
            # Compute weighted probabilities for each ordinal level
            categories <- levels(combined_df$response)
            
            
            # Compute cumulative weighted probabilities
            cum_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses <= c]) / sum(weights)
            })
            
            
            
            # Compute non-cumulative (category-specific) probabilities
            cat_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses == c]) / sum(weights)
            })
            
            
            # Posterior median decision rule (ordinal)
            median_index <- which(cum_probs >= 0.5)[1]
            
            # Ordinal responses have meaningful order but not equal numeric spacing between labels.
            # this is the minimizer of expected absolute error,
            pred_class <- categories[median_index]
            
            
            # Output
            pred_probs<- cat_probs
            
            
            
            
            weighted_average_prediction_prob[[i]] <- pred_probs
            
            weighted_average_prediction[i] <- as.numeric(as.character(pred_class))
            
            
          }else{
            
            ############
            
            interval_index_pos<-interval_index[[i]]
            new_point_pos<-new_point[i,]
            
            for(ind in 1:length(interval_index_pos)){
              
              index<-interval_index_pos[ind]
              
              # Step 3: Use the subset corresponding to the interval
              subset_new_point <- subsets[[index]]
              
              
              # Find the representative points (centroids) for each cluster in the current interval
              cluster_centroids <- sapply(1:length(clusters[[index]]$csize), function(i) {
                # Extract the points that belong to the i-th cluster
                cluster_points <- subsets[[index]][clusters[[index]]$membership == i, ]
                
                # Calculate the centroid of the i-th cluster (mean of the points)
                if(length(cluster_points)==ncol(X)){
                  centroid <- cluster_points
                  return(centroid)
                }else{
                  centroid <- colMeans(cluster_points)
                  
                  return(centroid)
                }
                
              })
              
              
              
              # Calculate the distance from the new point to each cluster centroid
              distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
                # Calculate the Euclidean distance between the new point and the centroid
                dist(rbind(centroid, new_point_pos))
              })
              
              # Find the index of the closest cluster
              distances_to_centroids[which(clusters[[index]]$csize==1)]<-99999 # avoid taking the centroid that has only one cluster point. puting large number so that minimum does not take that
              closest_cluster_idx <- which.min(distances_to_centroids)
              
              dist_closest_cluster_idx[ind] <- min(distances_to_centroids)
              
              # Assign the new point to the closest cluster
              closest_cluster <- closest_cluster_idx
              
              
              closest_cluster_for_j_and_j1[ind]<-closest_cluster
              
            }
            
            
            ############
            
            
            
            
            
            ###################
            
            closest_node<-paste0("I", interval_index_pos, "_C", closest_cluster_for_j_and_j1)
            
            
            
            # connected_rows <- subset(cross_interval_edges, node_from   == closest_node)
            connected_rows <- subset(cross_interval_edges, node_from %in% closest_node)
            
            
            
            
            # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
            # Initialize empty data frame
            cluster_points <- data.frame()
            
            # Loop over the pairs
            for (ii in seq_along(interval_index_pos)) {
              jj <- interval_index_pos[ii]
              cluster_id <- closest_cluster_for_j_and_j1[ii]
              
              subset_new_point<-subsets[[jj]]  ## this need to fixed in other code
              
              
              # Extract rows for this cluster
              tmp <- subset_new_point[clusters[[jj]]$membership == cluster_id, ]
              
              if (nrow(cluster_points) == 0) {
                # First iteration, just take all rows
                cluster_points <- tmp
              } else {
                # Keep only rows not already in cluster_points
                tmp_unique <- tmp[!apply(tmp, 1, function(row) {
                  any(apply(cluster_points, 1, function(existing) all(existing == row)))
                }), ]
                
                # Append unique rows
                cluster_points <- rbind(cluster_points, tmp_unique)
              }
            }
            
            
            
            
            # Find the row indices in X that correspond to each row in cluster_points
            row_indices_list <- list()
            
            for(ind in seq_along(interval_index_pos)){
              
              index <- interval_index_pos[ind]
              
              row_indices_list[[ind]] <-
                row_indices_sub[[index]][clusters[[index]]$membership == closest_cluster_for_j_and_j1[ind]]
              
            }
            
            row_indices <- unique(unlist(row_indices_list))
            
            
            # Get the corresponding response values for the selected cluster points
            response_values <-Y[row_indices]
            
            
            
            
            if( secondary ==T  ){
              
              
              
              node_to_list <- connected_rows$node_to
              
              different_node<-unique( setdiff(node_to_list, closest_node) )
              
              all_nodes <- lapply(different_node, function(node) {
                # Extract interval and cluster from the node string
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster <- as.numeric(parts[2])
                
                # Subset the nodes automatically
                subsets[[interval]][clusters[[interval]]$membership == cluster, ]
              })
              
              # Combine everything into one data frame
              all_nodes_combined <- unique(do.call(rbind, all_nodes))  #???? check later to avoid unique
              
              
              
              
              
              # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
              # Initialize empty data frame
              cluster_points <- data.frame()
              
              # Loop over the pairs
              for (ii in seq_along(interval_index_pos)) {
                jj <- interval_index_pos[ii]
                cluster_id <- closest_cluster_for_j_and_j1[ii]
                
                subset_new_point<-subsets[[jj]]  ## this need to fixed in other code
                
                
                # Extract rows for this cluster
                tmp <- subset_new_point[clusters[[jj]]$membership == cluster_id, ]
                
                if (nrow(cluster_points) == 0) {
                  # First iteration, just take all rows
                  cluster_points <- tmp
                } else {
                  # Keep only rows not already in cluster_points
                  tmp_unique <- tmp[!apply(tmp, 1, function(row) {
                    any(apply(cluster_points, 1, function(existing) all(existing == row)))
                  }), ]
                  
                  # Append unique rows
                  cluster_points <- rbind(cluster_points, tmp_unique)
                }
              }
              
              
              
              
              # Convert matrices to data frames
              all_nodes_df <- as.data.frame(all_nodes_combined)
              cluster_points_df <- as.data.frame(cluster_points)
              
              
              cluster_points<- as.matrix(unique(rbind(all_nodes_df, cluster_points_df)))
              
              
              
              # Find the row indices in X that correspond to each row in cluster_points
              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })
              
              
              # Get the corresponding response values for the selected cluster points
              response_values <-Y[row_indices]
              
              
              neighbor_2nd_count<-neighbor_2nd_count+1
              
              
            }
            
            
            
            
            combined_df <- data.frame(response = response_values, cluster_points)
            
            # Convert your vector into a one-row data frame
            new_point_df <- as.data.frame(as.list(new_point_pos))
            
            
            
            
            ###################################
            
            # Compute Euclidean distances from the new point
            distances <- apply(combined_df[, -1], 1, function(x) sum((x - unlist(new_point_df))^2))

            
            # Extract neighbor responses and distances
            neighbor_responses <- combined_df$response
            
            # Compute weights (inverse distance)
            weights <- 1 / (distances + 1e-6)  # avoid divide by zero
            
            # Normalize so weights sum to 1
            weights <- weights / sum(weights)
            
            
            # Compute weighted probabilities for each ordinal level
            categories <- levels(combined_df$response)
            
            
            # Compute cumulative weighted probabilities
            cum_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses <= c]) / sum(weights)
            })
            
            
            # Compute non-cumulative (category-specific) probabilities
            cat_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses == c]) / sum(weights)
            })
            
            
            # Posterior median decision rule (ordinal)
            median_index <- which(cum_probs >= 0.5)[1]
            
            
            # Ordinal responses have meaningful order but not equal numeric spacing between labels.
            # this is the minimizer of expected absolute error,
            pred_class <- categories[median_index]
            
            
            # Output
            pred_probs<- cat_probs
            
            
            
            weighted_average_prediction_prob[[i]] <- pred_probs
            
            weighted_average_prediction[i] <- as.numeric(as.character(pred_class))
            
            ###################
            
            
            
          }
          
        }
        
        
        mapper_Prediction_alp[[alp]]<-weighted_average_prediction
        mapper_probability_prediction_alp[[alp]]<-weighted_average_prediction_prob
        
        

        
        saved_intervals[[alp]] <- intervals
        saved_subsets[[alp]] <- subsets
        saved_clusters[[alp]] <- clusters
        saved_cross_interval_edges[[alp]] <- cross_interval_edges
        saved_row_indices_sub[[alp]] <- row_indices_sub
        
        levels_Y <- as.numeric(levels(true_Y))
        pred_class_mapper <- factor(weighted_average_prediction, levels = levels_Y)
        y_numeric <- as.numeric(as.character(true_Y))
        accuracy_rate[alp]<- kappa2(data.frame(y_numeric, as.numeric(as.character(pred_class_mapper))),
                                    weight = "squared")$value
        
        
        
      }
      
      
      Position_accuracy_rate<- max(which(accuracy_rate == max(accuracy_rate, na.rm = TRUE)))
      # intervals_NUM
      
      mapper_Prediction_CV[[cv]]<-mapper_Prediction_alp[[Position_accuracy_rate]]
      mapper_probability_prediction[[cv]]<-mapper_probability_prediction_alp[[Position_accuracy_rate]]
      True_Y_CV[[cv]]<-true_Y
      
      opt_intervals_CV[cv]<-intervals_NUM[Position_accuracy_rate]
      accuracy_rate_CV[[cv]]<-accuracy_rate
      
      
      
      mat_mapper<- do.call(rbind,mapper_probability_prediction[[cv]])
      colnames(mat_mapper) <- levels(Y)
      
      
      y_numeric<-as.numeric(as.character(true_Y))
      qwk_mapper <- kappa2(data.frame(y_numeric, as.numeric(as.character(mapper_Prediction_CV[[cv]]))),
                           weight = "squared")$value
      
      mapper_qwk[cv]<-qwk_mapper
      
      
      
      
      
      intervals <- saved_intervals[[Position_accuracy_rate]]
      subsets <- saved_subsets[[Position_accuracy_rate]]
      clusters <- saved_clusters[[Position_accuracy_rate]]
      cross_interval_edges <- saved_cross_interval_edges[[Position_accuracy_rate]]
      row_indices_sub <- saved_row_indices_sub[[Position_accuracy_rate]]
      
      
      
      
      
      ############################
      #prediction with permuting
      ############################
      
      
      
      variable_importance <- numeric(ncol(X))
      
      P_mapper_qwk<-c()
      
      new_point_original<-as.matrix(test_data[,-1])
      
      
      for (r in 1:ncol(X)) {
        
        
        print(paste("R =", r))
        
        # Permute variable r in test set
        new_point<-new_point_original
        
        new_point[, r] <- sample(new_point[, r], replace = FALSE)

        
        filter_value <- drop(new_point %*% filter_beta) #For variable importance of the whole Mapper model, topology and prediction
        
        true_Y<-test_data[, 1]
        
        
        
        
        
        
        
        
        
        interval_index <- list()
        for (i in 1:length(filter_value)) {
          if (filter_value[i] < min(intervals$start)) {
            # Assign to the first interval if the filter value is smaller than the starting point
            interval_index[[i]] <- 1
          } else if (filter_value[i] > max(intervals$end)) {
            # Assign to the last interval if the filter value is larger than the maximum endpoint
            interval_index[[i]] <- nrow(intervals)
          } else {
            # Otherwise, find the appropriate interval
            interval_index[[i]] <- which(
              filter_value[i] >= intervals$start &
                filter_value[i] <= intervals$end
            )
          }
        }
        
        
        
        # red=1, blue=0
        
        P_weighted_average_prediction<-c()
        P_weighted_probability_prediction<-c()
        P_weighted_average_prediction_prob<-list()
        closest_cluster_for_j_and_j1<-c()
        dist_closest_cluster_idx<-c()
        
        SE_store<-list()
        
        
        for(i in 1:length(interval_index)){
          
          if(length(interval_index[[i]])==1){
            
            interval_index_pos<-interval_index[[i]]
            # Step 3: Use the subset corresponding to the interval
            subset_new_point <- subsets[[interval_index_pos]]
            
            new_point_pos<-new_point[i,]
            
            
            
            
            
            # Find the representative points (centroids) for each cluster in the current interval
            cluster_centroids <- sapply(1:length(clusters[[interval_index_pos]]$csize), function(i) {
              # Extract the points that belong to the i-th cluster
              cluster_points <- subsets[[interval_index_pos]][clusters[[interval_index_pos]]$membership == i, ]
              
              # Calculate the centroid of the i-th cluster (mean of the points)
              if(length(cluster_points)==ncol(X)){    
                centroid <- cluster_points
                return(centroid)
              }else{
                centroid <- colMeans(cluster_points)
                return(centroid)
              }
              
            })
            
            
            
            # Calculate the distance from the new point to each cluster centroid
            distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
              # Calculate the Euclidean distance between the new point and the centroid
              dist(rbind(centroid, new_point_pos))
            })
            
            # Find the index of the closest cluster
            distances_to_centroids[which(clusters[[interval_index_pos]]$csize==1)]<-99999 # avoid taking the centroid that has only one cluster point. puting large number so that minimum does not take that
            closest_cluster_idx <- which.min(distances_to_centroids)
            
            # Assign the new point to the closest cluster
            closest_cluster <- closest_cluster_idx
            
            
            
            
            cluster_points <- subset_new_point[clusters[[interval_index_pos]]$membership == closest_cluster, ]
            
            # Find the row indices in X that correspond to each row in cluster_points
            row_indices <- row_indices_sub[[interval_index_pos]][clusters[[interval_index_pos]]$membership == closest_cluster]
      
            
            # Get the corresponding response values for the selected cluster points
            response_values <-Y[row_indices]
            
            
            
            if( secondary ==T  ){
              
              
              
              
              #################
              # cross_interval_edges
              
              closest_node<-paste0("I", interval_index_pos, "_C", closest_cluster)
              
              
              connected_rows <- subset(cross_interval_edges, node_from   == closest_node)
              
              node_to_list <- connected_rows$node_to
              
              all_nodes <- lapply(node_to_list, function(node) {
                # Extract interval and cluster from the node string
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster <- as.numeric(parts[2])
                
                # Subset the nodes automatically
                subsets[[interval]][clusters[[interval]]$membership == cluster, ]
              })
              
              # Combine everything into one data frame
              all_nodes_combined <- unique(do.call(rbind, all_nodes))
              
              
              
              
              # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
              cluster_points <- subset_new_point[clusters[[interval_index_pos]]$membership == closest_cluster, ]
              
              
              
              # Convert matrices to data frames
              all_nodes_df <- as.data.frame(all_nodes_combined)
              cluster_points_df <- as.data.frame(cluster_points)
              
              
              
              cluster_points<- as.matrix(unique(rbind(all_nodes_df, cluster_points_df)))
              
              
              
              #################
              
              
              
              # Find the row indices in X that correspond to each row in cluster_points
              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })
              
              
              # Get the corresponding response values for the selected cluster points
              response_values <-Y[row_indices]
              
              
              neighbor_2nd_count<-neighbor_2nd_count+1
              
              
            }
            
            
            
            combined_df <- data.frame(response = response_values, cluster_points)
            
            # Convert your vector into a one-row data frame
            new_point_df <- as.data.frame(as.list(new_point_pos))
            
            
            
            ###################################
            
            # Compute Euclidean distances from the new point
            distances <- apply(combined_df[, -1], 1, function(x) sum((x - unlist(new_point_df))^2))  # need to write why squared distance
            # distances <- apply(combined_df[, -1], 1, function(x) sqrt(sum((x - unlist(new_point_df))^2)))
            
            
            # Extract neighbor responses and distances
            neighbor_responses <- combined_df$response
            
            # Compute weights (inverse distance)
            weights <- 1 / (distances + 1e-6)  # avoid divide by zero
            
            # Normalize so weights sum to 1
            weights <- weights / sum(weights)
            
            
            # Compute weighted probabilities for each ordinal level
            categories <- levels(combined_df$response)
            
            
            # Compute cumulative weighted probabilities
            cum_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses <= c]) / sum(weights)
            })
            
            
            
            # Compute non-cumulative (category-specific) probabilities
            cat_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses == c]) / sum(weights)
            })
            
            
            # Posterior median decision rule (ordinal)
            median_index <- which(cum_probs >= 0.5)[1]
            
            
            # Ordinal responses have meaningful order but not equal numeric spacing between labels.
            # this is the minimizer of expected absolute error,
            pred_class <- categories[median_index]
            
            # Output
            pred_probs<- cat_probs
            
            
            
            
            P_weighted_average_prediction_prob[[i]] <- pred_probs
            
            P_weighted_average_prediction[i] <- as.numeric(as.character(pred_class))
            
            
          }else{
            
            ############
            
            interval_index_pos<-interval_index[[i]]
            new_point_pos<-new_point[i,]
            
            for(ind in 1:length(interval_index_pos)){
              
              index<-interval_index_pos[ind]
              
              # Step 3: Use the subset corresponding to the interval
              subset_new_point <- subsets[[index]]
              
              
              # Find the representative points (centroids) for each cluster in the current interval
              cluster_centroids <- sapply(1:length(clusters[[index]]$csize), function(i) {
                # Extract the points that belong to the i-th cluster
                cluster_points <- subsets[[index]][clusters[[index]]$membership == i, ]
                
                # Calculate the centroid of the i-th cluster (mean of the points)
                if(length(cluster_points)==ncol(X)){
                  centroid <- cluster_points
                  return(centroid)
                }else{
                  centroid <- colMeans(cluster_points)
                  
                  return(centroid)
                }
                
              })
              
              
              
              # Calculate the distance from the new point to each cluster centroid
              distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
                # Calculate the Euclidean distance between the new point and the centroid
                dist(rbind(centroid, new_point_pos))
              })
              
              # Find the index of the closest cluster
              distances_to_centroids[which(clusters[[index]]$csize==1)]<-99999 # avoid taking the centroid that has only one cluster point. puting large number so that minimum does not take that
              closest_cluster_idx <- which.min(distances_to_centroids)
              
              dist_closest_cluster_idx[ind] <- min(distances_to_centroids)
              
              # Assign the new point to the closest cluster
              closest_cluster <- closest_cluster_idx
              
              
              closest_cluster_for_j_and_j1[ind]<-closest_cluster
              
            }
            
            
            ############
            
            
            
            closest_node<-paste0("I", interval_index_pos, "_C", closest_cluster_for_j_and_j1)
            
            
            
            # connected_rows <- subset(cross_interval_edges, node_from   == closest_node)
            connected_rows <- subset(cross_interval_edges, node_from %in% closest_node)
            
            
            
            
            # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
            # Initialize empty data frame
            cluster_points <- data.frame()
            
            # Loop over the pairs
            for (ii in seq_along(interval_index_pos)) {
              jj <- interval_index_pos[ii]
              cluster_id <- closest_cluster_for_j_and_j1[ii]
              
              subset_new_point<-subsets[[jj]]  ## this need to fixed in other code
              
              
              # Extract rows for this cluster
              tmp <- subset_new_point[clusters[[jj]]$membership == cluster_id, ]
              
              if (nrow(cluster_points) == 0) {
                # First iteration, just take all rows
                cluster_points <- tmp
              } else {
                # Keep only rows not already in cluster_points
                tmp_unique <- tmp[!apply(tmp, 1, function(row) {
                  any(apply(cluster_points, 1, function(existing) all(existing == row)))
                }), ]
                
                # Append unique rows
                cluster_points <- rbind(cluster_points, tmp_unique)
              }
            }
            
            
            
            
            # Find the row indices in X that correspond to each row in cluster_points
            row_indices_list <- list()
            
            for(ind in seq_along(interval_index_pos)){
              
              index <- interval_index_pos[ind]
              
              row_indices_list[[ind]] <-
                row_indices_sub[[index]][clusters[[index]]$membership == closest_cluster_for_j_and_j1[ind]]
              
            }
            
            row_indices <- unique(unlist(row_indices_list))
            
            
            # Get the corresponding response values for the selected cluster points
            response_values <-Y[row_indices]
            
            
            
            
            if( secondary ==T  ){
              
              
              
              node_to_list <- connected_rows$node_to
              
              different_node<-unique( setdiff(node_to_list, closest_node) )
              
              all_nodes <- lapply(different_node, function(node) {
                # Extract interval and cluster from the node string
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster <- as.numeric(parts[2])
                
                # Subset the nodes automatically
                subsets[[interval]][clusters[[interval]]$membership == cluster, ]
              })
              
              # Combine everything into one data frame
              all_nodes_combined <- unique(do.call(rbind, all_nodes))  
              
              
              
              
              
              # Get the subset of points corresponding to the closest cluster (1st degree neighbor)
              # Initialize empty data frame
              cluster_points <- data.frame()
              
              # Loop over the pairs
              for (ii in seq_along(interval_index_pos)) {
                jj <- interval_index_pos[ii]
                cluster_id <- closest_cluster_for_j_and_j1[ii]
                
                subset_new_point<-subsets[[jj]]  ## this need to fixed in other code
                
                
                # Extract rows for this cluster
                tmp <- subset_new_point[clusters[[jj]]$membership == cluster_id, ]
                
                if (nrow(cluster_points) == 0) {
                  # First iteration, just take all rows
                  cluster_points <- tmp
                } else {
                  # Keep only rows not already in cluster_points
                  tmp_unique <- tmp[!apply(tmp, 1, function(row) {
                    any(apply(cluster_points, 1, function(existing) all(existing == row)))
                  }), ]
                  
                  # Append unique rows
                  cluster_points <- rbind(cluster_points, tmp_unique)
                }
              }
              
              
              
              
              # Convert matrices to data frames
              all_nodes_df <- as.data.frame(all_nodes_combined)
              cluster_points_df <- as.data.frame(cluster_points)
              
              
              cluster_points<- as.matrix(unique(rbind(all_nodes_df, cluster_points_df)))
              
              
              
              # Find the row indices in X that correspond to each row in cluster_points
              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })
              
              
              # Get the corresponding response values for the selected cluster points
              response_values <-Y[row_indices]
              
              
              neighbor_2nd_count<-neighbor_2nd_count+1
              
              
            }
            
            
            combined_df <- data.frame(response = response_values, cluster_points)
            
            # Convert your vector into a one-row data frame
            new_point_df <- as.data.frame(as.list(new_point_pos))
            
            
            
            
            ###################################
            
            # Compute Euclidean distances from the new point
            distances <- apply(combined_df[, -1], 1, function(x) sum((x - unlist(new_point_df))^2))

            
            # Extract neighbor responses and distances
            neighbor_responses <- combined_df$response
            
            # Compute weights (inverse distance)
            weights <- 1 / (distances + 1e-6)  # avoid divide by zero
            
            # Normalize so weights sum to 1
            weights <- weights / sum(weights)
            
            
            # Compute weighted probabilities for each ordinal level
            categories <- levels(combined_df$response)
            
            
            # Compute cumulative weighted probabilities
            cum_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses <= c]) / sum(weights)
            })
            
            
            # Compute non-cumulative (category-specific) probabilities
            cat_probs <- sapply(categories, function(c) {
              sum(weights[neighbor_responses == c]) / sum(weights)
            })
            
            
            # Posterior median decision rule (ordinal)
            median_index <- which(cum_probs >= 0.5)[1]
            
            # Ordinal responses have meaningful order but not equal numeric spacing between labels.
            # this is the minimizer of expected absolute error,
            pred_class <- categories[median_index]
            
            
            # Output
            pred_probs<- cat_probs
            
            
            
            P_weighted_average_prediction_prob[[i]] <- pred_probs
            
            P_weighted_average_prediction[i] <- as.numeric(as.character(pred_class))
            
            ###################
            
            
            
          }
          
        }
        
        
        
        #########################
        ## K fold cross validation
        #########################
        
        
        P_mat_mapper<- do.call(rbind,P_weighted_average_prediction_prob)
        colnames(P_mat_mapper) <- levels(Y)
        
        
        
        
        y_numeric<-as.numeric(as.character(true_Y))
        P_qwk_mapper <- kappa2(data.frame(y_numeric, as.numeric(as.character(P_weighted_average_prediction))),
                               weight = "squared")$value
        
        P_mapper_qwk[r]<-P_qwk_mapper
        
        
        
        
        
      }
      
      P_mapper_qwk_cv[[cv]]<-P_mapper_qwk
    }
    
    
    
    
    
    # Compute % decrease in accuracy
    percent_decrease <- sapply(1:max_cv, function(i) {
      (mapper_qwk[i] - P_mapper_qwk_cv[[i]]) / mapper_qwk[i] * 100
    })
    
    
    mean_decrease_accuracy <- rowMeans(percent_decrease)
    VAR_NAMES <- colnames(X)
    
    
  
  
  ##############################
  ## Return results
  ##############################
  
  return(list(
    variable_importance = data.frame(
      VAR_NAMES = VAR_NAMES,
      mean_decrease_accuracy = mean_decrease_accuracy
    )
  ))
}