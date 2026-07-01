
# === Fully Automated Theoretical Mapper Partition ===
# Inputs:
#   n : total number of observations
#   a, b : constants controlling bias-variance tradeoff (use a=b=1 if unknown)
# Output:
#   l*, S*, q*, lambda*

optimal_mapper_auto <- function(n, alpha = .5) {

  # Optimal number of intervals
  # l_star <- ((8 * b * n) / a)^(1/5)
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




mapper_cv_function <- function(data, secondary = FALSE, max_cv = 10, seed = 12345,
                               max_cluster_no = 5) {

  # ============================================================
  # Storage objects across CV runs
  # ============================================================
  mapper_Prediction_CV <- list()              # predicted class for each CV split
  mapper_probability_prediction <- list()     # predicted probabilities for each CV split
  True_Y_CV <- list()                         # true test labels for each CV split

  opt_intervals_CV <- c()                     # selected number of intervals per CV
  accuracy_rate_CV <- list()                  # accuracy/QWK for each alpha per CV

  # ============================================================
  # Cross-validation loop
  # ============================================================
  for (cv in 1:max_cv) {

    cat("CV =", cv, "\n")

    # ------------------------------------------------------------
    # Randomly split data into training and test sets
    # 20% test, 80% train
    # ------------------------------------------------------------
    set.seed(seed + cv)
    kk <- sample(1:nrow(data), size = round(0.2 * nrow(data)))

    train_data <- data[-kk, ]
    test_data  <- data[kk, ]

    # X = predictors, Y = ordinal response
    X <- as.matrix(train_data[, -1])
    Y <- train_data[, 1]

    # ============================================================
    # Filter function: Ordinal LDA / ordASDA
    # This learns a 1D supervised projection for Mapper
    # ============================================================
    filter_olda_beta <- function(X, Y) {
      X <- as.matrix(X)
      Y <- ordered(Y)

      fit <- accSDA::ordASDA(
        Xt = X,
        Yt = as.numeric(Y)
      )

      # return the coefficient vector for the filter projection
      fit$beta[1:ncol(X), , drop = FALSE]
    }

    # Learned filter coefficients
    filter_beta <- filter_olda_beta(X, Y)

    # Project training points onto the 1D filter axis
    filter_values <- drop(X %*% filter_beta)
    sorted_values <- sort(filter_values)

    # ============================================================
    # Balanced cover construction
    # Builds overlapping intervals along the filter function
    # ============================================================
    create_intervals <- function(sorted_values, alpha = 0.5) {
      n <- length(sorted_values)

      # optimal_mapper_auto() chooses cover size parameters
      opt <- optimal_mapper_auto(n, alpha)
      l <- opt$l_star
      S <- floor(opt$S_star)

      # adaptive stride so the intervals cover all points
      q_adapt <- (n - S) / (l - 1)

      intervals_idx <- vector("list", l)

      for (j in 1:l) {
        start_idx <- 1 + floor((j - 1) * q_adapt)
        end_idx   <- min(start_idx + S - 1, n)
        intervals_idx[[j]] <- start_idx:end_idx
      }

      # convert index intervals into actual filter-value intervals
      data.frame(
        start = sapply(intervals_idx, function(idx) min(sorted_values[idx])),
        end   = sapply(intervals_idx, function(idx) max(sorted_values[idx]))
      )
    }

    # ============================================================
    # Alpha tuning loop
    # Different alpha values give different Mapper covers
    # ============================================================
    intervals_NUM <- c()
    mapper_Prediction_alp <- list()
    mapper_probability_prediction_alp <- list()
    accuracy_rate <- c()

    alpha_seq <- seq(0.9, 0.1, -0.1)

    for (alp in seq_along(alpha_seq)) {

      # ----------------------------------------------------------
      # Build intervals for this alpha
      # ----------------------------------------------------------
      intervals <- create_intervals(sorted_values, alpha = alpha_seq[alp])
      intervals_NUM[alp] <- nrow(intervals)

      # skip repeated numbers of intervals
      if (nrow(intervals) %in% intervals_NUM[-alp]) next

      # ==========================================================
      # Subset training data into overlapping interval bins
      # ==========================================================
      subsets <- vector("list", nrow(intervals))
      row_indices_sub <- vector("list", nrow(intervals))

      for (j in seq_len(nrow(intervals))) {

        # find which training rows fall into interval j
        idx_j <- which(filter_values >= intervals$start[j] &
                         filter_values <= intervals$end[j])

        row_indices_sub[[j]] <- idx_j
        subsets[[j]] <- X[idx_j, , drop = FALSE]
      }

      # ==========================================================
      # Cluster within each interval
      # This forms the Mapper nodes
      # ==========================================================
      clusters <- list()
      all_cluster_assignments <- list()

      for (j in seq_along(subsets)) {

        # if interval has 0 or 1 point, force one trivial cluster
        if (nrow(subsets[[j]]) <= 1) {
          cluster_assignments <- rep(1, nrow(subsets[[j]]))
          clusters[[j]] <- list(
            membership = cluster_assignments,
            csize = table(factor(cluster_assignments, levels = 1))
          )

          all_cluster_assignments[[j]] <- data.frame(
            row_id = row_indices_sub[[j]],
            cluster = cluster_assignments,
            interval = j
          )
          next
        }

        # pairwise Euclidean distances
        d_matrix <- parallelDist::parDist(subsets[[j]], method = "euclidean")

        # hierarchical clustering
        hc <- fastcluster::hclust(d_matrix, method = "ward.D2")

        # candidate number of clusters
        k_vals <- 2:max_cluster_no
        k_vals <- k_vals[k_vals < nrow(subsets[[j]])]

        if (length(k_vals) == 0) {
          optimal_k <- 1
        } else {
          # choose k using average silhouette width
          sil_vals <- sapply(k_vals, function(k) {
            cl <- cutree(hc, k)
            mean(cluster::silhouette(cl, d_matrix)[, 3])
          })
          optimal_k <- if (all(is.na(sil_vals))) 1 else k_vals[which.max(sil_vals)]
        }

        # final cluster labels in interval j
        cluster_assignments <- cutree(hc, k = optimal_k)

        clusters[[j]] <- list(
          membership = cluster_assignments,
          csize = table(cluster_assignments)
        )

        # keep row IDs so we can connect clusters across intervals
        all_cluster_assignments[[j]] <- data.frame(
          row_id = row_indices_sub[[j]],
          cluster = cluster_assignments,
          interval = j
        )
      }

      # ==========================================================
      # Build adjacency across neighboring intervals
      # Two clusters are connected if they share at least one point
      # due to overlap in the cover
      # ==========================================================
      cross_interval_edges <- data.frame(
        interval_from = integer(),
        cluster_from  = integer(),
        interval_to   = integer(),
        cluster_to    = integer()
      )

      for (j in seq_along(all_cluster_assignments)) {
        df_j <- all_cluster_assignments[[j]]

        # connect to previous interval
        if (j > 1) {
          df_prev <- all_cluster_assignments[[j - 1]]
          for (c1 in unique(df_j$cluster)) {
            rows_c1 <- df_j$row_id[df_j$cluster == c1]
            for (c2 in unique(df_prev$cluster)) {
              rows_c2 <- df_prev$row_id[df_prev$cluster == c2]
              if (length(intersect(rows_c1, rows_c2)) > 0) {
                cross_interval_edges <- rbind(
                  cross_interval_edges,
                  data.frame(
                    interval_from = j,
                    cluster_from  = c1,
                    interval_to   = j - 1,
                    cluster_to    = c2
                  )
                )
              }
            }
          }
        }

        # connect to next interval
        if (j < length(all_cluster_assignments)) {
          df_next <- all_cluster_assignments[[j + 1]]
          for (c1 in unique(df_j$cluster)) {
            rows_c1 <- df_j$row_id[df_j$cluster == c1]
            for (c2 in unique(df_next$cluster)) {
              rows_c2 <- df_next$row_id[df_next$cluster == c2]
              if (length(intersect(rows_c1, rows_c2)) > 0) {
                cross_interval_edges <- rbind(
                  cross_interval_edges,
                  data.frame(
                    interval_from = j,
                    cluster_from  = c1,
                    interval_to   = j + 1,
                    cluster_to    = c2
                  )
                )
              }
            }
          }
        }
      }

      # create text labels for nodes like I3_C2
      if (nrow(cross_interval_edges) > 0) {
        cross_interval_edges$node_from <- paste0(
          "I", cross_interval_edges$interval_from, "_C", cross_interval_edges$cluster_from
        )
        cross_interval_edges$node_to <- paste0(
          "I", cross_interval_edges$interval_to, "_C", cross_interval_edges$cluster_to
        )
      } else {
        cross_interval_edges$node_from <- character(0)
        cross_interval_edges$node_to <- character(0)
      }

      # ==========================================================
      # Predict test observations
      # ==========================================================
      new_point <- as.matrix(test_data[, -1])
      true_Y    <- test_data[, 1]

      # project test points using the same learned filter
      filter_value <- drop(new_point %*% filter_beta)

      # ----------------------------------------------------------
      # Find which interval(s) each test point belongs to
      # because intervals overlap, some points can fall in >1 interval
      # ----------------------------------------------------------
      interval_index <- vector("list", length(filter_value))

      for (i in seq_along(filter_value)) {
        if (filter_value[i] < min(intervals$start)) {
          interval_index[[i]] <- 1
        } else if (filter_value[i] > max(intervals$end)) {
          interval_index[[i]] <- nrow(intervals)
        } else {
          interval_index[[i]] <- which(
            filter_value[i] >= intervals$start &
              filter_value[i] <= intervals$end
          )
        }
      }

      weighted_average_prediction <- c()
      weighted_average_prediction_prob <- list()

      # ==========================================================
      # Loop over test points
      # ==========================================================
      for (i in seq_along(interval_index)) {

        new_point_pos <- new_point[i, ]
        interval_index_pos <- interval_index[[i]]

        closest_cluster_for_j_and_j1 <- c()

        # ========================================================
        # CASE 1: Test point falls into exactly one interval
        # ========================================================
        if (length(interval_index_pos) == 1) {

          idx_int <- interval_index_pos
          subset_new_point <- subsets[[idx_int]]

          if (nrow(subset_new_point) == 0) {
            weighted_average_prediction[i] <- NA
            weighted_average_prediction_prob[[i]] <- NA
            next
          }

          # ------------------------------------------------------
          # Find centroid of each cluster in this interval
          # ------------------------------------------------------
          cluster_ids <- unique(clusters[[idx_int]]$membership)

          cluster_centroids <- sapply(cluster_ids, function(cl_id) {
            cluster_points <- subset_new_point[
              clusters[[idx_int]]$membership == cl_id, , drop = FALSE
            ]
            if (nrow(cluster_points) == 1) cluster_points else colMeans(cluster_points)
          })

          if (is.null(dim(cluster_centroids))) {
            cluster_centroids <- matrix(cluster_centroids, ncol = 1)
          }

          # distance from test point to each cluster centroid
          distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
            sqrt(sum((centroid - new_point_pos)^2))
          })

          # avoid singleton clusters by giving them huge distance
          csize_vec <- as.numeric(clusters[[idx_int]]$csize)
          distances_to_centroids[which(csize_vec == 1)] <- 99999

          # choose the nearest cluster
          closest_cluster <- which.min(distances_to_centroids)

          # cluster points and corresponding row IDs
          cluster_points <- subset_new_point[
            clusters[[idx_int]]$membership == closest_cluster, , drop = FALSE
          ]

          row_indices <- row_indices_sub[[idx_int]][
            clusters[[idx_int]]$membership == closest_cluster
          ]

          response_values <- Y[row_indices]

          # ------------------------------------------------------
          # If secondary = TRUE, add neighboring connected nodes
          # from adjacent intervals
          # ------------------------------------------------------
          if (secondary && nrow(cross_interval_edges) > 0) {
            closest_node <- paste0("I", idx_int, "_C", closest_cluster)

            connected_rows <- subset(cross_interval_edges, node_from == closest_node)
            node_to_list <- connected_rows$node_to

            if (length(node_to_list) > 0) {
              all_nodes <- lapply(node_to_list, function(node) {
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster  <- as.numeric(parts[2])

                subsets[[interval]][clusters[[interval]]$membership == cluster, , drop = FALSE]
              })

              all_nodes_combined <- unique(do.call(rbind, all_nodes))

              # combine original closest cluster with connected neighbors
              cluster_points <- as.matrix(unique(rbind(all_nodes_combined, cluster_points)))

              # map back to original training rows
              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })

              row_indices <- row_indices[!is.na(row_indices)]
              response_values <- Y[row_indices]
            }
          }

        } else {

          # ======================================================
          # CASE 2: Test point falls into multiple overlapping intervals
          # ======================================================
          cluster_points <- NULL
          row_indices_list <- list()

          for (ind in seq_along(interval_index_pos)) {
            idx_int <- interval_index_pos[ind]
            subset_new_point <- subsets[[idx_int]]

            # find nearest cluster inside each candidate interval
            cluster_ids <- unique(clusters[[idx_int]]$membership)

            cluster_centroids <- sapply(cluster_ids, function(cl_id) {
              pts <- subset_new_point[
                clusters[[idx_int]]$membership == cl_id, , drop = FALSE
              ]
              if (nrow(pts) == 1) pts else colMeans(pts)
            })

            if (is.null(dim(cluster_centroids))) {
              cluster_centroids <- matrix(cluster_centroids, ncol = 1)
            }

            distances_to_centroids <- apply(cluster_centroids, 2, function(centroid) {
              sqrt(sum((centroid - new_point_pos)^2))
            })

            csize_vec <- as.numeric(clusters[[idx_int]]$csize)
            distances_to_centroids[which(csize_vec == 1)] <- 99999

            closest_cluster_for_j_and_j1[ind] <- which.min(distances_to_centroids)
          }

          # combine all nearest clusters from the overlapping intervals
          for (ii in seq_along(interval_index_pos)) {
            jj <- interval_index_pos[ii]
            cluster_id <- closest_cluster_for_j_and_j1[ii]

            tmp <- subsets[[jj]][clusters[[jj]]$membership == cluster_id, , drop = FALSE]

            if (is.null(cluster_points)) {
              cluster_points <- tmp
            } else {
              tmp_unique <- tmp[!apply(tmp, 1, function(row) {
                any(apply(cluster_points, 1, function(existing) all(existing == row)))
              }), , drop = FALSE]
              cluster_points <- rbind(cluster_points, tmp_unique)
            }

            row_indices_list[[ii]] <- row_indices_sub[[jj]][
              clusters[[jj]]$membership == cluster_id
            ]
          }

          row_indices <- unique(unlist(row_indices_list))
          response_values <- Y[row_indices]

          # ------------------------------------------------------
          # If secondary = TRUE, add graph-neighbor nodes
          # ------------------------------------------------------
          if (secondary && nrow(cross_interval_edges) > 0) {
            closest_node <- paste0("I", interval_index_pos, "_C", closest_cluster_for_j_and_j1)
            connected_rows <- subset(cross_interval_edges, node_from %in% closest_node)

            node_to_list <- connected_rows$node_to
            different_node <- unique(setdiff(node_to_list, closest_node))

            if (length(different_node) > 0) {
              all_nodes <- lapply(different_node, function(node) {
                parts <- strsplit(node, "_C")[[1]]
                interval <- as.numeric(sub("I", "", parts[1]))
                cluster  <- as.numeric(parts[2])

                subsets[[interval]][clusters[[interval]]$membership == cluster, , drop = FALSE]
              })

              all_nodes_combined <- unique(do.call(rbind, all_nodes))

              cluster_points <- as.matrix(unique(rbind(all_nodes_combined, cluster_points)))

              row_indices <- apply(cluster_points, 1, function(row) {
                match(TRUE, apply(X, 1, function(x_row) all(x_row == row)))
              })

              row_indices <- row_indices[!is.na(row_indices)]
              response_values <- Y[row_indices]
            }
          }
        }

        # ========================================================
        # Weighted ordinal prediction
        # Use inverse-distance weights from the local neighborhood
        # ========================================================
        combined_df <- data.frame(response = response_values, cluster_points)
        new_point_df <- as.data.frame(as.list(new_point_pos))

        # squared Euclidean distance from test point to each neighbor
        distances <- apply(combined_df[, -1, drop = FALSE], 1, function(x) {
          sum((x - unlist(new_point_df))^2)
        })

        neighbor_responses <- combined_df$response

        # inverse-distance weights
        weights <- 1 / (distances + 1e-6)
        weights <- weights / sum(weights)

        categories <- levels(combined_df$response)

        # cumulative ordinal probabilities
        cum_probs <- sapply(categories, function(c) {
          sum(weights[neighbor_responses <= c]) / sum(weights)
        })

        # category-specific probabilities
        cat_probs <- sapply(categories, function(c) {
          sum(weights[neighbor_responses == c]) / sum(weights)
        })

        # posterior median decision rule for ordinal outcome
        median_index <- which(cum_probs >= 0.5)[1]
        pred_class <- categories[median_index]

        weighted_average_prediction_prob[[i]] <- cat_probs
        weighted_average_prediction[i] <- as.numeric(as.character(pred_class))
      }

      # store results for this alpha
      mapper_Prediction_alp[[alp]] <- weighted_average_prediction
      mapper_probability_prediction_alp[[alp]] <- weighted_average_prediction_prob

      # ----------------------------------------------------------
      # Evaluate prediction using weighted kappa
      # ----------------------------------------------------------
      levels_Y <- as.numeric(levels(true_Y))
      pred_class_mapper <- factor(weighted_average_prediction, levels = levels_Y)
      y_numeric <- as.numeric(as.character(true_Y))

      accuracy_rate[alp] <- irr::kappa2(
        data.frame(y_numeric, as.numeric(as.character(pred_class_mapper))),
        weight = "squared"
      )$value
    }

    # ============================================================
    # Select the alpha with best performance in this CV split
    # ============================================================
    Position_accuracy_rate <- max(which(accuracy_rate == max(accuracy_rate, na.rm = TRUE)))

    mapper_Prediction_CV[[cv]] <- mapper_Prediction_alp[[Position_accuracy_rate]]
    mapper_probability_prediction[[cv]] <- mapper_probability_prediction_alp[[Position_accuracy_rate]]
    True_Y_CV[[cv]] <- true_Y
    opt_intervals_CV[cv] <- intervals_NUM[Position_accuracy_rate]
    accuracy_rate_CV[[cv]] <- accuracy_rate
  }

  # ==============================================================
  # Return all outputs
  # ==============================================================
  return(list(
    mapper_prediction_cv = mapper_Prediction_CV,
    mapper_probability_prediction = mapper_probability_prediction,
    true_y_cv = True_Y_CV,
    opt_intervals_cv = opt_intervals_CV,
    accuracy_rate_cv = accuracy_rate_CV
  ))
}




