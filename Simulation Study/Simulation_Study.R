# Install pacman if not already installed
if (!require(pacman)) install.packages("pacman")

# Load all required packages
pacman::p_load(
  TDA, ggplot2, plotly, FNN, cluster, matrixStats, dbscan, igraph,
  rgl, mappeR, grid, ks, tidyr, devtools, fastcluster,
  DescTools, pROC, MASS, fclust, umap, mclust, NbClust,
  proxy, boot, pls, dplyr, infotheo, sigclust,
  randomForest, irr, accSDA, brant, RColorBrewer,
  factoextra, nnet, ordinalForest, survival, parallelDist
)

RNGkind(sample.kind = "Rejection")




# path <- " "  #please set the path accordingly
# setwd(path)
#
path<-"C:\\Users\\ahsanm8\\Desktop\\Dr. Nitai_Final_code\\Github_Material_Mapper"
source(file.path(path, "Mapper_Prediction_function.R"))
source(file.path(path, "Compititive_models_function.R"))













#
simulate_ordinal_data <- function(n = 250, noise = 0.2, scenario = 1, seed = 123) {

  set.seed(seed)

  ## ---------------------------
  ## Correlated predictors
  ## ---------------------------
  rho <- 0.6
  mean_vec <- c(0, 0)
  sd_vec   <- c(0.5, 0.3)

  Sigma <- matrix(c(sd_vec[1]^2, rho*sd_vec[1]*sd_vec[2],
                    rho*sd_vec[1]*sd_vec[2], sd_vec[2]^2), nrow = 2)

  normals <- MASS::mvrnorm(n, mu = mean_vec, Sigma = Sigma)
  x1 <- exp(normals[,1])
  x2 <- exp(normals[,2])

  ## ---------------------------
  ## Branch assignment
  ## ---------------------------
  branch <- sample(c(1, 2, 3), size = n, replace = TRUE,
                   prob = c(0.30, 0.40, 0.30))

  ## ---------------------------
  ## Nonlinear x3
  ## ---------------------------
  x3 <- dplyr::case_when(
    branch == 1 ~ cos(1.2 * x1) + 0.6 * log1p(x2)^2,
    branch == 2 ~ cos(x1 * 0.6) * sqrt(x2) + 0.3 * x2,
    branch == 3 ~ tanh(x1 * 1.2) + 1.4 * sin(x2)
  ) + rnorm(n, 0, 0.1)

  ## ---------------------------
  ## Scenario logic
  ## ---------------------------
  eta <- numeric(n)
  Y <- character(n)

  if (scenario == 1) {
    ## ---------------------------
    ## Scenario 1: Branch-specific ordinal structure
    ## ---------------------------

    for (b in 1:3) {
      idx <- which(branch == b)

      if (length(idx) > 0) {

        if (b == 1) {
          eta[idx] <- 1.5 * x1[idx] + 0.04 * x2[idx] + 0.6 * x3[idx] +
            rnorm(length(idx), 0, noise)

        } else if (b == 2) {
          eta[idx] <- 0.04 * x1[idx] + 0.4 * x2[idx] + 1.5 * x3[idx] +
            rnorm(length(idx), 0, noise)

        } else if (b == 3) {
          eta[idx] <- 1 * x1[idx] + 0.8 * x2[idx] +
            1.0 * log1p(x3[idx]) +
            rnorm(length(idx), 0, noise)
        }

        cuts <- quantile(eta[idx],
                         probs = c(0.25, 0.45, 0.75, 0.9),
                         na.rm = TRUE)

        Y[idx] <- as.character(cut(
          eta[idx],
          breaks = c(-Inf, cuts, Inf),
          labels = paste0("Stage ", 1:5),
          include.lowest = TRUE,
          ordered_result = TRUE
        ))
      }
    }

  } else if (scenario == 2) {
    ## ---------------------------
    ## Scenario 2: Global ordinal structure
    ## ---------------------------

    eta <- 0.6 * x1 + 0.8 * x2 + 0.6 * sin(x3) +
      rnorm(n, 0, noise)

    cuts <- quantile(eta,
                     probs = c(0.25, 0.45, 0.75, 0.9),
                     na.rm = TRUE)

    Y <- as.character(cut(
      eta,
      breaks = c(-Inf, cuts, Inf),
      labels = paste0("Stage ", 1:5),
      include.lowest = TRUE,
      ordered_result = TRUE
    ))

  } else {
    stop("scenario must be either 1 or 2")
  }

  ## ---------------------------
  ## Final formatting
  ## ---------------------------
  Y <- factor(Y, ordered = TRUE)

  df <- data.frame(
    x1 = x1,
    x2 = x2,
    x3 = x3,
    branch = factor(branch),
    response = Y
  ) %>%
    dplyr::filter(!is.na(response))

  X <- df[, 1:3]

  Y_ord <- as.numeric(df$response) - 1
  Y_final <- factor(Y_ord, ordered = TRUE)

  data <- cbind(Y = Y_final, X)

  return(data)
}





## --------------------------------
## user choice of n,noise,scenario

data<-simulate_ordinal_data(n=250, noise = 0.2, scenario = 1)

## --------------------------------




mapper_result <- mapper_cv_function(data = data, secondary = FALSE, max_cv = 10)

names(mapper_result)

True_Y_CV<-mapper_result$true_y_cv
mapper_probability_prediction<-mapper_result$mapper_probability_prediction
mapper_Prediction_CV<-mapper_result$mapper_prediction_cv


results_competitive_models <- competitive_models_cv(data = data, max_cv = 10)

names(results_competitive_models)

multi_logistic_class_CV<-results_competitive_models$multinomial$class
ordinal_logistic_class_CV<-results_competitive_models$ordinal_logistic$class
rf_class_CV<-results_competitive_models$random_forest$class
ordfor_class_CV<-results_competitive_models$ordinal_forest$class









all_mat_list <- lapply(mapper_probability_prediction, function(x) {
  do.call(rbind, lapply(x, as.numeric))
})

all_mat_list <- lapply(all_mat_list, function(mat) {
  colnames(mat) <- levels(data[,1])
  mat
})


# If column names are category labels:
pred_class_mapper_nominal_list <- lapply(all_mat_list, function(mat) {
  apply(mat, 1, function(x) {
    as.numeric(colnames(mat)[which.max(x)])
  })
})







# True response as numeric
True_Y_numeric <- lapply(True_Y_CV, function(y) {
  as.numeric(as.character(y))
})



qwk_O_mapper_values <- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      as.numeric(as.character(mapper_Prediction_CV[[i]]))
    ),
    weight = "squared"
  )$value
})

mean_qwk_O_mapper <- mean(qwk_O_mapper_values)
se_qwk_O_mapper   <- sd(qwk_O_mapper_values) /
  sqrt(length(qwk_O_mapper_values))


qwk_mapper_nominal_values <- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      pred_class_mapper_nominal_list[[i]]
    ),
    weight = "squared"
  )$value
})

mean_qwk_mapper_nominal <- mean(qwk_mapper_nominal_values)
se_qwk_mapper_nominal <- sd(qwk_mapper_nominal_values) /
  sqrt(length(qwk_mapper_nominal_values))






O_Mapper_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- mapper_Prediction_CV[[i]]

  concordance(y_true ~ y_pred)$concordance
})

O_Mapper_c_index_mean <- mean(O_Mapper_c_index)
O_Mapper_c_index_se <- sd(O_Mapper_c_index) / sqrt(length(O_Mapper_c_index))


M_Mapper_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- pred_class_mapper_nominal_list[[i]]

  concordance(y_true ~ y_pred)$concordance
})

M_Mapper_c_index_mean <- mean(M_Mapper_c_index)
M_Mapper_c_index_se <- sd(M_Mapper_c_index) / sqrt(length(M_Mapper_c_index))











###################

qwk_multi_logistic <- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      as.numeric(as.character(multi_logistic_class_CV[[i]]))
    ),
    weight = "squared"
  )$value
})

mean_qwk_multi_logistic <- mean(qwk_multi_logistic)
se_qwk_multi_logistic   <- sd(qwk_multi_logistic) /
  sqrt(length(qwk_multi_logistic))


qwk_ordinal_logistic <- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      ordinal_logistic_class_CV[[i]]
    ),
    weight = "squared"
  )$value
})

mean_qwk_ordinal_logistic <- mean(qwk_ordinal_logistic)
se_qwk_ordinal_logistic <- sd(qwk_ordinal_logistic) /
  sqrt(length(qwk_ordinal_logistic))






qwk_RF<- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      as.numeric(as.character(rf_class_CV[[i]]))
    ),
    weight = "squared"
  )$value
})

mean_qwk_RF <- mean(qwk_RF)
se_qwk_RF   <- sd(qwk_RF) /
  sqrt(length(qwk_RF))



qwk_O_RF<- sapply(seq_along(True_Y_numeric), function(i) {
  kappa2(
    data.frame(
      True_Y_numeric[[i]],
      as.numeric(as.character(ordfor_class_CV[[i]]))
    ),
    weight = "squared"
  )$value
})

mean_qwk_O_RF <- mean(qwk_O_RF)
se_qwk_O_RF   <- sd(qwk_O_RF) /
  sqrt(length(qwk_O_RF))









MLR_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- as.numeric(as.character(multi_logistic_class_CV[[i]]))

  concordance(y_true ~ y_pred)$concordance
})

MLR_c_index_mean <- mean(MLR_c_index)
MLR_c_index_se <- sd(MLR_c_index) / sqrt(length(MLR_c_index))


OLR_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- as.numeric(as.character(ordinal_logistic_class_CV[[i]]))

  concordance(y_true ~ y_pred)$concordance
})

OLR_c_index_mean <- mean(OLR_c_index)
OLR_c_index_se <- sd(OLR_c_index) / sqrt(length(OLR_c_index))



RF_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- as.numeric(as.character(rf_class_CV[[i]]))

  concordance(y_true ~ y_pred)$concordance
})

RF_c_index_mean <- mean(RF_c_index)
RF_c_index_se <- sd(RF_c_index) / sqrt(length(RF_c_index))


ORF_c_index <- sapply(seq_along(True_Y_numeric), function(i) {

  y_true <- True_Y_numeric[[i]]
  y_pred <- as.numeric(as.character(ordfor_class_CV[[i]]))

  concordance(y_true ~ y_pred)$concordance
})

ORF_c_index_mean <- mean(ORF_c_index)
ORF_c_index_se <- sd(ORF_c_index) / sqrt(length(ORF_c_index))











results_table <- data.frame(
  Method = c(
    "Mapper (Ordinal)",
    "Mapper (Nominal)",
    "Multinomial Logistic",
    "Ordinal Logistic",
    "Random Forest",
    "Ordinal Random Forest"
  ),


  QWK = sprintf("%.3f (%.3f)",
                c(mean_qwk_O_mapper,
                  mean_qwk_mapper_nominal,
                  mean_qwk_multi_logistic,
                  mean_qwk_ordinal_logistic,
                  mean_qwk_RF,
                  mean_qwk_O_RF),

                c(se_qwk_O_mapper,
                  se_qwk_mapper_nominal,
                  se_qwk_multi_logistic,
                  se_qwk_ordinal_logistic,
                  se_qwk_RF,
                  se_qwk_O_RF)),


  C_index = sprintf("%.3f (%.3f)",
                    c(O_Mapper_c_index_mean,
                      M_Mapper_c_index_mean,
                      MLR_c_index_mean,
                      OLR_c_index_mean,
                      RF_c_index_mean,
                      ORF_c_index_mean),

                    c(O_Mapper_c_index_se,
                      M_Mapper_c_index_se,
                      MLR_c_index_se,
                      OLR_c_index_se,
                      RF_c_index_se,
                      ORF_c_index_se))
)


print(results_table)




