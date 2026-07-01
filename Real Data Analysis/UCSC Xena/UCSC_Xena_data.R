
path<-"C:\\Users\\ahsanm8\\Desktop\\Dr. Nitai_Final_code\\Github_Material_Mapper"
source(file.path(path, "Mapper_Prediction_function.R"))
source(file.path(path, "Compititive_models_function.R"))

if (!require("pacman")) install.packages("pacman"); pacman::p_load(TDA, ggplot2, plotly, FNN, cluster, matrixStats, dbscan, igraph, rgl, mappeR, grid, ks, tidyr, devtools, fastcluster, DescTools, pROC, MASS, fclust, umap, mclust, NbClust, proxy, boot, pls,
                                                                   dplyr, infotheo, sigclust, randomForest, irr, accSDA, brant, RColorBrewer, factoextra, nnet, ordinalForest, survival, parallelDist,readxl)


# Set your file path
file_path <- "C:/Users/ahsanm8/Desktop/Dr. Nitai_Final_code/Github_Material_Mapper/Real Data Analysis/UCSC Xena"

clinical_file <- file.path(file_path, "TCGA.GBMLGG.sampleMap_GBMLGG_clinicalMatrix")
expr_file     <- file.path(file_path, "HiSeqV2")

# Read as a tab-delimited text file
clinical_data <- read.delim(clinical_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)


clin_selected <- clinical_data %>%
  dplyr::select(
    sampleID,
    histological_type,
    neoplasm_histologic_grade,
    
    Age_at_diag = age_at_initial_pathologic_diagnosis,
    KPS=karnofsky_performance_score
  )


# 01=solid tumor
# Extract the 4th field (sample type) from the TCGA barcode
sample_type <- sapply(strsplit(clin_selected$sampleID, split = "-"), `[`, 4)

# Keep only rows where sample type is "01"
clin_selected <- clin_selected[sample_type == "01", ]


# Assign G4 to GBM types
clin_selected$grade3[grepl("GBM", clin_selected$histological_type, ignore.case = TRUE)] <- "G4"

clin_selected$grade3[clin_selected$neoplasm_histologic_grade == "G2"] <- "G2"

clin_selected$grade3[clin_selected$neoplasm_histologic_grade == "G3"] <- "G3"



# Convert to ordered factor
clin_selected$grade3 <- factor(clin_selected$grade3,levels = c("G2", "G3", "G4"), ordered = TRUE)








# Load gene expression data
expr <- read.table(expr_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


# Transpose and add sample ID column
expr_t <- as.data.frame(t(expr))
expr_t$sampleID <- rownames(expr_t)

# Merge with clinical data
merged_data <- inner_join(clin_selected, expr_t, by = "sampleID")


Age_at_diag<-merged_data$Age_at_diag
KPS<-merged_data$KPS

merged_data<-merged_data[, !(colnames(merged_data) %in% c("Age_at_diag","KPS"))]


rows_with_na <- which(!complete.cases(merged_data))

Age_at_diag<-Age_at_diag[-rows_with_na]
KPS<-KPS[-rows_with_na]

merged_data<-drop_na(merged_data)


Y<-merged_data$grade3


# Drop non-gene columns and keep only gene expression features
X <- merged_data[, !(colnames(merged_data) %in%
                       c("sampleID", "histological_type", "neoplasm_histologic_grade", "grade3","Age_at_diag","KPS"))]





# keep top 10% data with most variability
gene_var <- colVars(as.matrix(X))
keep <- gene_var > quantile(gene_var, 0.9)   
X_var <- X[, keep]
X<-X_var

data<-data.frame(cbind(Y,X))

data$Y <- factor(as.numeric(data$Y),ordered=T)

Y<-data[,1]
X<-data[,-1]




#####################################
#Variable selection by Random forest
#####################################
set.seed(123)

rf_model <- randomForest(
  x = X,
  y = Y,
  ntree = 500,
  importance = TRUE,
  proximity = FALSE
)

rf_importance <- randomForest::importance(rf_model, type = 1)   # type=1 = MeanDecreaseAccuracy

rf_imp_df <- data.frame(
  Variable = rownames(rf_importance),
  MeanDecreaseAccuracy = rf_importance[,1]
)

rf_imp_df <- rf_imp_df %>%
  arrange(desc(MeanDecreaseAccuracy))


selected_vars <- rf_imp_df %>%
  filter(MeanDecreaseAccuracy >= 3)

selected_genes <- selected_vars$Variable

data1<-data
data <- data1[, c("Y", selected_genes)]

Y<-data[,1]
X<-data[,-1]

save(data,file = file.path(file_path, "UCSC_Xena_data.RData"))
##########################






mapper_result <- mapper_cv_function(data = data, secondary = FALSE, max_cv = 10)


True_Y_CV<-mapper_result$true_y_cv
mapper_probability_prediction<-mapper_result$mapper_probability_prediction
mapper_Prediction_CV<-mapper_result$mapper_prediction_cv


results_competitive_models <- competitive_models_cv(data = data, max_cv = 10)


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












