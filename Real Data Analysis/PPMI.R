# Data is available here: https://www.ppmi-info.org

path<-"C:\\Users\\ahsanm8\\Desktop\\Dr. Nitai_Final_code\\Github_Material_Mapper"
source(file.path(path, "Mapper_Prediction_function.R"))
source(file.path(path, "Compititive_models_function.R"))
source(file.path(path, "Variable_Importance_Rank_Function.R"))

if (!require("pacman")) install.packages("pacman"); pacman::p_load(TDA, ggplot2, plotly, FNN, cluster, matrixStats, dbscan, igraph, rgl, mappeR, grid, ks, tidyr, devtools, fastcluster,
                                                                   DescTools, pROC, MASS, fclust, umap, mclust, NbClust, proxy, boot, pls, dplyr, infotheo, sigclust, randomForest, irr, 
                                                                   accSDA, brant, RColorBrewer, factoextra, nnet, ordinalForest, survival, parallelDist,readxl)


Curated_data<-read_excel("C:\\Users\\ahsanm8\\Desktop\\PPMI_Data\\PPMI_Curated_Data_Cut_Public_20250321.xlsx",
                         sheet = "20250310")

Curated_data_BL<-Curated_data[Curated_data$EVENT_ID=="BL",]

Curated_data_BL$PATNO<-as.numeric(Curated_data_BL$PATNO)


data<-Curated_data_BL%>%
  dplyr::select(PATNO, Y = hy,
                moca, upsit, quip, ess, gds, scopa,lns,
                abeta, tau, ptau, asyn,
                hvlt_discrimination,hvlt_immediaterecall,hvlt_retention,HVLTFPRL,HVLTRDLY,HVLTREC,
                VLTANIM,SDMTOTAL,bjlot,stai,rem,
                updrs1_score,updrs2_score,
                age,EDUCYRS,APOE_e4
  )

# Recode Y so that 2 and 3 are collapsed into 2
data$Y <- ifelse(data$Y == 3, 2, data$Y)


data<-data%>%
  dplyr::select(-PATNO)

data<-data.frame(data)


data<-drop_na(data)

data$Y <- factor(data$Y, levels = c(0,1,2), ordered = TRUE)



X<-data[,-1]
Y<-data[,1]


dim(data)


# #To get the importance score using MAPLE
# Important_Score<-MAPLE_Importance_score(data = data, secondary = FALSE, max_cv = 2)
# 
# variable_importance<-Important_Score$variable_importance
# VAR_NAMES<-variable_importance$VAR_NAMES
# 
# median_decrease_accuracy<-variable_importance$mean_decrease_accuracy
# INDX<-order(median_decrease_accuracy,decreasing = T)
# (median_decrease_accuracy_ordered <- round(median_decrease_accuracy[INDX],3))
# (Var_name_ordered<-VAR_NAMES[INDX])
# 
# data1<-data
# data <- data1[, c("Y", Var_name_ordered[1:4])] #based on the score 4 variables were exceded importance 5



# These are the variables selected by MAPLE
data <- data[, c("Y", "updrs2_score", "upsit", "rem", "stai")]
X<-data[,-1]
Y<-data[,1]







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












