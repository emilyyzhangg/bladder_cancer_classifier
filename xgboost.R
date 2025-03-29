library(survival)
library(rms)
library(pROC)
library(xgboost)
library(caret)

bladder_data <- readRDS("./data/UROMOL_TaLG.teachingcohort.rds")
knowles_data <- readRDS("./data/knowles_matched_TaLG_final.rds")

# remove constant vars
constant_vars <- names(bladder_data)[sapply(bladder_data, function(x) length(unique(x)) <= 1)]
bladder_data <- bladder_data[, !names(bladder_data) %in% constant_vars]

# also removing these
bladder_data <- bladder_data[, !(names(bladder_data) %in% c("Progression", "PFS_time"))]

# find intersection of features
common_features <- intersect(names(bladder_data), names(knowles_data))
bladder_data <- bladder_data[, common_features, drop=FALSE]
colnames(bladder_data) <- gsub("FOXP1-AS1", "FOXP1AS1", colnames(bladder_data))
knowles_data <- knowles_data[, common_features, drop=FALSE]
colnames(knowles_data) <- gsub("FOXP1-AS1", "FOXP1AS1", colnames(knowles_data))

# finding common genes in exprs matrices
bladder_exprs <- as.data.frame(bladder_data$exprs)
knowles_exprs <- as.data.frame(knowles_data$exprs)

common_genes <- intersect(colnames(bladder_exprs), colnames(knowles_exprs))
bladder_exprs <- bladder_exprs[, common_genes, drop=FALSE]
knowles_exprs <- knowles_exprs[, common_genes, drop=FALSE]

# combine clinical data and matrices
bladder_data <- cbind(bladder_data[, !names(bladder_data) %in% "exprs"], bladder_exprs )
knowles_data <- cbind(knowles_data[, !names(knowles_data) %in% "exprs"], knowles_exprs)

#KNN imputation
bladder_data <- kNN(bladder_data, k = 3, imp_var = FALSE, imp_suffix = FALSE)
knowles_data <- kNN(knowles_data, k = 3, imp_var = FALSE, imp_suffix = FALSE)

# converting chars to factor
char_cols <- sapply(bladder_data, is.character)
bladder_data[char_cols] <- lapply(bladder_data[char_cols], as.factor)
knowles_data[char_cols] <- lapply(knowles_data[char_cols], as.factor)

bladder_data <- subset(bladder_data, select = -`PFS_time.`)
knowles_data <- subset(knowles_data, select = -`PFS_time.`)

# univariate cox analysis to get significant vars
univariate_results <- data.frame(Variable = character(), p_value = numeric())
for (var in setdiff(names(bladder_data), c("RFS_time", "Recurrence"))) {
  single_cox_model <- coxph(as.formula(paste("surv_obj ~", paste0("`", var, "`"))), data = bladder_data)
  p_value <- summary(single_cox_model)$coefficients[, "Pr(>|z|)"]
  univariate_results <- rbind(univariate_results, data.frame(Variable = var, p_value = p_value))
}

# Select significant variables (feature selection)
significant_vars <- univariate_results$Variable[univariate_results$p_value < 0.05]

vars_to_keep <- c("RFS_time", "Recurrence", significant_vars)
bladder_data_filtered <- bladder_data[, vars_to_keep]
bladder_data_filtered <- subset(bladder_data_filtered, select = -Concomitant.CIS)
bladder_data_filtered <- subset(bladder_data_filtered, select = -BCG)
bladder_data_filtered <- subset(bladder_data_filtered, select = -FUtime_days.)
bladder_data_filtered <- subset(bladder_data_filtered, select = -PFS_time.)
bladder_data_filtered <- bladder_data_filtered %>%
  mutate(UROMOL2021.classification = recode(UROMOL2021.classification, "Class 2a" = "Class 2b"))

vars_to_keep_knowles <- intersect(vars_to_keep, names(knowles_data))
knowles_data_filtered <- knowles_data[, vars_to_keep_knowles]
knowles_data_filtered <- subset(knowles_data_filtered, select = -Concomitant.CIS)
knowles_data_filtered <- subset(knowles_data_filtered, select = -BCG)
knowles_data_filtered <- subset(knowles_data_filtered, select = -FUtime_days.)
knowles_data_filtered <- subset(knowles_data_filtered, select = -PFS_time.)
knowles_data_filtered <- knowles_data_filtered %>%
  mutate(UROMOL2021.classification = recode(UROMOL2021.classification, 
                                            "Class_1" = "Class 1",
                                            "Class_2b" = "Class 2b",
                                            "Class_3" = "Class 3"))
median_rfs_time <- median(bladder_data_filtered$RFS_time[bladder_data_filtered$Recurrence == 0], na.rm = TRUE)
# replace NA values in RFS_time in knowles_data_filtered
knowles_data_filtered <- knowles_data_filtered %>%
  mutate(RFS_time = ifelse(is.na(RFS_time), median_rfs_time, RFS_time))

predictors <- setdiff(names(bladder_data_filtered), c("RFS_time", "Recurrence"))

bladder_data_filtered <- bladder_data_filtered %>%
  mutate(across(all_of(predictors), ~ if(is.character(.)) as.factor(.) else .))

clean_column_names <- function(df) {
  colnames(df) <- gsub("`", "", colnames(df))  # Remove backticks
  colnames(df) <- gsub(" ", "", colnames(df))  # Remove spaces
  return(df)
}


# Prepare data for XGBoost
train_data <- model.matrix(~ . - 1, data = bladder_data_filtered[, predictors])
train_data <- clean_column_names(train_data)
train_labels <- bladder_data_filtered$RFS_time
train_events <- bladder_data_filtered$Recurrence

# Combine time and event as required for Cox model
dtrain <- xgb.DMatrix(data = train_data, label = train_labels, weight = train_events)

# Prepare Knowles dataset
complete_knowles <- knowles_data_filtered[complete.cases(knowles_data_filtered[, predictors]), ]
complete_knowles <- complete_knowles %>%
  mutate(across(all_of(predictors), ~ if(is.character(.)) as.factor(.) else .))

knowles_matrix <- model.matrix(~ . - 1, data = complete_knowles[, predictors])
knowles_matrix <- clean_column_names(knowles_matrix )

# Hyperparameter tuning grid
# xgb_grid <- expand.grid(
#   nrounds = 400,
#   eta = c(0.01, 0.1, 0.3),
#   max_depth = c(3, 6, 10),
#   gamma = c(0, 0.1, 0.2, 0.3),
#   colsample_bytree = c(0.3, 0.5, 0.7),
#   min_child_weight = c(1, 5),
#   subsample = c(0.5, 0.7, 0.8, 1)
# )
# 
# # Train model with caret
# train_control <- trainControl(method = "cv", number = 5, search = "grid", verboseIter = TRUE)
# xgb_model_caret <- train(
#   x = train_data,
#   y = train_labels,
#   weights = train_events,
#   method = "xgbTree",
#   trControl = train_control,
#   tuneGrid = xgb_grid,
#   objective = "survival:cox"
# )

# Best hyperparameters from the grid search
#print(xgb_model_caret$bestTune)
#best_xgb_model <- xgb_model_caret$finalModel

params <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik",
  max_depth = 3,
  eta = 0.01,
  gamma = 0.3,
  colsample_bytree = 0.3,
  min_child_weight = 5,
  subsample = 0.5
)

# best model
nrounds <- 400
best_xgb_model <- xgb.train(params = params, data = dtrain, nrounds = nrounds)

# risk scores (negative log hazard)
risk_scores <- predict(best_xgb_model, newdata = knowles_matrix)
complete_knowles$risk_score <- - risk_scores

# C-index
sink("./model_evaluation.txt")

c_index <- concordance(Surv(complete_knowles$RFS_time, complete_knowles$Recurrence) ~ complete_knowles$risk_score)
print(paste("C-index on Knowles data:", c_index$concordance))

# AUC
roc_curve <- roc(complete_knowles$Recurrence, complete_knowles$risk_score)
auc_value <- auc(roc_curve)
print(paste("AUC on Knowles data:", auc_value))

sink()

# Risk stratification into groups (Low, High)
complete_knowles$risk_group <- cut(
  complete_knowles$risk_score, 
  breaks = quantile(complete_knowles$risk_score, probs = c(0, 0.5, 1), na.rm = TRUE), 
  labels = c("Low", "High")
)

# Kaplan-Meier survival analysis
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = complete_knowles)
ggsurvplot(km_fit, data = complete_knowles, pval = TRUE, conf.int = FALSE,
           title = "Kaplan-Meier Curves by Risk Group (Knowles Data)")

xgb.save(best_xgb_model, "./best_xgb_model.model")
