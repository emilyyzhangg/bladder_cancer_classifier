library(survival)
library(randomForestSRC)
library(rms)
library(pROC)

bladder_data <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/UROMOL_TaLG.teachingcohort.rds")
knowles_data <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/knowles_matched_TaLG_final.rds")

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
bladder_data_filtered <- subset(bladder_data_filtered, select = -UROMOL2021.classification.1)
bladder_data_filtered <- subset(bladder_data_filtered, select = -Concomitant.CIS)
bladder_data_filtered <- subset(bladder_data_filtered, select = -BCG)
bladder_data_filtered <- subset(bladder_data_filtered, select = -FUtime_days.)
bladder_data_filtered <- bladder_data_filtered %>%
  mutate(UROMOL2021.classification = recode(UROMOL2021.classification, "Class 2a" = "Class 2b"))

vars_to_keep_knowles <- intersect(vars_to_keep, names(knowles_data))
knowles_data_filtered <- knowles_data[, vars_to_keep_knowles]
knowles_data_filtered <- subset(knowles_data_filtered, select = -Concomitant.CIS)
knowles_data_filtered <- subset(knowles_data_filtered, select = -BCG)
knowles_data_filtered <- subset(knowles_data_filtered, select = -FUtime_days.)
knowles_data_filtered <- knowles_data_filtered %>%
  mutate(UROMOL2021.classification = recode(UROMOL2021.classification, 
                                            "Class_1" = "Class 1",
                                            "Class_2b" = "Class 2b",
                                            "Class_3" = "Class 3"))

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

rf_model <- rfsrc(Surv(RFS_time, Recurrence) ~ ., data = bladder_data_filtered[, c("RFS_time", "Recurrence", predictors)], ntree = 400, importance = "permute", nimpute = 1)

print(rf_model)

complete_knowles <- knowles_data_filtered[complete.cases(knowles_data_filtered[, predictors]), ]
chf_pred <- predict(rf_model, newdata = complete_knowles[, predictors])$chf
complete_knowles$risk_score <- -1 * chf_pred[, ncol(chf_pred)]

c_index <- concordance(Surv(complete_knowles$RFS_time, complete_knowles$Recurrence) ~ complete_knowles$risk_score)
print(paste("C-index on Knowles data:", c_index$concordance))

roc_curve <- roc(complete_knowles$Recurrence, complete_knowles$risk_score)
auc_value <- auc(roc_curve)
print(paste("AUC on Knowles data:", auc_value))

complete_knowles$risk_group <- cut(complete_knowles$risk_score, breaks = quantile(complete_knowles$risk_score, probs = c(0, 0.33, 0.67, 1)), labels = c("Low", "Medium", "High"))
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = complete_knowles)
ggsurvplot(km_fit, data = complete_knowles, pval = TRUE, conf.int = FALSE,
           title = "Kaplan-Meier Curves by Risk Group (Knowles Data)")