library(survival)
library(survminer)
library(pROC)
library(rms)
library(VIM)
library(glmnet)
library(dplyr)

# load files
bladder_datan <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/UROMOL_TaLG.teachingcohort.rds")
knowles_datan <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/knowles_matched_TaLG_final.rds")

# remove constant vars
constant_varsn <- names(bladder_data)[sapply(bladder_data, function(x) length(unique(x)) <= 1)]
bladder_datan <- bladder_data[, !names(bladder_data) %in% constant_vars]

# also removing these
bladder_datan <- bladder_data[, !(names(bladder_data) %in% c("Progression", "PFS_time"))]

bladder_datan <- bladder_data[!(is.na(bladder_data$Recurrence) & is.na(bladder_data$RFS_time)), ]
knowles_datan <- knowles_data[!(is.na(knowles_data$Recurrence) & is.na(knowles_data$RFS_time)), ]


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

# removing rows where target feature is not available
bladder_data <- bladder_data[!(is.na(bladder_data$Recurrence) & is.na(bladder_data$RFS_time)), ]
knowles_data <- knowles_data[!(is.na(knowles_data$Recurrence) & is.na(knowles_data$RFS_time)), ]

#KNN imputation
bladder_data <- kNN(bladder_data, k = 3)
knowles_data <- kNN(knowles_data, k = 3)

# survival variables
surv_obj <- Surv(time = bladder_data$RFS_time, event = bladder_data$Recurrence)

# univariate cox analysis to get significant vars
univariate_results <- data.frame(Variable = character(), p_value = numeric())
for (var in setdiff(names(bladder_data), c("RFS_time", "Recurrence"))) {
  single_cox_model <- coxph(as.formula(paste("surv_obj ~", paste0("`", var, "`"))), data = bladder_data)
  p_value <- summary(single_cox_model)$coefficients[, "Pr(>|z|)"]
  univariate_results <- rbind(univariate_results, data.frame(Variable = var, p_value = p_value))
}

# selecting significant variables (feature selection)
significant_vars <- univariate_results$Variable[univariate_results$p_value < 0.05]
significant_vars <- intersect(significant_vars, colnames(bladder_data))
#significant_vars <- gsub("FOXP1-AS1", "FOXP1AS1", significant_vars)

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
max_rfs_time <- max(bladder_data_filtered$RFS_time, na.rm = TRUE)
# replace NA values in RFS_time in knowles_data_filtered
knowles_data_filtered <- knowles_data_filtered %>%
  mutate(RFS_time = ifelse(is.na(RFS_time), max_rfs_time, RFS_time))

predictors <- setdiff(names(bladder_data_filtered), c("RFS_time", "Recurrence"))

bladder_data_filtered <- bladder_data_filtered %>%
  mutate(across(all_of(predictors), ~ if(is.character(.)) as.factor(.) else .))

complete_knowles <- knowles_data_filtered[complete.cases(knowles_data_filtered[, predictors]), ]
complete_knowles <- complete_knowles %>%
  mutate(across(all_of(predictors), ~ if(is.character(.)) as.factor(.) else .))

X <- as.matrix(bladder_data_filtered[, -c(Recurrence, RFS_time)])
y <- Surv(time = bladder_data_filtered$RFS_time, event = bladder_data_filtered$Recurrence)

# Fit lasso-Cox model w/ cv
cv_model <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)  # alpha=1 for LASSO

# predict risk scores
bladder_data_filtered$risk_score <- -predict(cv_model, newx = X, s = "lambda.min", type = "link")

cox_test <- coxph(Surv(RFS_time, Recurrence) ~ ., data = bladder_data_filtered[, -1])
ph_test <- cox.zph(cox_test)
print(ph_test)

c_index_train <- concordance(Surv(bladder_data_filtered$RFS_time,bladder_data_filtered$Recurrence) ~ bladder_data_filtered$risk_score)
print(c_index_train)$concordance

X_valid <- as.matrix(complete_knowles[, -c(Recurrence, RFS_time)])
complete_knowles$risk_score <- predict(cv_model, newx = X_valid, s = "lambda.min", type = "link")


# discrimination
c_index <- concordance(Surv(complete_knowles$RFS_time, complete_knowles$Recurrence) ~ complete_knowles$risk_score)
print(c_index)$concordance
#AUC
roc_curve <- roc(complete_knowles$Recurrence, complete_knowles$risk_score)
auc_value <- auc(roc_curve)

# calibration
cal_plot <- val.surv(multi_cox_model, times = c(3, 6))
plot(cal_plot, xlab = "Predicted Probability", ylab = "Observed Probability", 
     main = "Calibration Plot", col = "blue", pch = 16)
abline(0, 1, lty = 2, col = "red")

# output summaries
print(summary(multi_cox_model))
print(c_index)
print(auc_value)


# plot kaplain-meier curve
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ cut(risk_score, breaks = quantile(risk_score, probs = c(0, 0.5, 1))), data = complete_knowles)
ggsurvplot(km_fit, data = complete_knowles, pval = TRUE, conf.int = TRUE)





