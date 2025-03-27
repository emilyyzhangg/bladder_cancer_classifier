library(survival)
library(survminer)
library(pROC)
library(rms)
library(VIM)

# load files
bladder_data <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/UROMOL_TaLG.teachingcohort.rds")
knowles_data <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/knowles_matched_TaLG_final.rds")

# remove constant vars
constant_vars <- names(bladder_data)[sapply(bladder_data, function(x) length(unique(x)) <= 1)]
bladder_data <- bladder_data[, !names(bladder_data) %in% constant_vars]

# find intersection of features
common_features <- intersect(names(bladder_data), names(knowles_data))
bladder_data <- bladder_data[, common_features, drop=FALSE]
knowles_data <- knowles_data[, common_features, drop=FALSE]

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
bladder_data <- kNN(bladder_data, k = 5)
knowles_data <- kNN(knowles_data, k = 5)
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

# multivariate cox model
multi_cox_model <- coxph(as.formula(paste("surv_obj ~", paste(significant_vars, collapse = " + "))), data = bladder_data)
bladder_data$risk_score <- predict(multi_cox_model, type = "risk")


# validation
knowles_data$surv_obj <- Surv(time = knowles_data$RFS_time, event = knowles_data$Recurrence)
knowles_data$risk_score <- predict(cox_model_final, newdata = knowles_data, type = "risk")

# discrimination
c_index <- concordance.index(knowles_data$risk_score, knowles_data$RFS_time, knowles_data$Recurrence, method = "noether")

#AUC
roc_curve <- roc(knowles_data$Recurrence, knowles_data$risk_score)
auc_value <- auc(roc_curve)

# calibration
cal_plot <- val.surv(cox_model_final, times = c(3, 6, 12, 24) * 30.44)

# output summaries
print(summary(cox_model_final))
print(c_index)
print(auc_value)


# plot kaplain-meier curve
km_fit <- survfit(Surv(RFS_time, Recurrence) ~ cut(risk_score, breaks = quantile(risk_score, probs = c(0, 0.5, 1))), data = knowles)
ggsurvplot(km_fit, data = knowles_data, pval = TRUE, conf.int = TRUE)





