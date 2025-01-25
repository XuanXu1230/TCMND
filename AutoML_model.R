#### Build Machine Learning Models ####
library(tidymodels)
library(h2o)

# Prepare data
data <- data.frame(t(exp_sub))
data$group <- as.factor(Group)  # Convert target column to factor type

# Initialize H2O
h2o.init(nthreads = 20)

# Initialize data frames to store evaluation results and variable importance
evaluation_results <- data.frame(
  Module = character(),
  Model_ID = character(),
  AUC = numeric(),
  Logloss = numeric(),
  AUCPR = numeric(),
  Mean_Per_Class_Error = numeric(),
  RMSE = numeric(),
  MCC = numeric(),  # Matthews Correlation Coefficient
  F1_Score = numeric(),  # F1 Score
  Training_Time = numeric(),
  stringsAsFactors = FALSE
)

variable_importance <- data.frame(
  Module = character(),
  Gene = character(),
  Importance = numeric(),
  stringsAsFactors = FALSE
)

# Custom stratified sampling function
stratified_split <- function(data, target_col, ratio = 0.7, seed = 123) {
  set.seed(seed)
  target_levels <- h2o.levels(data[[target_col]])
  
  train_indices <- NULL
  test_indices <- NULL
  
  for (level in target_levels) {
    level_data <- data[data[[target_col]] == level, ]
    split <- h2o.splitFrame(level_data, ratios = ratio, seed = seed)
    train_indices <- if (is.null(train_indices)) split[[1]] else h2o.rbind(train_indices, split[[1]])
    test_indices <- if (is.null(test_indices)) split[[2]] else h2o.rbind(test_indices, split[[2]])
  }
  
  list(train = train_indices, test = test_indices)
}

# Helper function to calculate additional metrics
calculate_metrics <- function(perf) {
  confusion_matrix <- h2o.confusionMatrix(perf)
  labels <- h2o.levels(test$group)
  
  tp <- as.numeric(confusion_matrix[labels[1], labels[1]])
  tn <- as.numeric(confusion_matrix[labels[2], labels[2]])
  fp <- as.numeric(confusion_matrix[labels[2], labels[1]])
  fn <- as.numeric(confusion_matrix[labels[1], labels[2]])
  
  mcc <- ifelse((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) > 0,
                ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)),
                NA)
  f1_score <- ifelse((2 * tp + fp + fn) > 0, 2 * tp / (2 * tp + fp + fn), NA)
  
  list(mcc = mcc, f1_score = f1_score)
}

# Main loop to train models for each module
for (module_name in names(mod_hub)) {
  genes <- mod_hub[[module_name]]
  valid_genes <- intersect(genes, colnames(data))
  
  if (length(valid_genes) < 3) {
    warning(sprintf("Skipping module %s: insufficient valid genes", module_name))
    next
  }
  
  module_data <- data[, c(valid_genes, "group"), drop = FALSE]
  module_data <- as.h2o(module_data)
  
  if (!is.h2o(module_data)) {
    stop("Conversion to H2OFrame failed!")
  }
  
  splits <- stratified_split(module_data, target_col = "group", ratio = 0.7, seed = 123)
  train <- splits$train
  test <- splits$test
  
  aml <- h2o.automl(
    x = valid_genes,
    y = "group",
    training_frame = train,
    max_models = 50,
    nfolds = 10,
    seed = 123
  )
  
  leaderboard <- as.data.frame(h2o.get_leaderboard(aml, extra_columns = "ALL"))
  
  for (i in 1:nrow(leaderboard)) {
    model_id <- leaderboard$model_id[i]
    model <- h2o.getModel(model_id)
    perf <- h2o.performance(model, newdata = test)
    
    metrics <- calculate_metrics(perf)
    
    evaluation_results <- rbind(evaluation_results, data.frame(
      Module = module_name,
      Model_ID = model_id,
      AUC = leaderboard$auc[i],
      Logloss = leaderboard$logloss[i],
      AUCPR = leaderboard$aucpr[i],
      Mean_Per_Class_Error = leaderboard$mean_per_class_error[i],
      RMSE = leaderboard$rmse[i],
      MCC = metrics$mcc,
      F1_Score = metrics$f1_score,
      Training_Time = leaderboard$training_time_ms[i] / 1000
    ))
  }
  
  top_models <- leaderboard$model_id[1:5]
  weights <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  
  for (j in seq_along(top_models)) {
    model_id <- top_models[j]
    model <- h2o.getModel(model_id)
    
    varimp <- h2o.varimp(model)
    if (is.null(varimp) || nrow(varimp) == 0) {
      warning(sprintf("Model %s does not have variable importances. Skipping.", model_id))
      next
    }
    
    varimp <- as.data.frame(varimp)
    varimp$Importance <- varimp$relative_importance * weights[j]
    varimp$Module <- module_name
    
    variable_importance <- rbind(variable_importance, varimp[, c("Module", "variable", "Importance")])
  }
}

# Summarize variable importance across models
variable_importance_summary <- aggregate(
  Importance ~ Module + variable,
  data = variable_importance,
  sum
)
colnames(variable_importance_summary) <- c("Module", "Gene", "Importance")

# Sort genes by importance within each module
variable_importance_summary <- variable_importance_summary[order(variable_importance_summary$Module, -variable_importance_summary$Importance), ]

# Sort evaluation results by AUCPR
evaluation_results <- evaluation_results[order(-evaluation_results$AUCPR), ]

# Save results
output_dir <- "~/your_directory/"  # Replace with your desired output directory
save(mod_hub, evaluation_results, variable_importance_summary, file = paste0(output_dir, num, "-", name, "_deg0.1_wgcna_AutoML50.Rdata"))
write.csv(evaluation_results, paste0(output_dir, num, "-", name, "_deg0.1_wgcna_module_model_results_with_metrics.csv"), row.names = FALSE)

# Shut down H2O
h2o.shutdown(prompt = FALSE)