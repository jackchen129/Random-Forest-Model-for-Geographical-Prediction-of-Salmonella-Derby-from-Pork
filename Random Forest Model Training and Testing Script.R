# Geographical Classification (China vs. France vs. Mexico vs. United States) of
# Pork-Associated Salmonella Derby using a Random Forest Model
# ================================================================
# 01_main_workflow.R
# ================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, caret, ranger, readxl, here, janitor,
  ggplot2, pROC, RColorBrewer
)

set.seed(123)

# --- Global settings --------------------------------------------
CLASS_LEVELS <- c("China", "France", "Mexico", "USA")

output_dir <- here("results")
dir.create(output_dir, showWarnings = FALSE)

metadata_file <- "metadata.xlsx"
gene_file <- "gene_presence_absence.csv"

# --- Source scripts ---------------------------------------------
source(here("scripts/02_data_preprocessing.R"))
source(here("scripts/03_model_training.R"))
source(here("scripts/04_model_evaluation.R"))

# --- Workflow ---------------------------------------------------
message("Running geographical classification (China vs. France vs. Mexico vs. USA)...")

data <- load_data(metadata_file, gene_file, CLASS_LEVELS)
processed_data <- preprocess_data(data)

train_idx <- createDataPartition(processed_data$country, p = 0.61, list = FALSE)
train_data <- processed_data[train_idx, ]
test_data <- processed_data[-train_idx, ]

model <- train_model(train_data)
results <- evaluate_model(model, test_data, CLASS_LEVELS)

# --- Save outputs -----------------------------------------------
write_csv(results$predictions, file.path(output_dir, "predictions.csv"))
write_csv(results$variable_importance, file.path(output_dir, "feature_importance.csv"))
write_csv(as.data.frame(results$confusion_matrix$table),
          file.path(output_dir, "confusion_matrix.csv"))

saveRDS(model, file.path(output_dir, "rf_model.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))

sink(file.path(output_dir, "performance_summary.txt"))
cat("Geographical Classification: China vs. France vs. Mexico vs. USA\n")
cat("====================================\n\n")
cat("Accuracy:", round(results$confusion_matrix$overall["Accuracy"], 3), "\n")
cat("Kappa:", round(results$confusion_matrix$overall["Kappa"], 3), "\n")
cat("AUC:", round(results$auc, 3), "\n\n")
cat("Confusion Matrix:\n")
print(results$confusion_matrix$table)
sink()

message("✔ Analysis complete: ", normalizePath(output_dir))

# ================================================================
#  02_data_preprocessing.R
# ================================================================

load_data <- function(metadata_file, gene_file, CLASS_LEVELS) {
  tryCatch({
    stopifnot(file.exists(metadata_file))
    stopifnot(file.exists(gene_file))
    
    message("Loading metadata...")
    meta <- read_excel(metadata_file) %>%
      clean_names() %>%
      rename(Isolate = isolate) %>%
      mutate(country = case_when(
        country %in% c("China") ~ "China",
        country %in% c("France") ~ "France",
        country %in% c("Mexico") ~ "Mexico",
        country %in% c("USA", "United States") ~ "USA",
        TRUE ~ NA_character_
      )) %>%
      drop_na(country) %>%
      mutate(country = factor(country, levels = CLASS_LEVELS))
    
    message("Country distribution:")
    print(count(meta, country))
    
    message("Loading gene matrix...")
    acc_raw <- read_csv(gene_file, col_types = cols(.default = "c"))
    
    gene_names <- acc_raw$Gene
    
    acc_matrix <- acc_raw %>%
      select(-Gene) %>%
      mutate(across(everything(), ~ ifelse(. == "1", 1, 0))) %>%
      as.matrix()
    
    rownames(acc_matrix) <- gene_names
    
    acc_matrix_t <- as.data.frame(t(acc_matrix))
    colnames(acc_matrix_t) <- make_clean_names(gene_names)
    acc_matrix_t <- acc_matrix_t %>%
      rownames_to_column("Isolate")
    
    message("Merging...")
    combined <- inner_join(meta, acc_matrix_t, by = "Isolate")
    
    if (nrow(combined) == 0) stop("No matched isolates")
    
    message("Final sample size: ", nrow(combined))
    return(combined)
    
  }, error = function(e) {
    message("Error: ", e$message)
    return(NULL)
  })
}

preprocess_data <- function(data) {
  tryCatch({
    gene_cols <- setdiff(names(data), c("Isolate", "country"))
    
    nzv <- nearZeroVar(select(data, all_of(gene_cols)), names = TRUE)
    
    if (length(nzv) > 0) {
      message("Removing ", length(nzv), " near-zero variance genes")
      data <- select(data, -all_of(nzv))
    }
    
    data <- data[, colSums(is.na(data)) < nrow(data)]
    
    return(data)
    
  }, error = function(e) {
    message("Error: ", e$message)
    return(NULL)
  })
}

# ================================================================
# 03_model_training.R
# ================================================================

train_model <- function(train_data) {
  tryCatch({
    
    ctrl <- trainControl(
      method = "cv",
      number = 5,
      classProbs = TRUE,
      summaryFunction = multiClassSummary,
      savePredictions = "final"
    )
    
    model <- train(
      x = select(train_data, -Isolate, -country),
      y = train_data$country,
      method = "ranger",
      importance = "permutation",
      num.trees = 750,
      metric = "Accuracy",
      trControl = ctrl
    )
    
    return(model)
    
  }, error = function(e) {
    message("Training error: ", e$message)
    return(NULL)
  })
}

# ================================================================
# 04_model_evaluation.R
# ================================================================

evaluate_model <- function(model, test_data, CLASS_LEVELS) {
  tryCatch({
    
    preds_prob <- predict(model, test_data, type = "prob")
    preds_class <- predict(model, test_data)
    
    preds <- bind_cols(
      SampleID = test_data$Isolate,
      actual = test_data$country,
      preds_prob,
      predicted = preds_class
    )
    
    cm <- confusionMatrix(
      preds$predicted,
      preds$actual
    )
    
    var_imp <- varImp(model)$importance %>%
      rownames_to_column("Feature") %>%
      arrange(desc(Overall))
    
    roc_obj <- multiclass.roc(
      response = preds$actual,
      predictor = as.matrix(preds_prob),
      levels = CLASS_LEVELS
    )
    
    # --- Plots ---------------------------------------------------
    imp_plot <- var_imp %>%
      head(20) %>%
      ggplot(aes(reorder(Feature, Overall), Overall)) +
      geom_col(fill = "forestgreen") +
      coord_flip() +
      labs(title = "Top 20 Predictive Genes",
           x = "Gene", y = "Importance Score") +
      theme_minimal()
    
    ggsave("results/feature_importance.png", imp_plot, width = 8, height = 6)
    
    return(list(
      predictions = preds,
      confusion_matrix = cm,
      variable_importance = var_imp,
      auc = auc(roc_obj)
    ))
    
  }, error = function(e) {
    message("Evaluation error: ", e$message)
    return(NULL)
  })
}