# Random Forest Model for Geographic Classification of *Salmonella enterica* Serotype Derby from Pork

# Multiclass Geographic Classification (China vs. France vs. Mexico vs. United States)

This repository contains R scripts for predicting the country of origin (*China*, *France*, *Mexico*, or *United States*) of pork-associated *Salmonella* Derby isolates using a random forest model trained on pan-genome gene presence/absence data.

## 🧠 Workflow Overview

1. **Data Loading & Preprocessing**  
   Loads isolate metadata and gene presence/absence matrix, merges by isolate, and removes near-zero variance genes.

2. **Model Training**  
   Trains a random forest classifier (750 trees, 5-fold CV) using *caret* and *ranger*.

3. **Model Evaluation**  
   Evaluates model performance on a test set using the area under the receiver operating characteristic curve, confusion matrix, and feature importance analysis.

## 📂 Scripts

- `01_main_workflow.R` — Main pipeline execution  
- `02_data_preprocessing.R` — Data loading and cleaning  
- `03_model_training.R` — Random forest training  
- `04_model_evaluation.R` — Model testing and visualization  

## 📊 Outputs

Generated in `Multiclass_Country_Results/`:
- `predictions.csv`
- `feature_importance.csv`
- `roc_curve.png`
- `feature_importance.png`
- `performance_summary.txt`

## 🧩 Dependencies

R packages: `tidyverse`, `caret`, `ranger`, `readxl`, `janitor`, `pROC`, `ggplot2`, `here`, `RColorBrewer`.

---

**Citation:**  
This workflow was developed for multiclass geographic classification of pork-associated *S.* Derby isolates based on pan-genomic features.  
If used in a publication, please cite the corresponding methods description.
