# PART_AB_Subset_100k_Enhanced.R
# Subset cleanVD2 to 100,000 rows, run enhanced ML for Part A (GVPS regression) and Part B (Pathogenicity classification),
# output predictions, risk/class assignments, summary tables, and ready-to-plot data/plots for dashboard use.
# All outputs are labeled with 'subset' in their names.

library(tidyverse)
library(tidymodels)
library(stringr)
library(yardstick)
library(ranger)
library(lightgbm)
library(bonsai)
library(xgboost)
library(vip)
library(fastshap)
library(ggplot2)
library(pROC)
library(ROCR)
library(reshape2)


# --- Load and subset data (no NAs in modeling columns) ---
cleanVD2_full <- readRDS("cleanVD2.rds")

# Columns required for both Part A and B
required_cols <- c(
  "rsid", "SIFT_score", "FATHMM_score", "PolyPhen_score", "CADD_PHRED",
  "polyphen_max", "revel_max", "spliceai_ds_max", "phylop",
  "sift_label", "fathmm_label", "polyphen_label"
)

# Drop rows with any NA in required columns
cleanVD2_complete <- cleanVD2_full %>%
  filter(if_all(all_of(required_cols), ~ !is.na(.)))


# Sample 100,000 from complete cases, ensuring unique rsid
set.seed(42)
cleanVD2_complete_unique <- cleanVD2_complete %>% distinct(rsid, .keep_all = TRUE)
cleanVD2_subset <- cleanVD2_complete_unique %>% slice_sample(n = min(100000, nrow(cleanVD2_complete_unique)))

saveRDS(cleanVD2_subset, "cleanVD2_subset.rds")
write.csv(cleanVD2_subset, "cleanVD2_subset.csv", row.names = FALSE)

# --- Part A: GVPS Regression ---
gvps_data <- cleanVD2_subset %>%
  select(rsid, SIFT_score, FATHMM_score, PolyPhen_score, CADD_PHRED,
         polyphen_max, revel_max, spliceai_ds_max, phylop) %>%
  mutate(
    SIFT_score = 1 - SIFT_score,
    GVPS_target = rowMeans(cbind(SIFT_score, FATHMM_score, PolyPhen_score), na.rm = TRUE)
  ) %>%
  drop_na(GVPS_target)

# Enhanced preprocessing recipe

# --- Check for missing variables before GVPS model fitting ---
gvps_vars <- c("SIFT_score", "FATHMM_score", "PolyPhen_score", "CADD_PHRED", "polyphen_max", "revel_max", "spliceai_ds_max", "phylop")
missing_gvps <- setdiff(gvps_vars, names(gvps_data))
if (length(missing_gvps) > 0) {
  cat("[GVPS] Missing variables in gvps_data:", paste(missing_gvps, collapse=", "), "\n")
} else {
  cat("[GVPS] All required variables present in gvps_data.\n")
}

gvps_recipe <- recipe(GVPS_target ~ SIFT_score + FATHMM_score + PolyPhen_score + CADD_PHRED + 
                                     polyphen_max + revel_max + spliceai_ds_max + phylop, 
                      data = gvps_data) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_interact(~ SIFT_score:PolyPhen_score + FATHMM_score:CADD_PHRED) %>%
  step_poly(CADD_PHRED, degree = 2) %>%
  step_normalize(all_numeric_predictors())

# Split data
set.seed(123)
split1 <- initial_split(gvps_data, prop = 0.8, strata = GVPS_target)
temp_data <- training(split1)
test_data <- testing(split1)
split2 <- initial_split(temp_data, prop = 0.75, strata = GVPS_target)
train_data <- training(split2)
validation_data <- testing(split2)

# Models
gvps_models <- list(
  rf = rand_forest(trees = 500, mtry = 2, min_n = 5) %>% set_engine("ranger") %>% set_mode("regression"),
  xgb = boost_tree(trees = 300, learn_rate = 0.1, tree_depth = 6) %>% set_engine("xgboost") %>% set_mode("regression"),
  lgb = boost_tree(trees = 300, learn_rate = 0.1, tree_depth = 6) %>% set_engine("lightgbm") %>% set_mode("regression")
)

# Fit models and save

fits <- list()
for (m in names(gvps_models)) {
  wf <- workflow() %>% add_recipe(gvps_recipe) %>% add_model(gvps_models[[m]])
  fits[[m]] <- fit(wf, data = train_data)
  saveRDS(fits[[m]], paste0("gvps_model_", m, "_subset.rds"))
}

# --- GVPS SHAP summary (Random Forest) ---

# --- Check for missing variables before GVPS SHAP ---
shap_gvps_vars <- setdiff(names(gvps_data), c("rsid", "GVPS_target"))
missing_shap_gvps <- setdiff(shap_gvps_vars, names(gvps_data %>% select(-rsid, -GVPS_target)))
if (length(missing_shap_gvps) > 0) {
  cat("[GVPS SHAP] Missing variables:", paste(missing_shap_gvps, collapse=", "), "\n")
} else {
  cat("[GVPS SHAP] All required variables present for SHAP.\n")
}


# Predict on all subset data with RF (for dashboard)
gvps_pred <- predict(fits$rf, new_data = gvps_data) %>% bind_cols(gvps_data %>% select(rsid))
cleanVD2_subset$GVPS_pred <- gvps_pred$.pred

# --- Dashboard helper functions (Part A) ---
format_GVPS_with_risk_level <- function(scores) {
  case_when(
    scores < 0.3 ~ "Low Risk",
    scores < 0.7 ~ "Medium Risk",
    TRUE ~ "High Risk"
  )
}
create_gvps_percentile <- function(scores) {
  ecdf_fn <- ecdf(scores)
  round(ecdf_fn(scores) * 100, 1)
}
flag_gvps_outliers <- function(scores) {
  z <- scale(scores)
  abs(z) > 2.5
}

cleanVD2_subset <- cleanVD2_subset %>%
  mutate(
    GVPS_risk = format_GVPS_with_risk_level(GVPS_pred),
    GVPS_percentile = create_gvps_percentile(GVPS_pred),
    GVPS_outlier = flag_gvps_outliers(GVPS_pred)
  )

# --- Part A Plots ---
p1 <- ggplot(cleanVD2_subset, aes(GVPS_pred)) +
  geom_histogram(bins = 40, fill = "#2c3e50", color = "white") +
  labs(title = "GVPS distribution", x = "GVPS", y = "Count") + theme_minimal()
ggsave("gvps_distribution_subset.png", p1)

p2 <- ggplot(cleanVD2_subset, aes(GVPS_risk, fill = GVPS_risk)) +
  geom_bar() + labs(title = "GVPS Risk Stratification", x = "Risk Level", y = "Count") + theme_minimal()
ggsave("gvps_risk_stratification_subset.png", p2)

p3 <- ggplot(cleanVD2_subset, aes(GVPS_pred, GVPS_percentile)) +
  geom_point(alpha = 0.3) + labs(title = "GVPS Percentile Scatter", x = "GVPS", y = "Percentile") + theme_minimal()
ggsave("gvps_percentile_scatter_subset.png", p3)

p4 <- ggplot(cleanVD2_subset, aes(GVPS_pred, color = GVPS_outlier)) +
  geom_density() + labs(title = "GVPS Outlier Density", x = "GVPS", y = "Density") + theme_minimal()
ggsave("gvps_outlier_density_subset.png", p4)

# --- Part B: Pathogenicity Classification ---
label_to_sift <- function(label) {
  case_when(
    is.na(label) ~ NA_real_,
    label == "tolerated" ~ 0.0,
    label == "deleterious" ~ 1.0,
    TRUE ~ NA_real_
  )
}
label_to_polyphen <- function(label) {
  case_when(
    is.na(label) ~ NA_real_,
    label == "benign" ~ 0.0,
    label == "possibly_damaging" ~ 0.5,
    label == "probably_damaging" ~ 1.0,
    TRUE ~ NA_real_
  )
}
label_to_fathmm <- function(label) {
  case_when(
    is.na(label) ~ NA_real_,
    label == "tolerated" ~ 0.0,
    label == "deleterious" ~ 1.0,
    TRUE ~ NA_real_
  )
}


# Pathogenicity modeling data: ensure unique rsid and Pathogenicity_Level present
path_data <- cleanVD2_subset %>%
  select(rsid, SIFT_score, FATHMM_score, PolyPhen_score, sift_label, fathmm_label, polyphen_label,
         CADD_PHRED, polyphen_max, revel_max, spliceai_ds_max, phylop) %>%
  mutate(
    sift_sev = label_to_sift(sift_label),
    polyphen_sev = label_to_polyphen(polyphen_label),
    fathmm_sev = label_to_fathmm(fathmm_label),
    avg_sev = rowMeans(cbind(sift_sev, polyphen_sev, fathmm_sev), na.rm = TRUE),
    Pathogenicity_Level = case_when(
      is.na(avg_sev) ~ NA_character_,
      avg_sev < 0.15 ~ "Benign",
      avg_sev < 0.50 ~ "Possibly_Deleterious",
      avg_sev < 0.85 ~ "Probably_Deleterious",
      TRUE ~ "Highly_Deleterious"
    ),
    Pathogenicity_Level = factor(Pathogenicity_Level,
      levels = c("Benign", "Possibly_Deleterious", "Probably_Deleterious", "Highly_Deleterious"))
  ) %>%
  select(-sift_sev, -polyphen_sev, -fathmm_sev, -avg_sev, -sift_label, -fathmm_label, -polyphen_label) %>%
  drop_na(Pathogenicity_Level) %>%
  distinct(rsid, .keep_all = TRUE)

# Enhanced preprocessing
# --- Print column names for debugging Pathogenicity ---
cat("[Pathogenicity] Names in path_data:\n"); print(names(path_data))

# --- Check for missing variables before Pathogenicity model fitting ---
path_vars <- c("SIFT_score", "FATHMM_score", "PolyPhen_score", "CADD_PHRED", "polyphen_max", "revel_max", "spliceai_ds_max", "phylop")
missing_path <- setdiff(path_vars, names(path_data))
if (length(missing_path) > 0) {
  cat("[Pathogenicity] Missing variables in path_data:", paste(missing_path, collapse=", "), "\n")
} else {
  cat("[Pathogenicity] All required variables present in path_data.\n")
}

path_recipe <- recipe(Pathogenicity_Level ~ SIFT_score + FATHMM_score + PolyPhen_score + 
                                            CADD_PHRED + polyphen_max + revel_max + 
                                            spliceai_ds_max + phylop,
                      data = path_data) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_interact(~ SIFT_score:PolyPhen_score + FATHMM_score:CADD_PHRED) %>%
  step_poly(CADD_PHRED, degree = 2) %>%
  step_normalize(all_numeric_predictors())

# Split data
set.seed(123)
split1 <- initial_split(path_data, prop = 0.8, strata = Pathogenicity_Level)
temp_data <- training(split1)
test_data <- testing(split1)
split2 <- initial_split(temp_data, prop = 0.75, strata = Pathogenicity_Level)
train_data <- training(split2)
validation_data <- testing(split2)

# --- Print column names for debugging Pathogenicity train_data ---
cat("[Pathogenicity] Names in train_data:\n"); print(names(train_data))

# Models
path_models <- list(
  rf = rand_forest(trees = 500, mtry = 2, min_n = 5) %>% set_engine("ranger") %>% set_mode("classification"),
  xgb = boost_tree(trees = 300, learn_rate = 0.1, tree_depth = 6) %>% set_engine("xgboost") %>% set_mode("classification"),
  lgb = boost_tree(trees = 300, learn_rate = 0.1, tree_depth = 6) %>% set_engine("lightgbm") %>% set_mode("classification")
)

# Fit models and save
path_fits <- list()
for (m in names(path_models)) {
  wf <- workflow() %>% add_recipe(path_recipe) %>% add_model(path_models[[m]])
  path_fits[[m]] <- fit(wf, data = train_data)
  saveRDS(path_fits[[m]], paste0("pathogenicity_model_", m, "_subset.rds"))
}

# --- Pathogenicity SHAP summary (Random Forest) ---

# --- Check for missing variables before Pathogenicity SHAP ---
shap_path_vars <- setdiff(names(path_data), c("rsid", "Pathogenicity_Level"))
missing_shap_path <- setdiff(shap_path_vars, names(path_data %>% select(-rsid, -Pathogenicity_Level)))
if (length(missing_shap_path) > 0) {
  cat("[Pathogenicity SHAP] Missing variables:", paste(missing_shap_path, collapse=", "), "\n")
} else {
  cat("[Pathogenicity SHAP] All required variables present for SHAP.\n")
}


# Predict on all subset data with RF (for dashboard)
path_pred <- predict(path_fits$rf, new_data = path_data, type = "prob")
path_pred_class <- predict(path_fits$rf, new_data = path_data)

# Ensure df_pred has unique rsid for joining
df_pred <- bind_cols(path_data %>% select(rsid), path_pred, .pred_class = path_pred_class$.pred_class) %>%
  distinct(rsid, .keep_all = TRUE)
cleanVD2_subset <- left_join(cleanVD2_subset, df_pred, by = "rsid")

# --- Part B Plots ---

# For plotting, use path_data (guaranteed to have Pathogenicity_Level)
pb1 <- ggplot(path_data, aes(Pathogenicity_Level)) +
  geom_bar(fill = "#2980b9") + labs(title = "Pathogenicity Class Distribution", x = "Class", y = "Count") + theme_minimal()
ggsave("pathogenicity_class_distribution_subset.png", pb1)



# Confusion matrix (on test set)
test_pred_class <- predict(path_fits$rf, new_data = test_data)
test_results <- tibble(
  truth = test_data$Pathogenicity_Level,
  prediction = test_pred_class$.pred_class
)
cm <- conf_mat(test_results, truth = truth, estimate = prediction)
pb2 <- autoplot(cm, type = "heatmap") + ggtitle("Confusion Matrix (Test Set)")
ggsave("pathogenicity_confusion_matrix_subset.png", pb2)

# ROC/PR curves (one-vs-all, test set)
roc_list <- list()
for (cls in levels(test_data$Pathogenicity_Level)) {
  prob_col <- paste0(".pred_", cls)
  if (prob_col %in% colnames(path_pred)) {
    roc_obj <- roc(as.numeric(test_data$Pathogenicity_Level == cls), predict(path_fits$rf, new_data = test_data, type = "prob")[[prob_col]])
    roc_list[[cls]] <- roc_obj
    png(paste0("roc_curve_", cls, "_subset.png"))
    plot(roc_obj, main = paste("ROC Curve -", cls))
    dev.off()
  }
}

# Save final combined subset with predictions
saveRDS(cleanVD2_subset, "cleanVD2_GVPS_Pathogenicity_Combined_subset.rds")
write.csv(cleanVD2_subset, "cleanVD2_GVPS_Pathogenicity_Combined_subset.csv", row.names = FALSE)
