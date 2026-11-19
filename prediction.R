###############################################################################
#  0. SETUP
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(caret)    # nested cross-validation (train, createFolds)
library(glmnet)   # elastic net
library(edgeR)    # RNA-seq filtering + TMM normalisation
library(here)     # project-relative file paths

set.seed(123)                                 # reproducibility of folds / CV

###############################################################################
#  1. ANALYSIS DATASET: one row per participant
###############################################################################

### How many rows per participant? (more than 1 = several vaccine-history rows)
mon_dem_clean %>% count(participant_id) %>% count(n, name = "n_participants")

### Collapse to one row per participant
pred_dat <- mon_dem_clean %>%
  group_by(participant_id, study_accession) %>%
  summarise(
    monocytes    = first(value_clean),     # day-7 classical monocytes (% of CD45+)
    age          = first(age),             # CHECK column name in dem
    sex          = first(sex),             # CHECK column name in dem
    # 1 if ANY recorded prior influenza vaccination, otherwise 0
    vacc_history = as.integer(any(vacc_binary == 1, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(participant_id)                  # same order as the columns of logcpm

nrow(pred_dat)                             # 77 (fewer if some lack vaccine history)
pred_dat %>% count(study_accession)        # 53 SDY2867, 24 2024_UGA
summary(pred_dat)                          # check for NAs

### Match the expression matrix to the participants in pred_dat
ids <- pred_dat$participant_id
stopifnot(all(ids %in% colnames(logcpm)))  # every participant has RNA-seq
logcpm_pred <- logcpm[, ids]               # same participants, same order

###############################################################################
#  2. OUTCOME: day-7 classical monocyte frequency (not log-transformed)
###############################################################################

### Distribution by study
p_outcome <- ggplot(pred_dat, aes(monocytes)) +
  geom_histogram(bins = 15) +
  facet_wrap(~ study_accession, scales = "free_y") +
  labs(x = "Day-7 classical monocytes (% of CD45+ cells)", y = "Participants",
       title = "Distribution of the outcome by study")
ggsave(file.path(fig_folder, "outcome_distribution.png"), p_outcome,
       width = 9, height = 5, dpi = 300)
p_outcome

### Standardise the pooled outcome (mean 0, SD 1)
y_vec <- as.numeric(scale(pred_dat$monocytes))

###############################################################################
#  3. PREDICTORS
###############################################################################

### 3a. Genes: participants (rows) x genes (columns), all genes after filterByExpr
X_genes <- t(logcpm_pred)

# remove genes with zero variance (they cannot be standardised)
X_genes <- X_genes[, apply(X_genes, 2, sd) > 0]

# standardise every gene (mean 0, SD 1), as in the practical
X_genes_s <- scale(X_genes)
dim(X_genes_s)                                   # 77 x number of genes

### 3b. Covariates: age, sex, vaccine history, study
# model.matrix turns sex and study into 0/1 columns
X_cov <- model.matrix(~ age + sex + vacc_history + study_accession,
                      data = pred_dat)[, -1, drop = FALSE] %>%
  scale()

# model.matrix silently drops rows with NA: check nothing was lost
stopifnot(nrow(X_cov) == length(y_vec))
colnames(X_cov)

### 3c. Three predictor sets
X_A <- X_cov                       # A: covariates only (reference)
X_B <- X_genes_s                   # B: genes only
X_C <- cbind(X_cov, X_genes_s)     # C: genes + covariates (main model)

### 3d. Penalty factor for model C
# 0 = covariate is never penalised (always kept in the model)
# 1 = gene is penalised as usual (can be shrunk to 0)
pf_C <- c(rep(0, ncol(X_cov)), rep(1, ncol(X_genes_s)))

n <- length(y_vec)

###############################################################################
#  4. HELPER FUNCTIONS
###############################################################################

### 4a. Performance metrics (from the practical)
# sRMSE: RMSE divided by SD of the outcome. < 1 = better than predicting the mean.
# R2   : squared Pearson correlation between observed and predicted values.
calculate_metrics <- function(true, predicted) {
  sRMSE <- sqrt(mean((true - predicted)^2)) / sd(true)
  R2    <- cor(true, predicted, method = "pearson")^2
  list(sRMSE = sRMSE, R2 = R2)
}

### 4b. Observed vs predicted plot (from the practical, coloured by study)
plot_predictions <- function(true, predicted, study, title) {
  rng  <- range(c(true, predicted))
  pad  <- diff(rng) * 0.02
  lims <- c(rng[1] - pad, rng[2] + pad)
  
  ggplot(data.frame(true, predicted, study), aes(true, predicted, colour = study)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = title, x = "Observed (standardised)", y = "Predicted", colour = "Study")
}

### 4c. Nested cross-validation for any caret model
# Outer loop : each fold is held out once as test data -> unbiased predictions.
# Inner loop : 5-fold CV on the training data only -> chooses hyperparameters.
# "..." passes extra arguments to the model (e.g. penalty.factor for glmnet).
nested_cv <- function(X, y, folds, method, grid, ...) {
  pred      <- rep(NA_real_, length(y))
  best_tune <- vector("list", length(folds))
  selected  <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    # tune and fit on the training participants only
    fit <- train(
      x         = X[train_idx, , drop = FALSE],
      y         = y[train_idx],
      method    = method,
      tuneGrid  = grid,
      trControl = trainControl(method = "cv", number = 5),
      metric    = "RMSE",
      ...
    )
    
    # predict the held-out participants
    pred[test_idx] <- predict(fit, newdata = X[test_idx, , drop = FALSE])
    best_tune[[i]] <- fit$bestTune
    
    # for glmnet models: store variables with non-zero coefficients
    if (method == "glmnet") {
      cf <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))
      selected[[i]] <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
    }
  }
  
  list(pred = pred, best_tune = bind_rows(best_tune), selected = selected)
}

###############################################################################
#  5. OUTER CROSS-VALIDATION FOLDS
###############################################################################

folds <- createFolds(factor(pred_dat$study_accession), k = 10, list = TRUE)

# check: both studies in every fold
sapply(folds, function(idx) table(pred_dat$study_accession[idx]))

###############################################################################
#  6. REFERENCE MODEL: covariates only (linear regression)
###############################################################################

pred_ref <- rep(NA_real_, n)

for (i in seq_along(folds)) {
  test_idx  <- folds[[i]]
  train_idx <- setdiff(seq_len(n), test_idx)
  
  df_train <- data.frame(y = y_vec[train_idx], X_A[train_idx, , drop = FALSE])
  df_test  <- data.frame(X_A[test_idx, , drop = FALSE])
  
  fit_lm <- lm(y ~ ., data = df_train)
  pred_ref[test_idx] <- predict(fit_lm, newdata = df_test)
}

calculate_metrics(y_vec, pred_ref)

###############################################################################
#  7. MAIN MODEL: ELASTIC NET
###############################################################################

### Hyperparameter grid
# alpha: 0 = ridge, 1 = lasso, in between = elastic net
# lambda: penalty strength, spread over several orders of magnitude
lambda_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 50)
enet_grid   <- expand.grid(alpha = seq(0, 1, 0.25), lambda = lambda_grid)

### 7a. Elastic net: genes only
res_enet_B <- nested_cv(X_B, y_vec, folds, method = "glmnet", grid = enet_grid)

### 7b. Elastic net: genes + covariates (MAIN MODEL)
res_enet_C <- nested_cv(X_C, y_vec, folds, method = "glmnet", grid = enet_grid,
                        penalty.factor = pf_C)

### Check the tuned hyperparameters in each outer fold
# If lambda is often the largest value in the grid, extend the grid upwards.
res_enet_B$best_tune
res_enet_C$best_tune

###############################################################################
#  8. SENSITIVITY ANALYSES: OTHER MODELS (predictors: genes + covariates)
###############################################################################

### 8a. Ridge regression (alpha = 0): keeps all genes, shrinks coefficients
res_ridge <- nested_cv(X_C, y_vec, folds, method = "glmnet",
                       grid = expand.grid(alpha = 0, lambda = lambda_grid),
                       penalty.factor = pf_C)

### 8b. Lasso regression (alpha = 1): sets most gene coefficients exactly to 0
res_lasso <- nested_cv(X_C, y_vec, folds, method = "glmnet",
                       grid = expand.grid(alpha = 1, lambda = lambda_grid),
                       penalty.factor = pf_C)

### 8c. Principal component regression: regression on the first PCs of the predictors
ncomp_grid <- data.frame(ncomp = c(1, 2, 5, 10, 20))
res_pcr <- nested_cv(X_C, y_vec, folds, method = "pcr", grid = ncomp_grid)

### 8d. Partial least squares: components chosen to covary with the outcome
res_pls <- nested_cv(X_C, y_vec, folds, method = "pls", grid = ncomp_grid)

### 8e. Random forest: non-linear ensemble of decision trees
p <- ncol(X_C)
rf_grid <- expand.grid(
  mtry          = c(floor(sqrt(p)), floor(p / 10)),   # genes tried at each split
  splitrule     = "variance",                         # regression
  min.node.size = c(5, 10)
)
res_rf <- nested_cv(X_C, y_vec, folds, method = "ranger", grid = rf_grid)

###############################################################################
#  9. PERFORMANCE OF ALL MODELS
###############################################################################

### Collect the cross-validated predictions of every model
all_preds <- list(
  "Reference: age, sex, vaccine history, study" = pred_ref,
  "Elastic net: genes only"                     = res_enet_B$pred,
  "Elastic net: genes + covariates (main)"      = res_enet_C$pred,
  "Ridge: genes + covariates"                   = res_ridge$pred,
  "Lasso: genes + covariates"                   = res_lasso$pred,
  "PCR: genes + covariates"                     = res_pcr$pred,
  "PLS: genes + covariates"                     = res_pls$pred,
  "Random forest: genes + covariates"           = res_rf$pred
)

### Metrics table: one row per model
model_metrics <- all_preds %>%
  map(~ calculate_metrics(y_vec, .x)) %>%
  bind_rows(.id = "model") %>%
  mutate(across(c(sRMSE, R2), ~ round(.x, 3)))

model_metrics
write.csv(model_metrics, file.path(tab_folder, "prediction_model_metrics.csv"),
          row.names = FALSE)

### Observed vs predicted plot for every model, saved to output/figures
iwalk(all_preds, function(pred, model_name) {
  p <- plot_predictions(y_vec, pred, pred_dat$study_accession,
                        paste0(model_name, "\n(nested cross-validation)"))
  file_name <- paste0("pred_", gsub("[^A-Za-z0-9]+", "_", model_name), ".png")
  ggsave(file.path(fig_folder, file_name), p, width = 8, height = 8, dpi = 300)
})

### Tuned hyperparameters per model (useful for the methods/appendix)
best_tunes <- list(
  "Elastic net: genes only"     = res_enet_B$best_tune,
  "Elastic net: genes + covs"   = res_enet_C$best_tune,
  "Ridge"                       = res_ridge$best_tune,
  "Lasso"                       = res_lasso$best_tune,
  "PCR"                         = res_pcr$best_tune,
  "PLS"                         = res_pls$best_tune,
  "Random forest"               = res_rf$best_tune
) %>% bind_rows(.id = "model")

write.csv(best_tunes, file.path(tab_folder, "prediction_best_hyperparameters.csv"),
          row.names = FALSE)

###############################################################################
#  10. INTERPRETATION: GENES SELECTED BY THE MAIN ELASTIC NET
###############################################################################

### 10a. Stability selection: proportion of outer folds selecting each gene
selection_freq <- table(unlist(res_enet_C$selected)) / length(folds)

stability <- data.frame(gene = names(selection_freq),
                        prop = as.numeric(selection_freq)) %>%
  filter(!gene %in% colnames(X_cov)) %>%          # keep genes, drop covariates
  arrange(desc(prop))

# Optional: Ensembl ID -> gene symbol (needs Bioconductor package org.Hs.eg.db)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stability$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                            keys = stability$gene,
                                            keytype = "ENSEMBL", column = "SYMBOL",
                                            multiVals = "first")
} else {
  stability$symbol <- NA_character_
}

head(stability, 20)
write.csv(stability, file.path(tab_folder, "enet_gene_stability.csv"), row.names = FALSE)

# Genes selected in at least half of the folds
p_stab <- stability %>%
  filter(prop >= 0.5) %>%
  mutate(label = coalesce(symbol, gene)) %>%
  ggplot(aes(reorder(label, prop), prop)) +
  geom_col() +
  coord_flip() +
  labs(x = "Gene", y = "Proportion of outer folds selected",
       title = "Genes selected in ≥50% of folds (main elastic net)")
ggsave(file.path(fig_folder, "enet_gene_stability.png"), p_stab,
       width = 8, height = 8, dpi = 300)
p_stab

### 10b. Final model on all participants (coefficients for interpretation only)
final_enet <- train(
  x = X_C, y = y_vec,
  method    = "glmnet",
  tuneGrid  = enet_grid,
  trControl = trainControl(method = "cv", number = 5),
  metric    = "RMSE",
  penalty.factor = pf_C
)

final_enet$bestTune

final_coefs <- as.matrix(coef(final_enet$finalModel, s = final_enet$bestTune$lambda)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature") %>%
  rename(coefficient = 2) %>%
  filter(coefficient != 0, feature != "(Intercept)") %>%
  arrange(desc(abs(coefficient)))

head(final_coefs, 20)
write.csv(final_coefs, file.path(tab_folder, "enet_final_coefficients.csv"), row.names = FALSE)
# Positive coefficient: higher baseline expression -> higher day-7 monocyte frequency