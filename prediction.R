###############################################################################
#  PREDICTION: day-7 classical monocyte frequency from baseline gene expression
###############################################################################


###############################################################################
#  0. SETUP
###############################################################################
pacman::p_load(caret, glmnet, pls, ranger)

set.seed(123)

fig_folder <- file.path(output_folder, "figures")
tab_folder <- file.path(output_folder, "tables")
dir.create(fig_folder, showWarnings = FALSE)
dir.create(tab_folder, showWarnings = FALSE)

###############################################################################
#  1. ANALYSIS DATASET: one row per participant
###############################################################################
## performing some checks (one row per patient, missingness) and removing 
## vaccine history from the data since it is uninformative

pred_dat <- mon_dem_clean %>%
  transmute(
    participant_id,
    study_accession,
    monocytes    = value_clean,                   
    age,
    sex          = biological_sex,
    vacc_history = coalesce(vacc_any, 0L)         
  ) %>%
  arrange(participant_id)                         

### Checks
nrow(pred_dat)                                    
n_distinct(pred_dat$participant_id)               
pred_dat %>% count(study_accession)               
pred_dat %>% count(study_accession, vacc_history) 
summary(pred_dat)                                 

### Expression matrix in the same participant order
stopifnot(all(colnames(logcpm) == pred_dat$participant_id))


###############################################################################
#  2. OUTCOME: day-7 classical monocyte frequency (pooled, not log-transformed)
###############################################################################

### Distribution by study (for the descriptive results)
p_outcome <- ggplot(pred_dat, aes(monocytes)) +
  geom_histogram(bins = 15) +
  facet_wrap(~ study_accession, scales = "free_y") +
  labs(x = "Day-7 classical monocytes (% of CD45+ cells)", y = "Participants",
       title = "Distribution of the outcome by study")
ggsave(file.path(fig_folder, "outcome_distribution.png"), p_outcome,
       width = 9, height = 5, dpi = 300)
p_outcome

### Standardize the pooled outcome (mean 0, SD 1), as in the practical
# sRMSE can then be read directly: 1 = no better than predicting the mean
y_vec <- as.numeric(scale(pred_dat$monocytes))
n     <- length(y_vec)


###############################################################################
#  3. PREDICTORS
###############################################################################

### 3a. Genes: participants (rows) x genes (columns), all filtered genes
## dropping zero-variance genes, standardizing genes
X_genes <- t(logcpm)
X_genes <- X_genes[, apply(X_genes, 2, sd) > 0]  
X_genes_s <- scale(X_genes)                      
dim(X_genes_s)                                    

### 3b. Covariates: age, sex, study (sex and study become 0/1)
X_cov <- model.matrix(~ age + sex + study_accession,
                      data = pred_dat)[, -1, drop = FALSE]

# safety check: drop any covariate that is constant (cannot be standardised)
X_cov <- X_cov[, apply(X_cov, 2, sd) > 0, drop = FALSE]
X_cov <- scale(X_cov)

stopifnot(nrow(X_cov) == n)                       
colnames(X_cov)                                   

### 3c. Three predictor sets
## creating 3 predictor sets, A (covariates only), B(genes only), C(genes+covariates)
X_A <- X_cov                       
X_B <- X_genes_s                   
X_C <- cbind(X_cov, X_genes_s)     

### 3d. Penalty factor for glmnet models using X_C
# 0 = covariate never penalised (always kept); 1 = gene penalised as usual
pf_C <- c(rep(0, ncol(X_cov)), rep(1, ncol(X_genes_s)))


###############################################################################
#  4. HELPER FUNCTIONS
###############################################################################

### 4a. Performance metrics (from the practical)
calculate_metrics <- function(true, predicted) {
  sRMSE <- sqrt(mean((true - predicted)^2)) / sd(true)   # < 1 = better than the mean
  R2    <- cor(true, predicted, method = "pearson")^2    # variance explained
  list(sRMSE = sRMSE, R2 = R2)
}

### 4b. Observed vs predicted plot (from the practical, coloured by study)
plot_predictions <- function(true, predicted, study, title = NULL) {
  rng  <- range(c(true, predicted))
  pad  <- diff(rng) * 0.02
  lims <- c(floor(rng[1]), ceiling(rng[2])) 
  
  ggplot(data.frame(true, predicted, study), aes(true, predicted, colour = study)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    coord_equal(xlim = lims, ylim = lims) +
    scale_colour_discrete(labels = study_labels) +  
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + # nicer legend names
    labs(title = title, x = "Observed", y = "Predicted", colour = "Study") 
}

### 4c. Nested cross-validation for any caret model
# Outer loop: each fold held out once -> honest predictions for every participant
# Inner loop: 5-fold CV on the training part only -> picks hyperparameters
nested_cv <- function(X, y, folds, method, grid, ...) {
  pred      <- rep(NA_real_, length(y))
  best_tune <- vector("list", length(folds))
  selected  <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    fit <- train(
      x         = X[train_idx, , drop = FALSE],
      y         = y[train_idx],
      method    = method,
      tuneGrid  = grid,
      trControl = trainControl(method = "cv", number = 5),
      metric    = "RMSE",
      ...                                   # e.g. penalty.factor for glmnet
    )
    
    pred[test_idx] <- predict(fit, newdata = X[test_idx, , drop = FALSE])
    best_tune[[i]] <- fit$bestTune
    
    # glmnet only: variables with non-zero coefficients in this fold
    if (method == "glmnet") {
      cf <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))
      selected[[i]] <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
    }
  }
  
  list(pred = pred, best_tune = bind_rows(best_tune), selected = selected)
}


###############################################################################
#  5. OUTER CROSS-VALIDATION FOLDS (stratified by study)
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
  
  pred_ref[test_idx] <- predict(lm(y ~ ., data = df_train), newdata = df_test)
}

calculate_metrics(y_vec, pred_ref)


###############################################################################
#  7. MAIN MODEL: ELASTIC NET
###############################################################################

### Hyperparameter grid
lambda_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 50)
enet_grid   <- expand.grid(alpha = seq(0, 1, 0.25), lambda = lambda_grid)

### 7a. Genes only
res_enet_B <- nested_cv(X_B, y_vec, folds, method = "glmnet", grid = enet_grid)

### 7b. Genes + covariates (MAIN MODEL)
res_enet_C <- nested_cv(X_C, y_vec, folds, method = "glmnet", grid = enet_grid,
                        penalty.factor = pf_C)

### Tuned hyperparameters per outer fold
# lambda often at the largest grid value -> extend the grid upwards
res_enet_B$best_tune
res_enet_C$best_tune


###############################################################################
#  8. SENSITIVITY ANALYSES: OTHER MODELS (genes + covariates)
###############################################################################

### 8a. Ridge (alpha = 0)
res_ridge <- nested_cv(X_C, y_vec, folds, method = "glmnet",
                       grid = expand.grid(alpha = 0, lambda = lambda_grid),
                       penalty.factor = pf_C)

### 8b. Lasso (alpha = 1)
res_lasso <- nested_cv(X_C, y_vec, folds, method = "glmnet",
                       grid = expand.grid(alpha = 1, lambda = lambda_grid),
                       penalty.factor = pf_C)

### 8d. Random forest
p <- ncol(X_C)
rf_grid <- expand.grid(
  mtry          = c(floor(sqrt(p)), floor(p / 10)),
  splitrule     = "variance",
  min.node.size = c(5, 10)
)
res_rf <- nested_cv(X_C, y_vec, folds, method = "ranger", grid = rf_grid)


###############################################################################
#  9. PERFORMANCE OF ALL MODELS
###############################################################################

### 9a. Cross-validated predictions of every model
# Each vector holds one out-of-sample prediction per participant (from nested CV)
all_preds <- list(
  "Reference: age, sex, study"             = pred_ref,
  "Elastic net: genes only"                = res_enet_B$pred,
  "Main: elastic net: genes + covariates" = res_enet_C$pred,
  "Ridge: genes + covariates"              = res_ridge$pred,
  "Lasso: genes + covariates"              = res_lasso$pred,
  "Random forest: genes + covariates"      = res_rf$pred
)

# check: does any model predict the same value for everyone? (sd = 0 -> R2 = NA)
sapply(all_preds, sd)


### 9b. Metrics table: one row per model
# sRMSE < 1 = better than predicting the mean; R2 = proportion of variance explained
model_metrics <- all_preds %>%
  map(~ calculate_metrics(y_vec, .x)) %>%
  bind_rows(.id = "model") %>%
  mutate(across(c(sRMSE, R2), ~ round(.x, 3)))

model_metrics
write.csv(model_metrics, file.path(tab_folder, "prediction_model_metrics.csv"),
          row.names = FALSE)


### 9c. Observed vs predicted plot for every model: store, show and save
pred_plots <- imap(all_preds, function(pred, model_name) {
  plot_predictions(y_vec, pred, pred_dat$study_accession)
})

# show each plot in the Plots pane (use the arrows to browse)
walk(pred_plots, print)

# save each plot as its own file
iwalk(pred_plots, function(p, model_name) {
  file_name <- paste0("pred_", gsub("[^A-Za-z0-9]+", "_", model_name), ".png")
  ggsave(file.path(fig_folder, file_name), p, width = 8, height = 7, dpi = 300)
})


### 9d. All models in one figure (one panel per model, metrics in panel titles)
all_preds_long <- imap_dfr(all_preds, function(pred, model_name) {
  m <- calculate_metrics(y_vec, pred)
  data.frame(
    model     = sprintf("%s\nsRMSE = %.2f, R² = %.2f", model_name, m$sRMSE, m$R2),
    observed  = y_vec,
    predicted = pred,
    study     = pred_dat$study_accession
  )
})

# keep panels in the same order as all_preds (not alphabetical)
all_preds_long$model <- factor(all_preds_long$model,
                               levels = unique(all_preds_long$model))

p_all <- ggplot(all_preds_long, aes(observed, predicted, colour = study)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  facet_wrap(~ model, ncol = 3) +
  coord_equal() +
  scale_colour_discrete(labels = study_labels) +
  labs(x = "Observed", y = "Predicted", colour = "Study",
       title = "Observed vs predicted day-7 classical monocyte frequency") +
  theme(strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(fig_folder, "pred_all_models.png"), p_all,
       width = 14, height = 10, dpi = 300)
p_all


### 9e. Tuned hyperparameters per outer fold (for the methods / appendix)
# glmnet models have alpha + lambda; random forest has mtry, splitrule, min.node.size.
# Columns that don't apply to a model are NA.
best_tunes <- list(
  "Elastic net: genes only"         = res_enet_B$best_tune,
  "Elastic net: genes + covariates" = res_enet_C$best_tune,
  "Ridge"                           = res_ridge$best_tune,
  "Lasso"                           = res_lasso$best_tune,
  "Random forest"                   = res_rf$best_tune
) %>%
  map(~ mutate(.x, fold = row_number())) %>%      # which outer fold each row is from
  bind_rows(.id = "model") %>%
  relocate(model, fold)

best_tunes
write.csv(best_tunes, file.path(tab_folder, "prediction_best_hyperparameters.csv"),
          row.names = FALSE)
# this shows that the model prefered a ridge regression because the alpha was 0 for all of the elastic net runs. This shows that more genes have smaller effects.

### 9f. Within study performance
# Correlation between observed and predicted values WITHIN each study.
# Shows whether models capture differences between individuals,
# not just the difference between the two studies.

### Correlation per model and per study
within_study <- imap_dfr(all_preds, function(pred, model_name) {
  pred_dat %>%
    mutate(observed = y_vec, predicted = pred) %>%
    group_by(study_accession) %>%
    summarise(
      model    = model_name,
      n        = n(),
      r_within = cor(observed, predicted),   # Pearson correlation within study
      .groups  = "drop"
    )
}) %>%
  mutate(r_within = round(r_within, 2)) %>%
  select(model, study_accession, n, r_within) %>%
  tidyr::pivot_wider(names_from  = study_accession,
                     values_from = c(n, r_within))

within_study
write.csv(within_study, file.path(tab_folder, "prediction_within_study_correlation.csv"),
          row.names = FALSE)


### 95% confidence intervals for the main model
# (uncertainty is large with n = 24 in UGA 2024)
main_within_ci <- pred_dat %>%
  mutate(observed = y_vec, predicted = res_enet_C$pred) %>%
  group_by(study_accession) %>%
  summarise(
    n     = n(),
    r     = cor(observed, predicted),
    lower = cor.test(observed, predicted)$conf.int[1],
    upper = cor.test(observed, predicted)$conf.int[2],
    p     = cor.test(observed, predicted)$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(c(r, lower, upper), ~ round(.x, 2)),
         p = signif(p, 2))

main_within_ci
write.csv(main_within_ci, file.path(tab_folder, "main_model_within_study_ci.csv"),
          row.names = FALSE)

###############################################################################
#  10. SENSITIVITY ANALYSIS: STUDY EFFECT REMOVED (within-study centring)
###############################################################################
# PCA showed that study dominates gene expression (PC1 = 62% of variance),
# and monocyte levels also differ between studies. Here the study effect is
# removed from BOTH the outcome and the genes, so the model can only learn
# differences between individuals within each study.

### 10a. Outcome: centre within study, then standardise (each study now has mean 0)
y_within <- pred_dat %>%
  mutate(y = monocytes) %>%
  group_by(study_accession) %>%
  mutate(y = y - mean(y)) %>%          
  ungroup() %>%
  pull(y)

y_within <- as.numeric(scale(y_within))   # overall SD = 1 (comparable sRMSE)

# check: both studies centred on 0
tapply(y_within, pred_dat$study_accession, mean)


### 10b. Genes: centre each gene within study, then standardise
X_genes_c <- X_genes                      
for (s in unique(pred_dat$study_accession)) {
  idx <- pred_dat$study_accession == s
  X_genes_c[idx, ] <- scale(X_genes[idx, ], center = TRUE, scale = FALSE)
}
X_genes_c <- scale(X_genes_c)
X_genes_c <- X_genes_c[, colSums(is.na(X_genes_c)) == 0]
dim(X_genes_c)


### 10c. Covariates: age and sex only (study no longer needed: effect removed)
X_cov_within <- model.matrix(~ age + sex, data = pred_dat)[, -1, drop = FALSE] %>%
  scale()

X_within  <- cbind(X_cov_within, X_genes_c)
pf_within <- c(rep(0, ncol(X_cov_within)), rep(1, ncol(X_genes_c)))


### 10d. Reference for comparison: age and sex only, same folds
pred_ref_within <- rep(NA_real_, n)
for (i in seq_along(folds)) {
  test_idx  <- folds[[i]]
  train_idx <- setdiff(seq_len(n), test_idx)
  df_train  <- data.frame(y = y_within[train_idx], X_cov_within[train_idx, , drop = FALSE])
  df_test   <- data.frame(X_cov_within[test_idx, , drop = FALSE])
  pred_ref_within[test_idx] <- predict(lm(y ~ ., data = df_train), newdata = df_test)
}


### 10e. Elastic net with study effect removed (same folds and grid as main model)
res_enet_within <- suppressWarnings(
  nested_cv(X_within, y_within, folds, method = "glmnet",
            grid = enet_grid, penalty.factor = pf_within)
)

res_enet_within$best_tune   # check tuning


### 10f. Performance
within_metrics <- list(
  "Reference (age, sex), study removed"            = pred_ref_within,
  "Elastic net (genes + age, sex), study removed"  = res_enet_within$pred
) %>%
  map(~ calculate_metrics(y_within, .x)) %>%
  bind_rows(.id = "model") %>%
  mutate(across(c(sRMSE, R2), ~ round(.x, 3)))

within_metrics
write.csv(within_metrics, file.path(tab_folder, "sensitivity_study_removed_metrics.csv"),
          row.names = FALSE)


### 10g. Observed vs predicted plot
p_within <- plot_predictions(
  y_within, res_enet_within$pred, pred_dat$study_accession,
  "Elastic net, study effect removed: observed vs predicted"
) + theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(fig_folder, "pred_enet_study_removed.png"), p_within,
       width = 8, height = 8, dpi = 300)
p_within