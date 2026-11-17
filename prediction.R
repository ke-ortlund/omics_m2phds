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
theme_set(theme_minimal(base_size = 14))      # default look for all plots
dir.create(here("output"), showWarnings = FALSE)

###############################################################################
#  1. DATA CHECKS BEFORE MODELLING
###############################################################################

scale_check <- flow %>%
  filter(study_accession %in% c("2024_UGA", "SDY2867"),
         parent_population %in% c("CD45plus_cells", "CD45+ cells"),
         timepoint == "7") %>%
  group_by(study_accession, participant_id) %>%
  summarise(total = sum(value), .groups = "drop") %>%
  group_by(study_accession) %>%
  summarise(median_total = median(total),
            min_total    = min(total),
            max_total    = max(total))
scale_check   # print it to check
write.csv(scale_check, here("output", "flow_scale_check.csv"), row.names = FALSE)

### 2a. One day-7 value per participant? (should return 0 rows)
flow_mon_clean %>% count(participant_id) %>% filter(n > 1)

### 2b. Are outcome values on the same scale in both studies?
flow_mon_clean %>%
  group_by(study_accession) %>%
  summarise(min = min(value), median = median(value), max = max(value))

# Plot 
x <- ggplot(flow_mon_clean, aes(study_accession, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(x = "Study", y = "Day-7 classical monocytes (proportion of CD45+ cells)")

# Save
ggsave(x, path)

### 2c. RNA-seq: material, raw counts available, duplicates
rna_mon %>% count(study_accession, material)
rna_mon %>% group_by(study_accession) %>% summarise(na_counts = sum(is.na(counts)))
rna_mon %>% count(participant_id, ensembl_gene_id) %>% filter(n > 1)   # should be 0 rows


###############################################################################
#  3. ANALYSIS DATASET: OUTCOME + COVARIATES
###############################################################################

dat <- flow_mon_clean %>%
  select(participant_id, study_accession, mono_d7 = value) %>%
  left_join(dem_clean %>% select(participant_id, age, sex),   # CHECK column names in dem_clean
            by = "participant_id") %>%
  arrange(participant_id)                                     # fixed order for matching later

summary(dat)   # check for NAs in age / sex
nrow(dat)      # 77


###############################################################################
#  4. RNA-SEQ PREPROCESSING
###############################################################################

### 4a. Genes measured in BOTH studies
# A gene missing from one study would otherwise be filled with 0 for all its
# participants, creating a fake "study" signal.
common_genes <- rna_mon %>%
  distinct(study_accession, ensembl_gene_id) %>%
  count(ensembl_gene_id) %>%
  filter(n == n_distinct(rna_mon$study_accession)) %>%
  pull(ensembl_gene_id)

length(common_genes)   # report this number in the methods

### 4b. Raw count matrix: genes (rows) x participants (columns)
counts_mat <- rna_mon %>%
  filter(participant_id %in% dat$participant_id,
         ensembl_gene_id %in% common_genes) %>%
  select(participant_id, ensembl_gene_id, counts) %>%
  pivot_wider(names_from  = participant_id,
              values_from = counts,
              values_fn   = sum) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

counts_mat <- counts_mat[, dat$participant_id]   # same participant order as dat
dim(counts_mat)                                  # genes x 77

### 4c. Filter low-expressed genes, TMM-normalise, log2-CPM
dge    <- DGEList(counts = counts_mat, group = dat$study_accession)
keep   <- filterByExpr(dge, group = dat$study_accession)
dge    <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
logcpm <- cpm(dge, log = TRUE)

sum(keep)   # number of genes after filtering (report in methods)

### 4d. Keep the 1,000 most variable genes
# Unsupervised: uses gene variance only, never the outcome, so no leakage.
top_genes <- names(sort(apply(logcpm, 1, var), decreasing = TRUE))[1:1000]
X_genes   <- t(logcpm[top_genes, ])   # participants x genes

dim(X_genes)                                  # 77 x 1000
all(rownames(X_genes) == dat$participant_id)  # MUST be TRUE


###############################################################################
#  5. OUTCOME TRANSFORMATION + PREDICTOR SETS
###############################################################################

### 5a. Log-transform the outcome (proportions are right-skewed), then standardise
dat <- dat %>% mutate(y_log = log2(mono_d7))
hist(dat$y_log, main = "log2 day-7 classical monocytes", xlab = "")

y_vec <- as.numeric(scale(dat$y_log))   # mean 0, SD 1 -> sRMSE < 1 beats the mean

### 5b. Covariates as a numeric matrix (sex and study become 0/1 columns)
X_cov <- model.matrix(~ age + sex + study_accession, data = dat)[, -1] %>% scale()
stopifnot(nrow(X_cov) == length(y_vec))   # fails if NAs dropped rows

### 5c. Standardised genes
X_genes_s <- scale(X_genes)

### 5d. Three predictor sets to compare
X_A <- X_cov                      # A: reference (age, sex, study)
X_B <- X_genes_s                  # B: genes only
X_C <- cbind(X_cov, X_genes_s)    # C: genes + age, sex, study

n <- length(y_vec)


###############################################################################
#  6. HELPER FUNCTIONS
###############################################################################

### 6a. Performance metrics
calculate_metrics <- function(true, predicted) {
  list(
    sRMSE    = sqrt(mean((true - predicted)^2)) / sd(true),  # < 1 = better than mean
    R2       = cor(true, predicted, method = "pearson")^2,   # variance explained
    Spearman = cor(true, predicted, method = "spearman")     # ranking accuracy
  )
}

### 6b. Observed vs predicted plot, coloured by study
plot_predictions <- function(true, predicted, study, title) {
  lims <- range(c(true, predicted))
  ggplot(data.frame(true, predicted, study),
         aes(true, predicted, colour = study)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = title,
         x = "Observed (standardised log2)", y = "Predicted", colour = "Study")
}

### 6c. Nested cross-validation for the elastic net
# Outer loop: evaluates performance on unseen participants.
# Inner loop (5-fold, inside caret::train): tunes alpha and lambda on training data only.
nested_cv_enet <- function(X, y, folds, grid, penalty_factor = rep(1, ncol(X))) {
  n         <- length(y)
  pred      <- rep(NA_real_, n)
  best_tune <- vector("list", length(folds))
  selected  <- vector("list", length(folds))
  
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)
    
    fit <- train(
      x = X[train_idx, ], y = y[train_idx],
      method    = "glmnet",
      tuneGrid  = grid,
      trControl = trainControl(method = "cv", number = 5),
      metric    = "RMSE",
      penalty.factor = penalty_factor          # 0 = never penalised (always kept)
    )
    
    pred[test_idx] <- predict(fit, newdata = X[test_idx, ])
    best_tune[[i]] <- fit$bestTune
    
    # genes with non-zero coefficients in this fold (for stability selection)
    cf <- as.matrix(coef(fit$finalModel, s = fit$bestTune$lambda))
    selected[[i]] <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
  }
  
  list(pred = pred, best_tune = bind_rows(best_tune), selected = selected)
}


###############################################################################
#  7. OUTER FOLDS (stratified by study)
###############################################################################

# Stratifying on study puts both studies in every training and test set.
folds <- createFolds(factor(dat$study_accession), k = 10, list = TRUE)

sapply(folds, function(idx) table(dat$study_accession[idx]))   # check composition


###############################################################################
#  8. FIT THE THREE MODELS
###############################################################################

### Hyperparameter grid (lambda spread over orders of magnitude)
lambda_grid <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)
enet_grid   <- expand.grid(alpha = seq(0, 1, 0.25), lambda = lambda_grid)

### 8a. Model A: reference (linear regression, same outer folds, nothing to tune)
pred_A <- rep(NA_real_, n)
for (i in seq_along(folds)) {
  test_idx  <- folds[[i]]
  train_idx <- setdiff(seq_len(n), test_idx)
  df_train  <- data.frame(y = y_vec[train_idx], X_A[train_idx, ])
  df_test   <- data.frame(X_A[test_idx, ])
  pred_A[test_idx] <- predict(lm(y ~ ., data = df_train), newdata = df_test)
}

### 8b. Model B: genes only (elastic net)
res_B <- nested_cv_enet(X_B, y_vec, folds, enet_grid)

### 8c. Model C: genes + covariates (covariates unpenalised so always included)
pf_C  <- c(rep(0, ncol(X_cov)), rep(1, ncol(X_genes_s)))
res_C <- nested_cv_enet(X_C, y_vec, folds, enet_grid, penalty_factor = pf_C)

### 8d. Check tuned hyperparameters
# If alpha/lambda sit at the edge of the grid, extend the grid in that direction.
res_B$best_tune
res_C$best_tune


###############################################################################
#  9. COMPARE PERFORMANCE
###############################################################################

model_metrics <- bind_rows(
  "A: Reference (age, sex, study)"     = calculate_metrics(y_vec, pred_A),
  "B: Genes only (elastic net)"        = calculate_metrics(y_vec, res_B$pred),
  "C: Genes + reference (elastic net)" = calculate_metrics(y_vec, res_C$pred),
  .id = "model"
)
model_metrics
write.csv(model_metrics, here("output", "prediction_metrics.csv"), row.names = FALSE)

# Key comparison: C vs A. Does gene expression add anything beyond age, sex, study?

p_B <- plot_predictions(y_vec, res_B$pred, dat$study_accession,
                        "Elastic net (genes only): nested CV predictions")
p_C <- plot_predictions(y_vec, res_C$pred, dat$study_accession,
                        "Elastic net (genes + covariates): nested CV predictions")
ggsave(here("output", "pred_enet_B.pdf"), p_B, width = 8, height = 8)
ggsave(here("output", "pred_enet_C.pdf"), p_C, width = 8, height = 8)
p_C


###############################################################################
#  10. STABILITY SELECTION: which genes are chosen consistently?
###############################################################################

stability <- function(selected_list, exclude = colnames(X_cov)) {
  freq <- table(unlist(selected_list)) / length(selected_list)
  data.frame(gene = names(freq), prop = as.numeric(freq)) %>%
    filter(!gene %in% exclude) %>%     # drop covariates, keep genes
    arrange(desc(prop))
}

stab_C <- stability(res_C$selected)

### Map Ensembl IDs to gene symbols (if org.Hs.eg.db is installed)
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stab_C$symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = stab_C$gene, keytype = "ENSEMBL",
    column = "SYMBOL", multiVals = "first"
  )
} else {
  stab_C$symbol <- NA_character_
}

head(stab_C, 20)
write.csv(stab_C, here("output", "stability_selection_C.csv"), row.names = FALSE)

### Genes selected in at least half of the outer folds
p_stab <- stab_C %>%
  filter(prop >= 0.5) %>%
  mutate(label = coalesce(symbol, gene)) %>%
  ggplot(aes(reorder(label, prop), prop)) +
  geom_col() +
  coord_flip() +
  labs(x = "Gene", y = "Proportion of outer folds selected",
       title = "Genes selected in ≥50% of folds (model C)")
ggsave(here("output", "stability_C.pdf"), p_stab, width = 8, height = 8)
p_stab


###############################################################################
#  11. FINAL MODEL ON ALL DATA (interpretation only, NOT for reporting performance)
###############################################################################

final_C <- train(
  x = X_C, y = y_vec,
  method    = "glmnet",
  tuneGrid  = enet_grid,
  trControl = trainControl(method = "cv", number = 5),
  metric    = "RMSE",
  penalty.factor = pf_C
)

final_C$bestTune

cf_final <- as.matrix(coef(final_C$finalModel, s = final_C$bestTune$lambda))
final_coefs <- data.frame(feature = rownames(cf_final), coef = cf_final[, 1]) %>%
  filter(coef != 0, feature != "(Intercept)") %>%
  arrange(desc(abs(coef)))

head(final_coefs, 20)
# Positive coefficient: higher baseline expression -> higher day-7 monocytes


###############################################################################
#  12. SENSITIVITY ANALYSIS: train on SDY2867, test on 2024_UGA
###############################################################################

train_idx <- which(dat$study_accession == "SDY2867")
test_idx  <- which(dat$study_accession == "2024_UGA")

# Study is constant within the training set, so remove the study column
keep_cols   <- !grepl("study_accession", colnames(X_C))
X_C_nostudy <- X_C[, keep_cols]
pf_nostudy  <- pf_C[keep_cols]

fit_ext <- train(
  x = X_C_nostudy[train_idx, ], y = y_vec[train_idx],
  method    = "glmnet",
  tuneGrid  = enet_grid,
  trControl = trainControl(method = "cv", number = 5),
  metric    = "RMSE",
  penalty.factor = pf_nostudy
)

pred_ext <- predict(fit_ext, newdata = X_C_nostudy[test_idx, ])
calculate_metrics(y_vec[test_idx], pred_ext)
# Exploratory only: 24 test participants


