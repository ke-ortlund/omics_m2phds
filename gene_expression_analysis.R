###############################################################################
# STEP 1- VERIFY PARTICIPANT IDENTIFIERS
###############################################################################

# Check whether participant IDs appear in more than one study.

sample_info %>%
  count(participant_id) %>%
  filter(n > 1)

###############################################################################
# STEP 2 — VERIFY DATASET ALIGNMENT
###############################################################################

# Check the number of participants in each dataset
ncol(dge$counts)
nrow(sample_info)
n_distinct(mon_dem_clean$participant_id)

# Check that all RNA-seq participants have clinical data
all(sample_info$participant_id %in% mon_dem_clean$participant_id)

# Check whether clinical data contain duplicate participants
mon_dem_clean %>%
  count(participant_id) %>%
  filter(n > 1)

###############################################################################
# STEP 3 — EXPLORE AVAILABLE COVARIATES
###############################################################################

# List all available variables
names(mon_dem_clean)


###############################################################################
#### GENE EXPRESSION ANALYSIS
###############################################################################

# Study population:
# 77 participants from two studies: 2024_UGA and SDY2867.

# Data pre-work:
# RNA-seq filtering and TMM normalization have already been done
# in environment.R. Monocyte percentages have also been harmonized.
# These preprocessed data will be reused without repeating normalization.

# Additional considerations:
# PCA identified substantial study-related variation (PC1: 62.3%) and
# sex-related variation (PC2: 6.7%) among the 1,000 most variable genes.
# Batch effects will be investigated separately.
#
# Statistical significance threshold FDR < 0.05.

# Analysis plan:

          ## With voom–limma (model 0-4) ## identify genes whose baseline expression 
          ## is associated with subsequent monocyte responses:

# - Y: Baseline gene expression (Day 0).
# - X: Day 7 classical monocyte percentage.
# Method: TMM normalization → voom → limma → empirical Bayes.
# Multiple testing: Benjamini–Hochberg FDR.
# Model 0 — Unadjusted analysis
      #   Both studies, shared genes, no covariate adjustment.
# Model 1 — Adjusted analysis
      #   Both studies, shared genes, adjusted for study, sex, and age.
# Model 2 — Extended gene set analysis
      #   Exploratory analysis including genes excluded from the shared gene set.
      #   Population and covariates to be defined.
# Model 3 — Study-specific analysis: 2024_UGA
      #   Association analysis within the UGA cohort, adjusted for sex and age.
# Model 4 — Study-specific analysis: SDY2867
      #   Association analysis within the SDY2867 cohort, adjusted for sex and age.


          ## With   Gene-wise linear regression (model5) ##
          # Identify baseline genes whose expression is associated with
          # subsequent Day 7 classical monocyte responses.$

# - Y: Day 7 classical monocyte percentage.
# - X: Baseline gene expression (Day 0), TMM-normalized log2-CPM.
# Method: Gene-wise linear regression (lm), fitting one model per gene.
# Multiple testing: Benjamini–Hochberg FDR.

# Model 5 — Reverse regression analysis
      #   Both studies, shared genes, adjusted for study, sex, and age.
      #   Unlike Models 0–4, monocyte percentage is the dependent variable (Y),
      #   while baseline gene expression is the independent variable (X).
      #   This model estimates the change in Day 7 monocyte percentage
      #   associated with a one-unit increase in baseline gene expression.

#   Results will be compared with the voom–limma analyses to assess
#   the consistency of gene-level associations across approaches.


###############################################################################
################################################################################
##############################################################################
################################MODEL0
################################################################################

###############################################################################
#### STEP 1 — PREPARE SAMPLE METADATA
###############################################################################

# Create a dataset containing the harmonized Day 7 monocyte percentages
# and the participant characteristics needed for our regression models.

sample_metadata <- mon_dem_clean %>%
  select(
    participant_id,
    study_accession,
    value_clean,
    age,
    biological_sex,
    vacc_any
  ) %>%
  rename(monocyte_pct = value_clean)


###############################################################################
#### STEP 2 — ALIGN SAMPLE METADATA WITH GENE EXPRESSION MATRIX
###############################################################################

# Reorder participants in sample_metadata to match the columns of logcpm.
# This ensures that each participant's gene expression values correspond
# to their correct Day 7 monocyte percentage.

sample_metadata <- sample_metadata %>%
  slice(match(colnames(logcpm), participant_id))

# Verify that participants are correctly aligned.
stopifnot(
  identical(sample_metadata$participant_id, colnames(logcpm))
)

# Display the dimensions of our metadata and gene expression matrix.
dim(sample_metadata)
dim(logcpm)

###############################################################################
#### STEP 3 — CHECK MISSING VALUES
###############################################################################

#Check missing values in participant characteristics.
# This helps us determine whether all 77 participants can be included
# in Model 0 and which variables are complete for the adjusted models

colSums(is.na(sample_metadata))

# Check whether gene expression data contain missing values.
sum(is.na(logcpm))


###############################################################################
#### STEP 4 — MODEL 0 CREATE THE DESIG MATRIX
###############################################################################

# 4.1 Define the unadjusted regression model.
# Y = Baseline gene expression (Day 0).
# X = Day 7 classical monocyte percentage.
# No adjustment for study, sex, age, or other covariates.

design_0 <- model.matrix(
  ~ monocyte_pct,
  data = sample_metadata
)

# 4.2 Assign participant IDs as row names.
# This ensures alignment with the RNA-seq count matrix.

rownames(design_0) <- sample_metadata$participant_id

# 4.3 Verify that participants are correctly aligned.

stopifnot(
  nrow(design_0) == ncol(dge$counts),
  identical(rownames(design_0), colnames(dge$counts))
)

# 4.4 Inspect the design matrix.

dim(design_0)
head(design_0)


###############################################################################
#### STEP 5 — MODEL 0: APPLY VOOM
###############################################################################

# 5.1 Apply voom to the TMM-normalized RNA-seq data.
# voom transforms RNA-seq counts into log2-CPM expression values
# and estimates precision weights for each gene and participant.
# We use the existing dge object from environment.R,
# which has already been filtered and TMM-normalized.

voom_0 <- voom(
  dge,
  design = design_0,
  plot = TRUE
)

# 5.2 Inspect the resulting voom object.
dim(voom_0$E)
dim(voom_0$weights)

###### VOOM MEAN–VARIANCE TREND
# Gene expression variability generally decreases as expression increases.
# The red curve captures this trend, supporting the use of voom
# precision weights in the subsequent limma analysis.


###############################################################################
#### STEP 6 — MODEL 0: FIT LINEAR MODELS WITH LIMMA
###############################################################################

# 6.1 Fit a weighted linear regression for each gene.
# Y = Baseline gene expression (Day 0).
# X = Day 7 classical monocyte percentage.
#
# The model is:
# Gene expression = beta0 + beta1 * monocyte_pct + error
#
# lmFit() uses the precision weights estimated by voom.

fit_0 <- lmFit(
  voom_0,
  design = design_0
)

# 6.2 Apply empirical Bayes moderation.
# This improves variance estimation by borrowing information
# across genes, particularly useful with high-dimensional RNA-seq data.

fit_0 <- eBayes(fit_0)

# 6.3 Inspect the fitted model.

dim(fit_0$coefficients)
head(fit_0$coefficients)


###############################################################################
######## STEP 7 — MODEL 0: EXTRACT REGRESSION RESULTS
###############################################################################

# 7.1 Extract the regression results for all 17,165 genes.
# We focus on the coefficient associated with Day 7 monocyte percentage.
#
# The results include:
# - logFC: Regression coefficient (change in log2 gene expression
#          per 1-percentage-point increase in monocyte percentage).
# - AveExpr: Average log2 gene expression.
# - t: Moderated t-statistic.
# - P.Value: Raw p-value.
# - adj.P.Val: Benjamini–Hochberg adjusted p-value (FDR).

results_0 <- topTable(
  fit_0,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# 7.2 Inspect the dimensions and first results.

dim(results_0)
head(results_0)

###############################################################################
#### STEP 8 — MODEL 0: COUNT SIGNIFICANT GENES
###############################################################################

# 8.1 Count genes significantly associated with Day 7 monocyte percentage.
# Statistical significance is defined as FDR < 0.05.

sum(results_0$adj.P.Val < 0.05)

# 8.2 Count significant genes with positive regression coefficients.
# Positive coefficient = higher monocyte percentage is associated
# with higher baseline gene expression.

sum(results_0$adj.P.Val < 0.05 & results_0$logFC > 0)

# 8.3 Count significant genes with negative regression coefficients.
# Negative coefficient = higher monocyte percentage is associated
# with lower baseline gene expression.

sum(results_0$adj.P.Val < 0.05 & results_0$logFC < 0)

# 8.4 Calculate the percentage of significant genes.

round(mean(results_0$adj.P.Val < 0.05) * 100, 2)

###############################################################################
#### STEP 9 — MODEL 0: CHECK MONOCYTE DISTRIBUTIONS BY STUDY
###############################################################################

# 9.1 Summarise Day 7 monocyte percentages separately for each study.
# This helps investigate whether study membership could confound
# the associations identified in the unadjusted Model 0.

sample_metadata %>%
  group_by(study_accession) %>%
  summarise(
    n = n(),
    mean_monocytes = mean(monocyte_pct),
    sd_monocytes = sd(monocyte_pct),
    median_monocytes = median(monocyte_pct),
    min_monocytes = min(monocyte_pct),
    max_monocytes = max(monocyte_pct),
    .groups = "drop"
  )

### MODEL 0 — INTERPRETATION
# Of 17,165 shared genes, 11,460 (66.76%) were significantly associated with 
#Day 7 monocyte percentages (FDR < 0.05).
# Associations were balanced between positive (5,657) and negative (5,803) coefficients.
# The high number of significant genes may reflect confounding by study differences, as Model 0 is unadjusted.

###############################################################################
################################################################################
##############################################################################
################################ MODEL1 #####################################
################################################################################

###############################################################################
#### STEP 10 — MODEL 1: CREATE THE ADJUSTED DESIGN MATRIX
###############################################################################

# Covariates: Adjust for differences between studies, sex, age.
# Model:
# Gene expression = beta0 + beta1 * monocyte_pct +
#                   beta2 * study + beta3 * sex +
#                   beta4 * age + error

design_1 <- model.matrix(
  ~ monocyte_pct + study_accession + biological_sex + age,
  data = sample_metadata
)

# 10.2 Assign participant IDs as row names.
# This ensures alignment with the RNA-seq count matrix.

rownames(design_1) <- sample_metadata$participant_id

# 10.3 Verify alignment and check for a full-rank design matrix.

stopifnot(
  nrow(design_1) == ncol(dge$counts),
  identical(rownames(design_1), colnames(dge$counts)),
  qr(design_1)$rank == ncol(design_1)
)

# 10.4 Inspect the adjusted design matrix.

dim(design_1)
head(design_1)
colnames(design_1)

###############################################################################
#### STEP 11 — MODEL 1: APPLY VOOM
###############################################################################

# 11.1 Apply voom using the adjusted design matrix.
# The RNA-seq data have already been filtered and TMM-normalized.
#
# voom estimates the mean–variance relationship and calculates
# precision weights while accounting for the adjusted model design.

voom_1 <- voom(
  dge,
  design = design_1,
  plot = TRUE
)

# 11.2 Verify the dimensions of the expression and weight matrices.

dim(voom_1$E)
dim(voom_1$weights)


###############################################################################
#### STEP 12 — MODEL 1: FIT ADJUSTED LINEAR MODELS USING LIMMA
###############################################################################

# 12.1 Fit a weighted linear regression for each gene.
# Y = Baseline gene expression (Day 0).
# X = Day 7 classical monocyte percentage.
#
# Adjusted for:
# - Study (2024_UGA vs SDY2867).
# - Biological sex (Female vs Male).
# - Age (continuous).
#
# Model:
# Gene expression = beta0 + beta1 * monocyte_pct +
#                   beta2 * study + beta3 * sex +
#                   beta4 * age + error
#
# lmFit() uses the precision weights estimated by voom.

fit_1 <- lmFit(
  voom_1,
  design = design_1
)

# 12.2 Apply empirical Bayes moderation.
# This improves variance estimation by borrowing information across genes.

fit_1 <- eBayes(fit_1)

# 12.3 Inspect the fitted model.

dim(fit_1$coefficients)
head(fit_1$coefficients)


###############################################################################
#### STEP 13 — MODEL 1: EXTRACT ADJUSTED REGRESSION RESULTS
###############################################################################

# 13.1 Extract results for the association between baseline gene
# expression and Day 7 classical monocyte percentage.
#
# Results are adjusted for study, biological sex, and age.
#
# The output includes:
# - logFC: Adjusted regression coefficient.
# - AveExpr: Average log2 gene expression.
# - t: Moderated t-statistic.
# - P.Value: Raw p-value.
# - adj.P.Val: Benjamini–Hochberg adjusted p-value (FDR).
# - B: Log-odds of differential expression under the limma model.

results_1 <- topTable(
  fit_1,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# 13.2 Inspect the dimensions and first results.

dim(results_1)
head(results_1)

###############################################################################
#### STEP 14 — MODEL 1: COUNT SIGNIFICANT GENES
###############################################################################

# 14.1 Count significant genes in Model 1 (FDR < 0.05).

sum(results_1$adj.P.Val < 0.05)

# 14.2 Count significant genes with positive coefficients.

sum(results_1$adj.P.Val < 0.05 & results_1$logFC > 0)

# 14.3 Count significant genes with negative coefficients.

sum(results_1$adj.P.Val < 0.05 & results_1$logFC < 0)

# 14.4 Calculate the percentage of significant genes.

round(mean(results_1$adj.P.Val < 0.05) * 100, 2)

# 14.5 Compare the total number of significant genes between models.

data.frame(
  Model = c("Model 0 - Unadjusted", "Model 1 - Adjusted"),
  Significant_genes = c(
    sum(results_0$adj.P.Val < 0.05),
    sum(results_1$adj.P.Val < 0.05)
  )
)


### MODEL 1 — INTERPRETATION
# After adjusting for study, sex, and age, 1,871 genes (10.90%) remained significantly
#associated with Day 7 monocyte percentages (FDR < 0.05).
# Associations were balanced between positive (987) and negative (884) coefficients.
# Compared with Model 0 (11,460 genes), the 83.7% reduction suggests substantial 
#confounding in the unadjusted analysis.

###############################################################################
#### STEP 15 — SAVE MODEL 0 AND MODEL 1 RESULTS
###############################################################################

# 15.1 Save complete results for Model 0.
# Includes gene IDs, regression coefficients, p-values and FDR.

write.csv(
  results_0,
  file = file.path(output_folder, "model_0_results.csv"),
  row.names = TRUE
)

# 15.2 Save complete results for Model 1.

write.csv(
  results_1,
  file = file.path(output_folder, "model_1_results.csv"),
  row.names = TRUE
)

# 15.3 Verify that both files were created.

file.exists(file.path(output_folder, "model_0_results.csv"))
file.exists(file.path(output_folder, "model_1_results.csv"))




##############################################################################
##################################################################################

###############################################################################
#### STEP 16 — MODELS 0 AND 1: CREATE SUMMARY TABLE
###############################################################################

# 16.1 Create a summary table comparing the unadjusted and adjusted models.
# Statistical significance is defined as FDR < 0.05.

summary_models_0_1 <- data.frame(
  Model = c("Model 0 - Unadjusted", "Model 1 - Adjusted"),
  
  Total_genes = c(
    nrow(results_0),
    nrow(results_1)
  ),
  
  Significant_genes = c(
    sum(results_0$adj.P.Val < 0.05),
    sum(results_1$adj.P.Val < 0.05)
  ),
  
  Positive_associations = c(
    sum(results_0$adj.P.Val < 0.05 & results_0$logFC > 0),
    sum(results_1$adj.P.Val < 0.05 & results_1$logFC > 0)
  ),
  
  Negative_associations = c(
    sum(results_0$adj.P.Val < 0.05 & results_0$logFC < 0),
    sum(results_1$adj.P.Val < 0.05 & results_1$logFC < 0)
  ),
  
  Percentage_significant = c(
    round(mean(results_0$adj.P.Val < 0.05) * 100, 2),
    round(mean(results_1$adj.P.Val < 0.05) * 100, 2)
  )
)

# 16.2 Display the summary table.

print(summary_models_0_1)

# 16.3 Save the summary table in the existing output folder.

write.csv(
  summary_models_0_1,
  file = file.path(output_folder, "summary_models_0_1.csv"),
  row.names = FALSE
)

# 16.4 Verify that the file was saved.

file.exists(file.path(output_folder, "summary_models_0_1.csv"))


###############################################################################
#### STEP 17 — MODELS 0 AND 1: TOP 20 GENES TABLES
###############################################################################

# 17.1 Extract complete results with 95% confidence intervals.
# logFC represents the regression coefficient (beta), not a
# conventional fold change, because monocyte_pct is continuous.

results_0_ci <- topTable(
  fit_0,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none",
  confint = 0.95
)

results_1_ci <- topTable(
  fit_1,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none",
  confint = 0.95
)


# 17.2 Prepare results for ranking.
# Rename logFC to beta and calculate absolute effect size.

prepare_top_genes <- function(results) {
  results %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    rename(beta = logFC) %>%
    mutate(abs_beta = abs(beta))
}

genes_0 <- prepare_top_genes(results_0_ci)
genes_1 <- prepare_top_genes(results_1_ci)


# 17.3 Select Top 20 genes by statistical significance (FDR).

top20_fdr_0 <- genes_0 %>%
  arrange(adj.P.Val, P.Value) %>%
  slice_head(n = 20)

top20_fdr_1 <- genes_1 %>%
  arrange(adj.P.Val, P.Value) %>%
  slice_head(n = 20)


# 17.4 Select Top 20 genes by absolute effect size.
# Only genes with FDR < 0.05 are considered.

top20_effect_0 <- genes_0 %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs_beta)) %>%
  slice_head(n = 20)

top20_effect_1 <- genes_1 %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs_beta)) %>%
  slice_head(n = 20)


# 17.5 Display the four tables.

print(top20_fdr_0)
print(top20_fdr_1)

print(top20_effect_0)
print(top20_effect_1)


# 17.6 Save all four tables in the existing output folder.

write.csv(
  top20_fdr_0,
  file.path(output_folder, "top20_fdr_model_0.csv"),
  row.names = FALSE
)

write.csv(
  top20_fdr_1,
  file.path(output_folder, "top20_fdr_model_1.csv"),
  row.names = FALSE
)

write.csv(
  top20_effect_0,
  file.path(output_folder, "top20_effect_model_0.csv"),
  row.names = FALSE
)

write.csv(
  top20_effect_1,
  file.path(output_folder, "top20_effect_model_1.csv"),
  row.names = FALSE
)


# 17.7 Verify that all four files were saved.

file.exists(file.path(output_folder, c(
  "top20_fdr_model_0.csv",
  "top20_fdr_model_1.csv",
  "top20_effect_model_0.csv",
  "top20_effect_model_1.csv"
)))


###############################################################################
#### STEP 18 — MODELS 0 AND 1: VOLCANO PLOTS
###############################################################################

# 18.1 Prepare results for plotting.
# logFC represents the regression coefficient (beta).
# Genes are classified according to FDR < 0.05 and effect direction.

prepare_volcano <- function(results) {
  
  results %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    mutate(
      beta = logFC,
      neg_log10_p = -log10(pmax(P.Value, .Machine$double.xmin)),
      
      significance = case_when(
        adj.P.Val < 0.05 & beta > 0 ~ "Positive association",
        adj.P.Val < 0.05 & beta < 0 ~ "Negative association",
        TRUE ~ "Not significant"
      )
    )
}

volcano_0 <- prepare_volcano(results_0)
volcano_1 <- prepare_volcano(results_1)


# 18.2 Define a function to generate volcano plots.

create_volcano <- function(data, title) {
  
  ggplot(data, aes(x = beta, y = neg_log10_p, color = significance)) +
    
    geom_point(alpha = 0.6, size = 1) +
    
    scale_color_manual(
      values = c(
        "Positive association" = "#D55E00",
        "Negative association" = "#0072B2",
        "Not significant" = "grey75"
      )
    ) +
    
    labs(
      title = title,
      x = "Regression coefficient (β)",
      y = expression(-log[10](italic(p))),
      color = "Gene classification"
    ) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}


# 18.3 Create volcano plots for both models.

volcano_plot_0 <- create_volcano(
  volcano_0,
  "Model 0 — Unadjusted"
)

volcano_plot_1 <- create_volcano(
  volcano_1,
  "Model 1 — Adjusted for study, sex and age"
)


# 18.4 Display both plots.

print(volcano_plot_0)
print(volcano_plot_1)


# 18.5 Save plots in the existing output folder.

ggsave(
  filename = file.path(output_folder, "volcano_model_0.png"),
  plot = volcano_plot_0,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(output_folder, "volcano_model_1.png"),
  plot = volcano_plot_1,
  width = 8,
  height = 6,
  dpi = 300
)


# 18.6 Verify that both files were saved.

file.exists(file.path(output_folder, c(
  "volcano_model_0.png",
  "volcano_model_1.png"
)))


###############################################################################
#### STEP 19 — MODEL 1: GENE LISTS FOR DAVID (GO ENRICHMENT)
###############################################################################

# 19.1 Prepare Model 1 results.
# Ensembl IDs are stored as row names; we remove any version suffix
# (e.g. ENSG00000123456.12 -> ENSG00000123456) because DAVID doesn't recognise it.

david_input <- results_1 %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  mutate(ensembl_gene_id = sub("\\..*$", "", ensembl_gene_id))


# 19.2 Significant genes split by direction (FDR < 0.05).
# Positive and negative associations are analysed separately because they
# usually reflect different biological processes.

genes_positive <- david_input %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  pull(ensembl_gene_id)

genes_negative <- david_input %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  pull(ensembl_gene_id)

# 19.3 Background: all genes tested in Model 1.
# Only genes expressed in blood and passing filterByExpr could become significant,
# so they (not the whole genome) are the correct reference for enrichment.

genes_background <- david_input$ensembl_gene_id

length(genes_positive)     # expect 987
length(genes_negative)     # expect 884
length(genes_background)   # all tested genes

write_lines(genes_positive,   file.path(output_folder, "david_model1_positive.txt"))
write_lines(genes_negative,   file.path(output_folder, "david_model1_negative.txt"))
write_lines(genes_background, file.path(output_folder, "david_model1_background.txt"))


