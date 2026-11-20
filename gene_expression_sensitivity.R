###############################################################################
## gene_expression_sensitivity.R
#### GENE EXPRESSION SENSITIVITY ANALYSES — MODELS 2, 3 AND 4
###############################################################################

# Model 2: Extended gene set sensitivity analysis
# Model 3: Study-specific analysis — 2024_UGA
# Model 4: Study-specific analysis — SDY2867

# Statistical approach:
# RNA-seq analysis using voom-limma.
# Study-specific models adjusted for biological sex and age.
# Benjamini-Hochberg correction for multiple testing (FDR < 0.05).




###############################################################################
###############################################################################
################################  MODEL 2
##############################################################################
###############################################################################
#### STEP 1 — LOAD ENVIRONMENT AND CHECK REQUIRED OBJECTS
###############################################################################

# 1.1 Load the existing environment script.
# This contains our data preparation and normalization pipeline.

source("environment.R")


# 1.2 Verify that the required objects are available.

required_objects <- c(
  "mon_dem_clean",
  "sample_info",
  "rna_mon",
  "common_genes",
  "dge",
  "logcpm",
  "output_folder"
)

sapply(required_objects, exists)


# 1.3 Check the number of participants in each study.

sample_info %>%
  count(study_accession)


# 1.4 Check the number of genes available in each study.
# This will help us understand which genes were excluded
# from the combined analysis.

rna_mon %>%
  group_by(study_accession) %>%
  summarise(
    n_genes = n_distinct(ensembl_gene_id),
    n_participants = n_distinct(participant_id),
    .groups = "drop"
  )


###############################################################################
#### STEP 2 — EXAMINE GENE OVERLAP BETWEEN STUDIES
###############################################################################

# 2.1 Extract unique gene identifiers from each study.

genes_uga <- rna_mon %>%
  filter(study_accession == "2024_UGA") %>%
  pull(ensembl_gene_id) %>%
  unique()

genes_sdy <- rna_mon %>%
  filter(study_accession == "SDY2867") %>%
  pull(ensembl_gene_id) %>%
  unique()


# 2.2 Identify shared and study-specific genes.

genes_shared <- intersect(genes_uga, genes_sdy)

genes_uga_only <- setdiff(genes_uga, genes_sdy)

genes_sdy_only <- setdiff(genes_sdy, genes_uga)


# 2.3 Summarise gene availability.

gene_overlap_summary <- data.frame(
  Category = c(
    "Genes in UGA",
    "Genes in SDY2867",
    "Shared genes",
    "UGA-only genes",
    "SDY2867-only genes"
  ),
  N_genes = c(
    length(genes_uga),
    length(genes_sdy),
    length(genes_shared),
    length(genes_uga_only),
    length(genes_sdy_only)
  )
)

print(gene_overlap_summary)


# 2.4 Compare shared genes with the existing common_genes object.

length(common_genes)
length(genes_shared)


###############################################################################
#### STEP 3 — EXAMINE EXPRESSION OF UGA-ONLY GENES
###############################################################################

# 3.1 Select genes found only in UGA.

rna_uga_only <- rna_mon %>%
  filter(
    study_accession == "2024_UGA",
    ensembl_gene_id %in% genes_uga_only
  )


# 3.2 Count genes and participants.

rna_uga_only %>%
  summarise(
    n_genes = n_distinct(ensembl_gene_id),
    n_participants = n_distinct(participant_id),
    n_observations = n()
  )


# 3.3 Examine the distribution of RNA-seq counts.

summary(rna_uga_only$counts)


# 3.4 Calculate the percentage of observations with zero counts.

mean(rna_uga_only$counts == 0, na.rm = TRUE) * 100


# 3.5 Count genes with at least one non-zero observation.

rna_uga_only %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    total_counts = sum(counts, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(
    genes_with_expression = sum(total_counts > 0),
    genes_without_expression = sum(total_counts == 0)
  )


###############################################################################
#### STEP 4 — CHECK UGA RNA-SEQ DATA STRUCTURE
###############################################################################

# 4.1 Select all UGA RNA-seq observations.

rna_uga <- rna_mon %>%
  filter(study_accession == "2024_UGA")


# 4.2 Check for duplicate gene-participant combinations.
# Each gene should have one count per participant.

duplicate_uga <- rna_uga %>%
  count(participant_id, ensembl_gene_id) %>%
  filter(n > 1)

nrow(duplicate_uga)


# 4.3 Check whether counts contain missing or negative values.

sum(is.na(rna_uga$counts))

sum(rna_uga$counts < 0, na.rm = TRUE)


# 4.4 Check the number of genes measured per participant.

rna_uga %>%
  group_by(participant_id) %>%
  summarise(
    n_genes = n_distinct(ensembl_gene_id),
    .groups = "drop"
  ) %>%
  summarise(
    min_genes = min(n_genes),
    max_genes = max(n_genes)
  )

###############################################################################
#### STEP 5 — CONSTRUCT UGA EXPRESSION MATRIX
###############################################################################

# 5.1 Convert RNA-seq data from long to wide format.
# Rows = genes; columns = participants; values = raw counts.

counts_uga_df <- rna_uga %>%
  select(ensembl_gene_id, participant_id, counts) %>%
  tidyr::pivot_wider(
    names_from = participant_id,
    values_from = counts
  )


# 5.2 Convert the data frame into a numeric matrix.

counts_uga <- counts_uga_df %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()


# 5.3 Check matrix dimensions.

dim(counts_uga)


# 5.4 Verify that all values are valid RNA-seq counts.

is.numeric(counts_uga)

sum(is.na(counts_uga))

sum(counts_uga < 0)


# 5.5 Check that participant IDs match our metadata.

all(colnames(counts_uga) %in% sample_info$participant_id)


# 5.6 Display the first five genes and participants.

counts_uga[1:5, 1:5]



###############################################################################
#### STEP 6 — PREPARE UGA METADATA AND DESIGN MATRIX
###############################################################################

# 6.1 Extract demographic and monocyte information for UGA participants.

metadata_uga <- mon_dem_clean %>%
  filter(study_accession == "2024_UGA") %>%
  select(
    participant_id,
    study_accession,
    value_clean,
    age,
    biological_sex
  ) %>%
  rename(monocyte_pct = value_clean)


# 6.2 Arrange metadata in exactly the same order as the count matrix.

metadata_uga <- metadata_uga %>%
  slice(match(colnames(counts_uga), participant_id))


# 6.3 Verify participant alignment and missing values.

stopifnot(
  identical(metadata_uga$participant_id, colnames(counts_uga))
)

nrow(metadata_uga)

colSums(is.na(metadata_uga))


# 6.4 Check the distribution of biological sex.

table(metadata_uga$biological_sex)


# 6.5 Create the adjusted design matrix.
# Study is excluded because all participants belong to UGA.

design_uga <- model.matrix(
  ~ monocyte_pct + biological_sex + age,
  data = metadata_uga
)

rownames(design_uga) <- metadata_uga$participant_id


# 6.6 Check design matrix dimensions and coefficients.

dim(design_uga)

colnames(design_uga)

qr(design_uga)$rank

ncol(design_uga)


###############################################################################
#### STEP 7 — UGA EXPRESSION FILTERING AND TMM NORMALIZATION
###############################################################################

# NOTE: TMM normalization was already performed in environment.R for Models 0–1.
# We recalculate it here because UGA has a different sample and gene set.


# 7.1 Create an edgeR DGEList object using raw RNA-seq counts.

dge_uga <- edgeR::DGEList(counts = counts_uga)


# 7.2 Identify sufficiently expressed genes using the model design.

keep_uga <- edgeR::filterByExpr(
  dge_uga,
  design = design_uga
)


# 7.3 Summarise filtering results.

filter_summary_uga <- data.frame(
  Category = c(
    "Genes before filtering",
    "Genes retained",
    "Genes excluded"
  ),
  N_genes = c(
    nrow(dge_uga),
    sum(keep_uga),
    sum(!keep_uga)
  )
)

print(filter_summary_uga)


# 7.4 Apply expression filtering.

dge_uga <- dge_uga[keep_uga, , keep.lib.sizes = FALSE]


# 7.5 Calculate UGA-specific TMM normalization factors.

dge_uga <- edgeR::calcNormFactors(
  dge_uga,
  method = "TMM"
)


# 7.6 Check final dimensions.

dim(dge_uga)


# 7.7 Inspect normalization factors.

dge_uga$samples


###############################################################################
#### STEP 8 — IDENTIFY RETAINED UGA-ONLY GENES
###############################################################################

# 8.1 Extract genes retained after UGA-specific expression filtering.

genes_uga_filtered <- rownames(dge_uga)


# 8.2 Identify retained genes absent from SDY2867.

genes_uga_only_filtered <- intersect(
  genes_uga_filtered,
  genes_uga_only
)


# 8.3 Identify retained genes shared between both studies.

genes_shared_uga_filtered <- intersect(
  genes_uga_filtered,
  genes_shared
)


# 8.4 Summarise retained genes by category.

gene_summary_uga <- data.frame(
  Category = c(
    "Total genes retained in UGA",
    "Retained shared genes",
    "Retained UGA-only genes"
  ),
  N_genes = c(
    length(genes_uga_filtered),
    length(genes_shared_uga_filtered),
    length(genes_uga_only_filtered)
  )
)

print(gene_summary_uga)


# 8.5 Verify that the categories sum to the total.

stopifnot(
  length(genes_shared_uga_filtered) +
    length(genes_uga_only_filtered) ==
    nrow(dge_uga)
)

###############################################################################
#### STEP 9 — MODEL 2: UGA-ONLY GENES (VOOM–LIMMA)
###############################################################################

# Model 2 investigates the 1,286 UGA-only genes excluded from Models 0–1.
# We use UGA-specific TMM normalization and adjust for sex and age.


# 9.1 Extract UGA-only genes from the normalized UGA dataset.

dge_model2 <- dge_uga[genes_uga_only_filtered, , keep.lib.sizes = TRUE]

dim(dge_model2)


# 9.2 Apply voom transformation.
# Keep the library sizes and TMM factors estimated using all eligible UGA genes.

voom_model2 <- limma::voom(
  dge_model2,
  design = design_uga,
  plot = TRUE
)


# 9.3 Fit gene-wise linear models.

fit_model2 <- limma::lmFit(
  voom_model2,
  design = design_uga
)


# 9.4 Apply empirical Bayes moderation.

fit_model2 <- limma::eBayes(fit_model2)


# 9.5 Extract results for Day 7 monocyte percentage.
# logFC represents the regression coefficient (beta).

results_model2 <- limma::topTable(
  fit_model2,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# 9.6 Summarise statistically significant genes (FDR < 0.05).

model_2_summary <- data.frame(
  Model = "Model 2 - UGA-only genes",
  Total_genes = nrow(results_model2),
  Significant_genes = sum(results_model2$adj.P.Val < 0.05),
  Positive_associations = sum(
    results_model2$adj.P.Val < 0.05 & results_model2$logFC > 0
  ),
  Negative_associations = sum(
    results_model2$adj.P.Val < 0.05 & results_model2$logFC < 0
  )
)

print(model_2_summary)


# 9.7 Display the 10 most statistically significant genes.

head(results_model2, 10)


###############################################################################
#### STEP 10 — SAVE MODEL 2 RESULTS
###############################################################################

# 10.1 Save complete gene-level results.

write.csv(
  results_model2,
  file = file.path(output_folder, "model2_uga_only_genes_results.csv"),
  row.names = TRUE
)


# 10.2 Save summary statistics.

write.csv(
  model_2_summary,
  file = file.path(output_folder, "model2_uga_only_genes_summary.csv"),
  row.names = FALSE
)


# 10.3 Save the voom diagnostic plot.

png(
  filename = file.path(output_folder, "model2_uga_only_genes_voom.png"),
  width = 1200,
  height = 900,
  res = 150
)

invisible(limma::voom(
  dge_model2,
  design = design_uga,
  plot = TRUE
))

dev.off()

### MODEL 2 INTERPRETATION
### Among 1,286 UGA-only genes, none were significantly associated with Day 7 
# classical monocyte percentage after adjusting for age, sex, and BH-FDR (FDR < 0.05).
### The lowest unadjusted p-value was 0.000203, but the corresponding adjusted 
# p-value was 0.168.
### These findings provide no FDR-significant evidence of association, 
# although the small sample size (n = 24) may limit statistical power.




###############################################################################
###############################################################################
############################### MODEL 3
##############################################################################

###############################################################################
#### STEP 1 — MODEL 3: UGA-SPECIFIC VOOM–LIMMA ANALYSIS
###############################################################################

# Model 3 investigates all 19,126 eligible UGA genes (n = 24).
# Outcome: baseline gene expression (Day 0).
# Predictor: Day 7 classical monocyte percentage.
# Covariates: biological sex and age.


# 1.1 Apply voom transformation.

voom_model3 <- limma::voom(
  dge_uga,
  design = design_uga,
  plot = TRUE
)


# 1.2 Fit gene-wise linear models.

fit_model3 <- limma::lmFit(
  voom_model3,
  design = design_uga
)


# 1.3 Apply empirical Bayes moderation.

fit_model3 <- limma::eBayes(fit_model3)


# 1.4 Extract results for Day 7 monocyte percentage.
# logFC represents the regression coefficient (beta).

results_model3 <- limma::topTable(
  fit_model3,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# 1.5 Summarise significant genes (FDR < 0.05).

model_3_summary <- data.frame(
  Model = "Model 3 - UGA-specific",
  Total_genes = nrow(results_model3),
  Significant_genes = sum(results_model3$adj.P.Val < 0.05),
  Positive_associations = sum(
    results_model3$adj.P.Val < 0.05 & results_model3$logFC > 0
  ),
  Negative_associations = sum(
    results_model3$adj.P.Val < 0.05 & results_model3$logFC < 0
  )
)

print(model_3_summary)


# 1.6 Display the 10 most statistically significant genes.

head(results_model3, 10)


###############################################################################
#### STEP 2 — SAVE MODEL 3 RESULTS
###############################################################################

# 2.1 Save complete gene-level results.

write.csv(
  results_model3,
  file = file.path(output_folder, "model3_uga_results.csv"),
  row.names = TRUE
)


# 2.2 Save summary statistics.

write.csv(
  model_3_summary,
  file = file.path(output_folder, "model3_uga_summary.csv"),
  row.names = FALSE
)


# 2.3 Save the voom diagnostic plot.

png(
  filename = file.path(output_folder, "model3_uga_voom.png"),
  width = 1200,
  height = 900,
  res = 150
)

invisible(limma::voom(
  dge_uga,
  design = design_uga,
  plot = TRUE
))

dev.off()


### MODEL 3 INTERPRETATION
### Among 19,126 UGA genes, none were significantly associated with Day 
# 7 classical monocyte percentage after adjusting for age, sex, and BH-FDR (FDR < 0.05).
### The lowest unadjusted p-value was 5.76e-06, but the corresponding adjusted p-value
# was 0.0655.
### Several genes showed suggestive associations, but the small sample size (n = 24)
# and multiple-testing burden may limit statistical power.

###############################################################################
#### STEP 3 — COMPARE MODEL 1 AND MODEL 3
###############################################################################

# 3.1 Extract genes significant in Model 1 (FDR < 0.05).

sig_genes_model1 <- rownames(results_1)[
  results_1$adj.P.Val < 0.05
]


# 3.2 Identify which significant genes are also analysed in Model 3.

genes_to_compare <- intersect(
  sig_genes_model1,
  rownames(results_model3)
)


# 3.3 Create a comparison table of regression coefficients and p-values.

comparison_1_3 <- data.frame(
  gene_id = genes_to_compare,
  beta_model1 = results_1[genes_to_compare, "logFC"],
  beta_model3 = results_model3[genes_to_compare, "logFC"],
  p_model1 = results_1[genes_to_compare, "P.Value"],
  p_model3 = results_model3[genes_to_compare, "P.Value"],
  FDR_model1 = results_1[genes_to_compare, "adj.P.Val"],
  FDR_model3 = results_model3[genes_to_compare, "adj.P.Val"]
)


# 3.4 Check whether associations have the same direction.

comparison_1_3$same_direction <- with(
  comparison_1_3,
  sign(beta_model1) == sign(beta_model3)
)


# 3.5 Summarise agreement between models.

comparison_summary <- data.frame(
  Significant_genes_model1 = length(sig_genes_model1),
  Genes_available_model3 = nrow(comparison_1_3),
  Same_direction = sum(comparison_1_3$same_direction),
  Same_direction_pct = round(
    mean(comparison_1_3$same_direction) * 100, 2
  ),
  Nominally_significant_model3 = sum(
    comparison_1_3$p_model3 < 0.05
  ),
  FDR_significant_model3 = sum(
    comparison_1_3$FDR_model3 < 0.05
  ),
  Beta_correlation = round(
    cor(comparison_1_3$beta_model1,
        comparison_1_3$beta_model3),
    3
  )
)

print(comparison_summary)


# 3.6 Display the first 10 genes.

head(comparison_1_3, 10)


### MODEL 1 VS MODEL 3 — INTERPRETATION
### Of 1,871 FDR-significant genes in Model 1, 1,802 were available in Model 3.
### Most associations had the same direction (93.84%), with strongly correlated 
# coefficients (r = 0.821).
### Although 343 genes had nominal p < 0.05 in UGA, none survived FDR correction.
### These findings are compatible with reduced statistical power in UGA-only (n = 24),
#but cohort differences cannot be excluded.





###############################################################################
###############################################################################
################################# MODEL 4
###############################################################################

###############################################################################
#### STEP 1 — MODEL 4: PREPARE SDY2867 DATA
###############################################################################

# Model 4 investigates baseline gene expression in SDY2867 participants.
# We will use all eligible SDY2867 genes, adjusting for age and sex.


# 1.1 Extract SDY2867 RNA-seq data.

rna_sdy <- rna_mon %>%
  filter(study_accession == "SDY2867")


# 1.2 Check the number of participants and genes.

cat("Participants:", n_distinct(rna_sdy$participant_id), "\n")
cat("Genes:", n_distinct(rna_sdy$ensembl_gene_id), "\n")


# 1.3 Check for duplicated participant–gene combinations.

sum(duplicated(
  rna_sdy[, c("participant_id", "ensembl_gene_id")]
))


# 1.4 Check for missing or negative RNA-seq counts.

sum(is.na(rna_sdy$counts))
sum(rna_sdy$counts < 0, na.rm = TRUE)


# 1.5 Check the number of genes measured per participant.

rna_sdy %>%
  count(participant_id, name = "n_genes") %>%
  summarise(
    min_genes = min(n_genes),
    max_genes = max(n_genes)
  )


###############################################################################
#### STEP 2 — MODEL 4: CREATE SDY2867 COUNT MATRIX
###############################################################################

# 2.1 Convert RNA-seq data into a wide-format table.

counts_sdy_df <- rna_sdy %>%
  select(ensembl_gene_id, participant_id, counts) %>%
  tidyr::pivot_wider(
    names_from = participant_id,
    values_from = counts
  )


# 2.2 Convert the table into a numeric count matrix.

counts_sdy <- counts_sdy_df %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()


# 2.3 Check matrix dimensions.

dim(counts_sdy)


# 2.4 Verify that counts are numeric.

is.numeric(counts_sdy)


# 2.5 Verify that all participants match sample_info.

all(colnames(counts_sdy) %in% sample_info$participant_id)


###############################################################################
#### STEP 3 — MODEL 4: PREPARE METADATA AND DESIGN MATRIX
###############################################################################

# 3.1 Extract SDY2867 participant metadata.

metadata_sdy <- mon_dem_clean %>%
  filter(study_accession == "SDY2867") %>%
  select(participant_id, study_accession, value_clean, age, biological_sex) %>%
  rename(monocyte_pct = value_clean)


# 3.2 Align metadata with the columns of the count matrix.

metadata_sdy <- metadata_sdy %>%
  slice(match(colnames(counts_sdy), participant_id))


# 3.3 Verify participant alignment and missing values.

stopifnot(identical(metadata_sdy$participant_id, colnames(counts_sdy)))

colSums(is.na(metadata_sdy))


# 3.4 Inspect sample characteristics.

summary(metadata_sdy[, c("monocyte_pct", "age")])

table(metadata_sdy$biological_sex)


# 3.5 Create the adjusted design matrix.

design_sdy <- model.matrix(
  ~ monocyte_pct + biological_sex + age,
  data = metadata_sdy
)

rownames(design_sdy) <- metadata_sdy$participant_id


# 3.6 Check dimensions, coefficients, and matrix rank.

dim(design_sdy)

colnames(design_sdy)

qr(design_sdy)$rank


###############################################################################
#### STEP 4 — MODEL 4: EXPRESSION FILTERING AND TMM NORMALIZATION
###############################################################################

# NOTE: TMM normalization was already performed in environment.R for Models 0–1.
# We recalculate it because SDY2867 has a different sample and gene set.


# 4.1 Create an edgeR DGEList object using raw RNA-seq counts.

dge_sdy <- edgeR::DGEList(counts = counts_sdy)


# 4.2 Identify sufficiently expressed genes.

keep_sdy <- edgeR::filterByExpr(
  dge_sdy,
  design = design_sdy
)


# 4.3 Summarise filtering results.

filter_summary_sdy <- data.frame(
  Category = c(
    "Genes before filtering",
    "Genes retained",
    "Genes excluded"
  ),
  N_genes = c(
    nrow(dge_sdy),
    sum(keep_sdy),
    sum(!keep_sdy)
  )
)

print(filter_summary_sdy)


# 4.4 Apply expression filtering.

dge_sdy <- dge_sdy[keep_sdy, , keep.lib.sizes = FALSE]


# 4.5 Calculate SDY2867-specific TMM normalization factors.

dge_sdy <- edgeR::calcNormFactors(
  dge_sdy,
  method = "TMM"
)


# 4.6 Check final dimensions.

dim(dge_sdy)


# 4.7 Inspect normalization factors.

dge_sdy$samples


###############################################################################
#### STEP 5 — MODEL 4: SDY2867-SPECIFIC VOOM–LIMMA ANALYSIS
###############################################################################

# Model 4 investigates all 17,973 eligible SDY2867 genes (n = 53).
# Outcome: baseline gene expression (Day 0).
# Predictor: Day 7 classical monocyte percentage.
# Covariates: biological sex and age.


# 5.1 Apply voom transformation.

voom_model4 <- limma::voom(
  dge_sdy,
  design = design_sdy,
  plot = TRUE
)


# 5.2 Fit gene-wise linear models.

fit_model4 <- limma::lmFit(
  voom_model4,
  design = design_sdy
)


# 5.3 Apply empirical Bayes moderation.

fit_model4 <- limma::eBayes(fit_model4)


# 5.4 Extract results for Day 7 monocyte percentage.
# logFC represents the regression coefficient (beta).

results_model4 <- limma::topTable(
  fit_model4,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# 5.5 Summarise significant genes (FDR < 0.05).

model_4_summary <- data.frame(
  Model = "Model 4 - SDY2867-specific",
  Total_genes = nrow(results_model4),
  Significant_genes = sum(results_model4$adj.P.Val < 0.05),
  Positive_associations = sum(
    results_model4$adj.P.Val < 0.05 & results_model4$logFC > 0
  ),
  Negative_associations = sum(
    results_model4$adj.P.Val < 0.05 & results_model4$logFC < 0
  )
)

print(model_4_summary)


# 5.6 Display the 10 most statistically significant genes.

head(results_model4, 10)


#############################################################################
### MODEL 4 INTERPRETATION
### Among 17,973 SDY2867 genes, 3,046 were significantly associated with Day 7
# classical monocyte percentage after adjusting for age, sex, and BH-FDR (FDR < 0.05).
### Of these, 1,351 showed positive associations and 1,695 showed negative associations.
### Unlike UGA (Model 3), SDY2867 identified numerous significant genes, potentially 
# reflecting greater statistical power and/or cohort-specific differences.


###############################################################################
#### STEP 6 — COMPARE MODEL 3 (UGA) AND MODEL 4 (SDY2867)
###############################################################################

# 6.1 Identify genes analysed in both models.

genes_common_3_4 <- intersect(
  rownames(results_model3),
  rownames(results_model4)
)


# 6.2 Create a comparison table.

comparison_3_4 <- data.frame(
  gene_id = genes_common_3_4,
  beta_uga = results_model3[genes_common_3_4, "logFC"],
  beta_sdy = results_model4[genes_common_3_4, "logFC"],
  p_uga = results_model3[genes_common_3_4, "P.Value"],
  p_sdy = results_model4[genes_common_3_4, "P.Value"],
  FDR_uga = results_model3[genes_common_3_4, "adj.P.Val"],
  FDR_sdy = results_model4[genes_common_3_4, "adj.P.Val"]
)


# 6.3 Check agreement in association direction.

comparison_3_4$same_direction <- with(
  comparison_3_4,
  sign(beta_uga) == sign(beta_sdy)
)


# 6.4 Summarise agreement between cohorts.

comparison_summary_3_4 <- data.frame(
  Genes_compared = nrow(comparison_3_4),
  Same_direction = sum(comparison_3_4$same_direction),
  Same_direction_pct = round(
    mean(comparison_3_4$same_direction) * 100, 2
  ),
  Beta_correlation = round(
    cor(comparison_3_4$beta_uga, comparison_3_4$beta_sdy),
    3
  ),
  Significant_UGA = sum(comparison_3_4$FDR_uga < 0.05),
  Significant_SDY = sum(comparison_3_4$FDR_sdy < 0.05),
  Significant_both = sum(
    comparison_3_4$FDR_uga < 0.05 &
      comparison_3_4$FDR_sdy < 0.05
  )
)

print(comparison_summary_3_4)

################################################################################
### MODEL 3 VS MODEL 4 — INTERPRETATION
### Among 16,747 genes analysed in both cohorts, 57.3% showed associations in the 
# same direction.
### Regression coefficients were weakly correlated (r = 0.161), suggesting limited 
# consistency between cohorts.
### SDY2867 identified 2,972 FDR-significant shared genes, while UGA identified none.
### These differences may reflect statistical power, cohort heterogeneity, or both.


###############################################################################
#### STEP 7 — COMPARE MODEL 1 (COMBINED) AND MODEL 4 (SDY2867)
###############################################################################

# 7.1 Extract genes significant in Model 1 (FDR < 0.05).

sig_genes_model1 <- rownames(results_1)[
  results_1$adj.P.Val < 0.05
]


# 7.2 Identify significant Model 1 genes also analysed in Model 4.

genes_to_compare_1_4 <- intersect(
  sig_genes_model1,
  rownames(results_model4)
)


# 7.3 Create a comparison table.

comparison_1_4 <- data.frame(
  gene_id = genes_to_compare_1_4,
  beta_model1 = results_1[genes_to_compare_1_4, "logFC"],
  beta_model4 = results_model4[genes_to_compare_1_4, "logFC"],
  p_model1 = results_1[genes_to_compare_1_4, "P.Value"],
  p_model4 = results_model4[genes_to_compare_1_4, "P.Value"],
  FDR_model1 = results_1[genes_to_compare_1_4, "adj.P.Val"],
  FDR_model4 = results_model4[genes_to_compare_1_4, "adj.P.Val"]
)


# 7.4 Check agreement in association direction.

comparison_1_4$same_direction <- with(
  comparison_1_4,
  sign(beta_model1) == sign(beta_model4)
)


# 7.5 Summarise agreement between models.

comparison_summary_1_4 <- data.frame(
  Significant_genes_model1 = length(sig_genes_model1),
  Genes_available_model4 = nrow(comparison_1_4),
  Same_direction = sum(comparison_1_4$same_direction),
  Same_direction_pct = round(
    mean(comparison_1_4$same_direction) * 100, 2
  ),
  Nominally_significant_model4 = sum(
    comparison_1_4$p_model4 < 0.05
  ),
  FDR_significant_model4 = sum(
    comparison_1_4$FDR_model4 < 0.05
  ),
  Beta_correlation = round(
    cor(
      comparison_1_4$beta_model1,
      comparison_1_4$beta_model4
    ),
    3
  )
)

print(comparison_summary_1_4)


# 7.6 Display the first 10 genes.

head(comparison_1_4, 10)

#################################################################################
### MODEL 1 VS MODEL 4 — INTERPRETATION
### Of 1,871 FDR-significant genes in Model 1, 1,846 were available in Model 4.
### Associations showed 99.95% directional agreement and strongly correlated 
# coefficients (r = 0.977).
### Among these genes, 1,527 remained FDR-significant in SDY2867.
### Combined-model associations closely resemble SDY2867 findings, but this is
# not independent replication.


###############################################################################
#### STEP 8 — STUDY × MONOCYTE INTERACTION ANALYSIS
###############################################################################

# Test whether gene–monocyte associations differ between UGA and SDY2867.
# Outcome: baseline gene expression (Day 0).
# Predictor: Day 7 classical monocyte percentage.
# Covariates: study, biological sex, age.
# Interaction: monocyte percentage × study.


# 8.1 Create the design matrix including the interaction.

design_interaction <- model.matrix(
  ~ monocyte_pct * study_accession + biological_sex + age,
  data = sample_metadata
)

rownames(design_interaction) <- sample_metadata$participant_id

# Check coefficients and matrix rank.

colnames(design_interaction)
dim(design_interaction)
qr(design_interaction)$rank


# 8.2 Apply voom transformation.

voom_interaction <- limma::voom(
  dge,
  design = design_interaction,
  plot = TRUE
)


# 8.3 Fit gene-wise linear models.

fit_interaction <- limma::lmFit(
  voom_interaction,
  design = design_interaction
)


# 8.4 Apply empirical Bayes moderation.

fit_interaction <- limma::eBayes(fit_interaction)


# 8.5 Identify the interaction coefficient.

interaction_coef <- grep(
  "monocyte_pct:study_accession",
  colnames(design_interaction),
  value = TRUE
)

print(interaction_coef)

stopifnot(length(interaction_coef) == 1)


# 8.6 Extract interaction results with BH-FDR correction.
# logFC represents the difference in slopes between studies.

results_interaction <- limma::topTable(
  fit_interaction,
  coef = interaction_coef,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


# 8.7 Summarise statistically significant interactions.

interaction_summary <- data.frame(
  Total_genes = nrow(results_interaction),
  Significant_interactions = sum(
    results_interaction$adj.P.Val < 0.05
  ),
  Positive_interactions = sum(
    results_interaction$adj.P.Val < 0.05 &
      results_interaction$logFC > 0
  ),
  Negative_interactions = sum(
    results_interaction$adj.P.Val < 0.05 &
      results_interaction$logFC < 0
  )
)

print(interaction_summary)


# 8.8 Display the 10 most significant interactions.

head(results_interaction, 10)


### MODEL 8 — STUDY × MONOCYTE INTERACTION ANALYSIS
### Among 17,165 genes, 15 showed significant study-by-monocyte interactions
# (BH-FDR < 0.05).
### Of these, 4 had positive interaction coefficients and 11 had negative coefficients.
### These results indicate that the association between baseline gene expression
# and Day 7 monocyte percentage differs significantly between cohorts for a small
# subset of genes.
### For the remaining genes, no statistically significant interaction was detected.
### However, nonsignificant interactions do not establish equivalence between cohorts,
# particularly given the smaller UGA sample size.

###############################################################################
###############################################################################
################################# MODEL 5
###############################################################################

###############################################################################
#### STEP 1 — MODEL 5: DATA PREPARATION
###############################################################################

# Model 5:
# Outcome: Day 7 classical monocyte percentage.
# Predictor: Day 0 gene expression (one gene at a time).
# Covariates: study, biological sex, age.
# Method: gene-wise linear regression.


# 1.1 Check dimensions of the expression matrix.

dim(logcpm)


# 1.2 Check dimensions of the participant metadata.

dim(sample_metadata)


# 1.3 Check that participants are correctly aligned.

identical(
  colnames(logcpm),
  sample_metadata$participant_id
)


# 1.4 Check for missing values in the outcome and covariates.

colSums(
  is.na(
    sample_metadata[
      , c(
        "monocyte_pct",
        "study_accession",
        "biological_sex",
        "age"
      )
    ]
  )
)


# 1.5 Check for missing or non-finite expression values.

sum(!is.finite(logcpm))


###############################################################################
#### STEP 2 — MODEL 5: GENE-WISE LINEAR REGRESSION
###############################################################################

# 2.1 Create the covariate design matrix.

design_model5 <- model.matrix(
  ~ study_accession + biological_sex + age,
  data = sample_metadata
)


# 2.2 Extract the outcome variable.

y_model5 <- sample_metadata$monocyte_pct


# 2.3 Fit one adjusted linear regression per gene.
# Gene expression is the predictor; monocyte percentage is the outcome.

model5_gene_results <- lapply(seq_len(nrow(logcpm)), function(i) {
  
  gene_expression <- as.numeric(logcpm[i, ])
  
  model <- lm.fit(
    x = cbind(design_model5, gene_expression),
    y = y_model5
  )
  
  # Calculate coefficient standard errors and p-values.
  
  df_residual <- model$df.residual
  
  residual_variance <- sum(model$residuals^2) / df_residual
  
  coefficient_variance <- residual_variance *
    diag(chol2inv(qr.R(model$qr)))
  
  beta <- tail(model$coefficients, 1)
  
  se <- sqrt(tail(coefficient_variance, 1))
  
  t_value <- beta / se
  
  p_value <- 2 * pt(
    abs(t_value),
    df = df_residual,
    lower.tail = FALSE
  )
  
  c(
    beta = beta,
    SE = se,
    t = t_value,
    P.Value = p_value
  )
})


# 2.4 Combine results into a table.

results_model5 <- as.data.frame(
  do.call(rbind, model5_gene_results)
)

rownames(results_model5) <- rownames(logcpm)


# 2.5 Apply BH-FDR correction.

results_model5$adj.P.Val <- p.adjust(
  results_model5$P.Value,
  method = "BH"
)


# 2.6 Sort genes by statistical significance.

results_model5 <- results_model5[
  order(results_model5$P.Value),
]


# 2.7 Summarise results.

model5_summary <- data.frame(
  Total_genes = nrow(results_model5),
  
  Significant_genes = sum(
    results_model5$adj.P.Val < 0.05
  ),
  
  Positive_associations = sum(
    results_model5$adj.P.Val < 0.05 &
      results_model5$beta > 0
  ),
  
  Negative_associations = sum(
    results_model5$adj.P.Val < 0.05 &
      results_model5$beta < 0
  )
)

print(model5_summary)


# 2.8 Display the 10 most significant genes.

head(results_model5, 10)


### MODEL 5 — INTERPRETATION
### Gene-wise linear regression identified 1,660 genes significantly associated
### with Day 7 classical monocyte percentage after adjustment for study, sex,
### and age (BH-FDR < 0.05).
### Among these, 800 showed positive associations and 860 negative associations.
### Results were broadly similar in the number of significant genes to Model 1,
### although the two models use different outcome definitions and methods.

###############################################################################
#### STEP 3 — MODEL 5: COMPARE WITH MODEL 1
###############################################################################

# 3.1 Clean up coefficient column names.

colnames(results_model5) <- c(
  "beta",
  "SE",
  "t",
  "P.Value",
  "adj.P.Val"
)


# 3.2 Identify significant genes in each model.

sig_model1 <- rownames(results_1)[
  results_1$adj.P.Val < 0.05
]

sig_model5 <- rownames(results_model5)[
  results_model5$adj.P.Val < 0.05
]


# 3.3 Identify genes significant in both models.

sig_both_1_5 <- intersect(
  sig_model1,
  sig_model5
)


# 3.4 Calculate overlap and directional agreement.

comparison_summary_1_5 <- data.frame(
  Significant_model1 = length(sig_model1),
  Significant_model5 = length(sig_model5),
  Significant_both = length(sig_both_1_5),
  Model1_recovered_pct = round(
    length(sig_both_1_5) / length(sig_model1) * 100,
    2
  ),
  Same_direction = sum(
    sign(results_1[sig_both_1_5, "logFC"]) ==
      sign(results_model5[sig_both_1_5, "beta"])
  )
)


# 3.5 Display results.

print(comparison_summary_1_5)

###############################################################################
#### STEP 4 — MODEL 5: SAVE RESULTS
###############################################################################

# 4.1 Save complete Model 5 results.

write.csv(
  results_model5,
  file = file.path(output_folder, "model5_gene_wise_results.csv"),
  row.names = TRUE
)


# 4.2 Save significant genes only.

write.csv(
  results_model5[
    results_model5$adj.P.Val < 0.05,
  ],
  file = file.path(output_folder, "model5_significant_genes.csv"),
  row.names = TRUE
)


# 4.3 Save Model 5 summary.

write.csv(
  model5_summary,
  file = file.path(output_folder, "model5_summary.csv"),
  row.names = FALSE
)


# 4.4 Save comparison with Model 1.

write.csv(
  comparison_summary_1_5,
  file = file.path(output_folder, "model1_vs_model5_summary.csv"),
  row.names = FALSE
)


# 4.5 Save genes significant in both models.

write.csv(
  data.frame(gene_id = sig_both_1_5),
  file = file.path(output_folder, "model1_model5_shared_genes.csv"),
  row.names = FALSE
)

### MODEL 5 COMPLETE



### MODEL 5 VS MODEL 1 — INTERPRETATION
### Model 1 identified 1,871 significant genes, while Model 5 identified 1,660.
### A total of 1,522 genes were significant in both models, representing
### 81.35% of the significant genes identified in Model 1.
### All 1,522 shared significant genes showed associations in the same direction.
### These findings indicate substantial agreement between the two approaches.
### However, the analyses use the same participants and are not independent validation.





