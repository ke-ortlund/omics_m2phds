###############################################################################
# GENE EXPRESSION ANALYSIS
###############################################################################
names(flow_mon_rna_clean)


# Check the units used for monocyte measurements.
# We need to confirm that the values represent percentages.
table(flow_mon_rna_clean$unit, useNA = "ifany")



###############################################################################
# 2. VERIFY MONOCYTE MEASUREMENTS
###############################################################################
# Check monocyte values by study.
# This helps identify differences in measurement scales
# before performing gene expression association analysis.
flow_mon_rna_clean %>%
  group_by(study_accession) %>%
  summarise(
    n_participants = n_distinct(participant_id),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE)
  )


### gene_analysis_data is now flow_mon_clean
### monocyte_pct is no value_clean 

###############################################################################
# 4. VERIFICATION OF MONOCYTE PERCENTAGE HARMONIZATION
###############################################################################
#### Compare monocyte percentages across the two studies.
# Use distinct() to avoid counting each participant repeatedly for every gene.

gene_analysis_data %>%
  distinct(study_accession, participant_id, monocyte_pct) %>%
  group_by(study_accession) %>%
  summarise(
    n_participants = n_distinct(participant_id),
    min_pct = min(monocyte_pct, na.rm = TRUE),
    max_pct = max(monocyte_pct, na.rm = TRUE),
    mean_pct = mean(monocyte_pct, na.rm = TRUE)
  )



###############################################################################
# 5. VERIFICATION OF RNA-SEQ DATA STRUCTURE
###############################################################################
#### Check the number of genes and participants in each study.
gene_analysis_data %>%
  group_by(study_accession) %>%
  summarise(
    n_participants = n_distinct(participant_id),
    n_genes = n_distinct(ensembl_gene_id),
    n_rows = n()
  ) %>%
  print()


###############################################################################
# 6. CHECK FOR DUPLICATED GENE MEASUREMENTS
###############################################################################

# Check whether the same participant has multiple measurements
# for the same gene within each study.
#
# Each participant-gene combination should appear only once
# if there is one RNA-seq sample per participant.

gene_analysis_data %>%
  count(study_accession, participant_id, ensembl_gene_id) %>%
  filter(n > 1) %>%
  summarise(n_duplicates = n()) %>%
  print()


###############################################################################
# 7. IDENTIFICATION OF GENES SHARED BETWEEN STUDIES
###############################################################################

# Extract the unique Ensembl gene IDs from each study.
genes_uga <- gene_analysis_data %>%
  filter(study_accession == "2024_UGA") %>%
  pull(ensembl_gene_id) %>%
  unique()

genes_sdy <- gene_analysis_data %>%
  filter(study_accession == "SDY2867") %>%
  pull(ensembl_gene_id) %>%
  unique()

# Identify genes present in both studies.
shared_genes <- intersect(genes_uga, genes_sdy)

# Display the number of shared genes.
print(length(shared_genes))
###Every gene measured in SDY2867 is also present in 2024_UGA. still there qre 24030 
# genes excluded from the combined analysis. we might want to analyze them separately (for 24 participants)


###############################################################################
# 8. VERIFICATION OF RAW RNA-SEQ COUNTS
###############################################################################

# Check whether RNA-seq counts contain missing values, negative values,
# or non-integer values.
#
# This is important because the voom-limma pipeline requires
# non-negative RNA-seq count data as input.

gene_analysis_data %>%
  group_by(study_accession) %>%
  summarise(
    n_missing = sum(is.na(counts)),
    n_negative = sum(counts < 0, na.rm = TRUE),
    n_non_integer = sum(counts != round(counts), na.rm = TRUE),
    min_counts = min(counts, na.rm = TRUE),
    max_counts = max(counts, na.rm = TRUE)
  ) %>%
  print()


###############################################################################
# 9. VERIFICATION OF RNA-SEQ LIBRARY SIZES
###############################################################################

# Calculate the total number of RNA-seq counts for each participant.
# Only genes shared between the two studies are included to ensure
# that library sizes are calculated from the same set of genes.

library_sizes <- gene_analysis_data %>%
  filter(ensembl_gene_id %in% shared_genes) %>%
  group_by(study_accession, participant_id) %>%
  summarise(
    total_counts = sum(counts),
    .groups = "drop"
  )

# Summarise library sizes separately for each study.
library_sizes %>%
  group_by(study_accession) %>%
  summarise(
    n_participants = n(),
    min_library = min(total_counts),
    median_library = median(total_counts),
    max_library = max(total_counts)
  ) %>%
  print()

#### Interpretation of library sizes and normalization

# Library sizes vary from 7.6 to 35.2 million counts (4.6-fold difference).
# Median library sizes also differ between studies:
# 2024_UGA: 17.44M | SDY2867: 22.49M.
#
# Higher counts may reflect sequencing depth rather than biological differences.
# Therefore, normalization is necessary before comparing gene expression.
#
# Why TMM instead of TPM?
# TPM adjusts for gene length and sequencing depth but not composition bias.
# TMM adjusts for RNA composition differences between samples.
# For voom-limma, we use raw counts -> TMM -> voom -> limma.
#
# Note: We still need to confirm that 'counts' are truly raw counts.


###############################################################################
# 10.1 PREPARATION OF THE RNA-SEQ COUNT MATRIX
###############################################################################

# Keep only genes measured in both studies.
# Create a unique sample ID by combining study and participant identifiers.
# This prevents confusion if participant IDs are repeated across studies.
# The original gene_analysis_data remains unchanged.

gene_counts_long <- gene_analysis_data %>%
  filter(ensembl_gene_id %in% shared_genes) %>%
  mutate(
    sample_id = paste(study_accession, participant_id, sep = "_")
  )

# Verify the number of samples and genes available for analysis.

gene_counts_long %>%
  summarise(
    n_samples = n_distinct(sample_id),
    n_genes = n_distinct(ensembl_gene_id)
  ) %>%
  print()


###############################################################################
# 10.2. CREATION OF THE RNA-SEQ COUNT MATRIX
###############################################################################

# Convert the dataset from long to wide format.
# Rows represent genes, columns represent participants,
# and values represent RNA-seq counts.

count_matrix <- gene_counts_long %>%
  select(ensembl_gene_id, sample_id, counts) %>%
  pivot_wider(
    names_from = sample_id,
    values_from = counts
  ) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

# Verify the dimensions of the matrix.
dim(count_matrix)

# Check whether any counts are missing.
sum(is.na(count_matrix))


###############################################################################
# 11. CREATION OF SAMPLE METADATA
###############################################################################

# Create one row per participant, containing the study identifier
# and Day 7 classical monocyte percentage.

sample_metadata <- gene_counts_long %>%
  distinct(sample_id, study_accession, monocyte_pct)

# Check that each sample has exactly one metadata row.
nrow(sample_metadata)

# Check for missing monocyte percentages.
sum(is.na(sample_metadata$monocyte_pct))

# Check whether any sample appears more than once.
sum(duplicated(sample_metadata$sample_id))



###############################################################################
# 12. ALIGN SAMPLE METADATA WITH THE COUNT MATRIX
###############################################################################

# Reorder metadata to match the columns of the RNA-seq count matrix.
sample_metadata <- sample_metadata[
  match(colnames(count_matrix), sample_metadata$sample_id),
]

# Verify that every sample is correctly aligned.
all(sample_metadata$sample_id == colnames(count_matrix))
###alignment because limma must know which monocyte percentage belongs
#to which participant's gene expression data.



###############################################################################
# 13. FILTERING LOWLY EXPRESSED GENES
###############################################################################

# Create the study-adjusted design matrix.
# Monocyte percentage is our variable of interest.
# Study is included to account for differences between cohorts.

design <- model.matrix(
  ~ monocyte_pct + study_accession,
  data = sample_metadata
)

# Identify genes with sufficient expression for analysis.
keep <- edgeR::filterByExpr(
  count_matrix,
  design = design
)

# Check how many genes are retained and excluded.
print(table(keep))
#### 18441 are the gene wer retain, 36461 are the genes with really low expression or no detecbtale expression


###############################################################################
# 14. CREATION OF THE FILTERED COUNT MATRIX
###############################################################################

# Retain only genes with sufficient expression.
# The original count_matrix remains unchanged.

filtered_counts <- count_matrix[keep, ]

# Verify the dimensions of the filtered matrix.
dim(filtered_counts)
## 18,441 genes × 77 participants ready for normalization



###############################################################################
##VERIFICATION OF ORIGINAL RNA-SEQ DATA
###############################################################################
### data are already normalized but for voomlimma might be necessary another type or normal
# Inspect the column names of the two original RNA-seq datasets.
# We want to identify whether raw counts and normalized values
# were provided separately.

names(read.csv("publicData_rnaseq_2024_UGA.csv", nrows = 5))

names(read.csv("publicData_rnaseq_SDY2867.csv", nrows = 5))

###############################################################################
#CHECK RNA-SEQ PREPROCESSING INFORMATION
###############################################################################

# Look for information about normalization or preprocessing
# in the comments provided with the original datasets.

gene_analysis_data %>%
  group_by(study_accession) %>%
  summarise(
    comments = paste(unique(na.omit(comments_rna)), collapse = "; "),
    .groups = "drop"
  ) %>%
  print(width = Inf)

readLines("README.md")

##############################################################################
#RNA-seq normalization assumption

# The datasets contain both counts and TPM values.
# TPM is already normalized, while counts appear to represent raw counts.
# Since preprocessing documentation is unavailable, we assume counts are raw.
# We therefore use TMM normalization followed by voom-limma.
###############################################################################


###############################################################################
# 15. TMM NORMALIZATION
###############################################################################

# Create an edgeR DGEList object containing the filtered raw counts.
# This object stores the count matrix and library-size information.

dge <- edgeR::DGEList(counts = filtered_counts)

# Calculate TMM normalization factors to adjust for differences
# in RNA composition between samples.

dge <- edgeR::calcNormFactors(dge, method = "TMM")

# Display the normalization factors for the first six samples.
head(dge$samples)


###############################################################################
# 16. VOOM TRANSFORMATION
###############################################################################

# Transform RNA-seq counts into log2-CPM using TMM-adjusted library sizes.
# Estimate the mean-variance relationship and calculate precision weights
# for the subsequent gene-wise linear models.

v <- limma::voom(
  dge,
  design = design,
  plot = TRUE
)
# Lowly expressed genes show greater variability, which decreases
# as expression increases and stabilizes at higher expression levels.
# voom models this relationship and assigns precision weights
# for the subsequent gene-wise linear models.


###############################################################################
# 17. GENE-WISE LINEAR MODELS USING LIMMA
###############################################################################

# Fit a weighted linear model for each gene.
# The design matrix includes monocyte percentage and study accession.
# voom precision weights are automatically used in the model.

fit <- limma::lmFit(v, design)

# Apply empirical Bayes moderation to improve variance estimation
# across genes.

fit <- limma::eBayes(fit)

# Check the coefficients included in our models.
colnames(fit$coefficients)


###############################################################################
# 18. EXTRACTION OF GENE-WISE ASSOCIATION RESULTS
###############################################################################

# Extract results for the association between baseline gene expression
# and Day 7 monocyte percentage, adjusting for study.
#
# BH correction controls the false discovery rate across genes.
# Results are ordered by statistical significance.

gene_results <- limma::topTable(
  fit,
  coef = "monocyte_pct",
  number = Inf,
  adjust.method = "BH"
)

# Display the first 10 genes.
head(gene_results, 10)

#### Interpretation of gene-wise results

# Several genes show significant associations after BH correction.
# Positive logFC: higher monocyte % is associated with higher gene expression.
# Negative logFC: higher monocyte % is associated with lower gene expression.
# Results indicate association, not causality.

###############################################################################
# 19. IDENTIFICATION OF SIGNIFICANT GENES
###############################################################################

# Identify genes significantly associated with Day 7 monocyte percentages
# using a Benjamini-Hochberg adjusted p-value threshold of 0.05.

significant_genes <- gene_results %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  filter(adj.P.Val < 0.05)

# Count the number of significant genes.
nrow(significant_genes)

# Count positive and negative associations.
significant_genes %>%
  summarise(
    total_significant = n(),
    positive_associations = sum(logFC > 0),
    negative_associations = sum(logFC < 0)
  ) %>%
  print()

#### Interpretation of significant genes

# 3,595 of 18,441 genes (19.5%) are significant at FDR < 0.05.
# 1,768 show positive and 1,827 show negative associations.
# Results suggest widespread associations between baseline gene
# expression and Day 7 monocyte percentages, adjusting for study.


###############################################################################
# 20. VOLCANO PLOT OF GENE-WISE ASSOCIATIONS
###############################################################################

# Prepare results and classify genes by significance and direction.
volcano_data <- gene_results %>%
  mutate(
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Positive",
      adj.P.Val < 0.05 & logFC < 0 ~ "Negative",
      TRUE ~ "Not significant"
    )
  )

# Visualize associations between baseline gene expression
# and Day 7 classical monocyte percentages.

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c(
      "Positive" = "red",
      "Negative" = "blue",
      "Not significant" = "grey"
    )
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Baseline Gene Expression and Day 7 Monocyte Percentage",
    subtitle = "Gene-wise associations adjusted for study",
    x = "Association coefficient (log2 expression per 1% monocyte increase)",
    y = "-log10(adjusted p-value)",
    color = "Association"
  ) +
  theme_minimal()

#### Volcano plot interpretation

# 3,595 genes are significant at FDR < 0.05.
# Positive and negative associations are approximately balanced.
# Most effect estimates are small, with a few larger associations.

###############################################################################
################################################################################
#auestions
#nromalization issue
#1 i used monocyte as X and gene expressione as Y, should we do the opposite?
#2 should we use other packages rather than voom-limma (with other normalization?)
#3 sensitivity analysis: reverse x-y, do the gene expr analysis on the two different studies separetly and co,pare results?
#4 what about the genes that where not matching in the two studies? should i do that and maybe use it as sensitivuty analysis?
#5 check why i needed to do ar,onization for monocyte percentage
#############################################################################



