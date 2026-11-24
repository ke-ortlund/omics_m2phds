###############################################################################
# PCA of day-0 (baseline) RNA-seq
###############################################################################
# Aim: identify the main sources of variation in baseline gene expression
# (study, sex) and check whether they relate to day-7 classical monocytes,
# to decide which covariates to include in the association models.
#
# Uses objects created in environment.R (filtering + TMM normalisation are done
# there, not repeated here):
#   logcpm         : TMM-normalised log2-CPM matrix, genes x participants
#   sample_info    : participant_id + study_accession, same order as logcpm columns
#   flow_mon_clean : day-7 classical monocytes, value_clean harmonised to %
#   dem_clean      : demographics (biological sex)
#   output_folder, study_labels

source("environment.R")   # adjust the path if environment.R is in another folder

### Check that everything needed from environment.R exists
needed <- c("logcpm", "sample_info", "flow_mon_clean", "dem_clean",
            "output_folder", "study_labels")
missing <- needed[!sapply(needed, exists)]
if (length(missing) > 0) stop("Not created by environment.R: ", paste(missing, collapse = ", "))

# Columns of logcpm must be in the same participant order as sample_info
stopifnot(all(colnames(logcpm) == sample_info$participant_id))


###############################################################################
# 1. Add the day-7 monocyte value to sample_info (one value per participant)
###############################################################################
# value_clean is already harmonised to % in environment.R (2024_UGA x 100).
# Mean in case a participant has more than one row.
sample_info <- sample_info %>%
  left_join(flow_mon_clean %>%
              group_by(participant_id) %>%
              summarise(monocyte = mean(value_clean, na.rm = TRUE)),
            by = "participant_id")

sum(is.na(sample_info$monocyte))   # should be 0


###############################################################################
# 2. PCA on the 1,000 most variable genes
###############################################################################
# Low-variance genes add noise rather than structure (standard for RNA-seq PCA).
# Centred but not scaled: log2-CPM is already on a common scale, and scaling
# would up-weight noisy low-variance genes.

top_genes <- names(sort(apply(logcpm, 1, var), decreasing = TRUE))[1:1000]   # 1,000 most variable genes
pca <- prcomp(t(logcpm[top_genes, ]), center = TRUE, scale. = FALSE)       # participants in rows

# % of variance explained by each PC
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1) # % variance explained by PC1-PC5 (a quick illustrative look)
var_expl[1:5]   # PC1-PC5

# PC scores: the position of each participant on the PCs, with participant info
pca_scores <- as_tibble(pca$x[, 1:10], rownames = "participant_id") %>% 
  left_join(sample_info, by = "participant_id")  #generating scores to posting each participant for the plot


###############################################################################
# 3. PCA coloured by study
###############################################################################
# Checks for a batch effect (the two studies were sequenced separately).

p_pca_study <- ggplot(pca_scores, aes(PC1, PC2, colour = study_accession)) +
  geom_point(size = 3) +
  scale_colour_discrete(labels = study_labels) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Study")

p_pca_study

# Result: the two studies are completely separated along PC1 (62.3% of the
# variance), so study must be included as a covariate.
# PC2 (6.7%) splits each study into two groups. This looks like a binary
# characteristic, most likely sex; checked in section 5.

ggsave(file.path(output_folder, "pca_by_study.png"),
       p_pca_study, width = 8, height = 6, dpi = 300)


###############################################################################
# 4. PCA coloured by day-7 monocytes
###############################################################################
# Shape = study, so we can see whether there is a monocyte gradient within
# each study and not only between studies.

p_pca_day7 <- ggplot(pca_scores, aes(PC1, PC2, colour = monocyte, shape = study_accession)) +
  geom_point(size = 3) +
  scale_colour_viridis_c() +
  scale_shape_discrete(labels = study_labels) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Day-7 monocytes (%)", shape = "Study") +
  theme(legend.key.size = unit(0.8, "cm"))   # smaller legend keys so the colour bar fits

p_pca_day7

# Result: monocyte levels are generally lower in UGA 2024 than in SDY 2867, so
# study is linked to both gene expression and the outcome (confirms the need to
# adjust for study). The two PC2 groups have similar monocyte levels within each
# study, and there is no clear monocyte gradient along PC1 or PC2 within studies.

ggsave(file.path(output_folder, "pca_by_day7.png"),
       p_pca_day7, width = 8, height = 6, dpi = 300)


###############################################################################
# 5. Does PC2 correspond to biological sex?
###############################################################################
# Add sex from dem_clean (distinct() so each participant appears once).
# PC2 splits into two groups at about 0 (see section 3), so 0 is used as the cut-off.

sex_check <- pca_scores %>%
  left_join(dem_clean %>% distinct(participant_id, biological_sex), by = "participant_id") %>%
  mutate(pc2_group = if_else(PC2 > 0, "PC2 high", "PC2 low"))

# Cross-table: if PC2 is sex, each PC2 group should be almost entirely one sex.
# useNA shows participants with missing sex.
table(sex_check$pc2_group, sex_check$biological_sex, useNA = "ifany")
# Result: each PC2 group is (almost) entirely one sex.

# Same table split by study: checks that both sexes are present in both studies
table(sex_check$pc2_group, sex_check$biological_sex, sex_check$study_accession)
# Result: both sexes are present in both studies.

# Plot: PCA coloured by sex
p_pca_sex <- ggplot(sex_check, aes(PC1, PC2, colour = biological_sex, shape = study_accession)) +
  geom_point(size = 3) +
  scale_shape_discrete(labels = study_labels) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Sex", shape = "Study") +
  theme(legend.key.size = unit(0.8, "cm"))   # smaller legend keys so both legends fit

p_pca_sex

# Result: PC2 corresponds to sex, so sex is also included as a covariate.

ggsave(file.path(output_folder, "pca_by_sex.png"),
       p_pca_sex, width = 8, height = 6, dpi = 300)

# Participants whose sex does not match the majority of their PC2 group
# (possible sample mix-ups or recording errors)
sex_check %>%
  add_count(pc2_group, biological_sex) %>%
  group_by(pc2_group) %>%
  filter(n < max(n)) %>%
  ungroup() %>%
  select(participant_id, study_accession, PC2, pc2_group, biological_sex)

# Result: one female falls below the PC2 = 0 cut-off (PC2 ≈ -10) but is clearly
# separate from the male cluster; not a sample mix-up, no action needed.
