###############################################################################
# PCA of day-0 (baseline) RNA-seq
###############################################################################
# Goal: To  identify the main sources of variation in baseline gene expression
#  and to decide which covariates to include in the association models.


source("environment.R")   # adjust the path if environment.R is in another folder

### Check that everything needed from environment.R exists
needed <- c("logcpm", "sample_info", "flow_mon_clean", "dem_clean",
            "output_folder", "study_labels")
missing <- needed[!sapply(needed, exists)]
if (length(missing) > 0) stop("Not created by environment.R: ", 
                              paste(missing, collapse = ", "))

# Columns of logcpm must be in the same participant order as sample_info
stopifnot(all(colnames(logcpm) == sample_info$participant_id))


###############################################################################
# 1. Participant Info
###############################################################################
# selecting columns needed for PCA plots
participant_info <- mon_dem_clean %>%
  select(participant_id, study_accession, biological_sex)

#Check if there is only one row per participant
stopifnot(!any(duplicated(participant_info$participant_id)))

# Every participant in the expression matrix should have this information
#should produce TRUE
all(colnames(logcpm) %in% participant_info$participant_id)

###############################################################################
# 2. PCA on the 1,000 most variable genes
###############################################################################
# Low-variance genes add noise rather than structure (standard for RNA-seq PCA).
# Centred but not scaled: log2-CPM is already on a common scale, and scaling
# would increase the weight of noisy low-variance genes.

#creating a vector of the 1000 most variable genes
top_genes <- names(sort(apply(logcpm, 1, var), decreasing = TRUE))[1:1000] 
# creating the pca
pca <- prcomp(t(logcpm[top_genes, ]), center = TRUE, scale. = FALSE)       


# % variance explained by PC1-PC5 (a quick illustrative look)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1) 
var_expl[1:5] 

#generating scores to position each participant for the plot
pca_scores <- as_tibble(pca$x[, 1:10], rownames = "participant_id") %>% 
  left_join(sample_info, by = "participant_id") 


###############################################################################
# 3. PCA coloured by study
###############################################################################
# Determine if there was a batch effect as studies were sequenced
#separately.

p_pca_study <- ggplot(pca_scores, aes(PC1, PC2, colour = study_accession)) +
  geom_point(size = 3) +
  scale_colour_discrete(labels = study_labels) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Study")

p_pca_study

#Plot result: The two studies are completely seperated along PC1. 
#Almost two thirds-of all variation comes form the study sample.
#Study must be inlcuded as a covariate

#PC2 also splits each study into two groups. The two way split appears as a 
#binary characteristic, most likely sex, but need to verify

ggsave(file.path(output_folder, "pca_by_study.png"),
       p_pca_study, width = 8, height = 6, dpi = 300)



###############################################################################
# 4. Verfication that PC2 corresponds to Sex
###############################################################################
# PC2 splits into two groups at about 0 (see section 3), so 0 is used as the cut-off.
sex_check <- pca_scores %>%
  mutate(pc2_group = if_else(PC2 > 0, "PC2 high", "PC2 low"))


# Cross-table: if PC2 is sex, each PC2 group should be almost entirely one sex.
# useNA shows participants with missing sex.
table(sex_check$pc2_group, sex_check$biological_sex, useNA = "ifany")

# Result: each PC2 group is entirely one sex.

# Same table split by study: checks that both sexes are present in both studies
table(sex_check$pc2_group, sex_check$biological_sex, sex_check$study_accession)

# Result: both sexes are present in both studies.

# PCA Plot coloured by Study and shape by sex
p_pca_sex <- ggplot(sex_check, aes(PC1, PC2, colour = biological_sex, 
                                   shape = study_accession)) +
  geom_point(size = 3) +
  scale_shape_discrete(labels = study_labels) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Study", shape = "Sex") +
  theme(legend.key.size = unit(0.8, "cm"))   

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
