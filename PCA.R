###############################################################################
# PCA of day-0 RNA-seq
###############################################################################


###### DELETE THIS SECTION, use flow_mon_clean, it was already completed in the environment script

### 1. Add day-7 monocyte value to sample_info (one value per participant)
# 2024_UGA x100 (proportion -> %), as above. Mean in case a participant has >1 row.
#sample_info <- sample_info %>%
 # left_join(flow_mon_clean %>%
#              mutate(value_clean = if_else(study_accession == "2024_UGA", value * 100, value)) %>%
#              group_by(participant_id) %>%
#              summarise(monocyte = mean(value_clean, na.rm = TRUE)),
#            by = "participant_id")


### 2. Run PCA on the 1000 most variable genes
# Low-variance genes add noise, not structure (standard for RNA-seq PCA).
# Centred but not scaled: logCPM is already on a common scale; scaling would
# up-weight noisy low-variance genes.
top_genes <- names(sort(apply(logcpm, 1, var), decreasing = TRUE))[1:1000] #creating a vector of the 1000 most variable genes
pca <- prcomp(t(logcpm[top_genes, ]), center = TRUE, scale. = FALSE) #creating the pca


var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
var_expl[1:5]   # % variance explained by PC1-PC5 (a quick illustrative look)

pca_scores <- as_tibble(pca$x[, 1:10], rownames = "participant_id") %>%
  left_join(sample_info, by = "participant_id") #generating scores to posting each participant for the plot



### 3. Plots
# By study: checks for a batch effect (studies sequenced separately)
p_pca_study <- ggplot(pca_scores, aes(PC1, PC2, colour = study_accession)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Study")

p_pca_study


#Plot result: The two studies are completely seperated along PC1. Almost two thirds-of all variation comes form the study sample. Study must be inlcuded as a covariate

#PC2 also splits each study into two groups. The two way split appears as a bianry charctertic, most likely sex, but needto verify

ggsave(file.path(output_folder, "pca_by_study.png"), #saving output
       p_pca_study, width = 8, height = 6, dpi = 300)


##### Doesn't run!
# By day-7 monocytes: shape = study to see the gradient within each study
p_pca_day7 <- ggplot(pca_scores, aes(PC1, PC2, colour = monocyte, shape = study_accession)) +
  geom_point(size = 3) +
  scale_colour_viridis_c() +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = "Day-7 monocytes (%)", shape = "Study") +
  theme(legend.key.size = unit(0.8, "cm")) # 3 cm from global theme is too big for a colour bar

p_pca_day7

#Plot result: Monocyte levels also differ by study. However PC2 doesn't seem related to monocyte levels. the upper and lower levels seem to have simialr colours in both studies. Confirms need to adjust for study

ggsave(file.path(output_folder, "pca_by_day7.png"), #save output
       p_pca_day7, width = 8, height = 6, dpi = 300)

#### Doesn'r run!
### 4. Does PC2 match biological sex?
# Add sex from dem (distinct() so each participant appears once).
# PC2 splits into two groups at about 0 (see PCA plot), so we use 0 as the cut-off.
sex_check <- pca_scores %>%
  left_join(dem_clean %>% distinct(participant_id, biological_sex), by = "participant_id") %>%
  mutate(pc2_group = if_else(PC2 > 0, "PC2 high", "PC2 low"))

# Cross-table: if PC2 is sex, each PC2 group should be almost entirely one sex.
# useNA shows participants with missing sex.
table(sex_check$pc2_group, sex_check$biological_sex, useNA = "ifany") #yes basically every PC2 is entirely one sex

# Same table split by study: also checks that both sexes are present in both studies
table(sex_check$pc2_group, sex_check$biological_sex, sex_check$study_accession) #both sexes are present in both studies

# Plot: PCA coloured by sex
p_pca_sex <- ggplot(sex_check, aes(PC1, PC2, colour = biological_sex, shape = study_accession)) +
  geom_point(size = 3) +
  labs(colour = "Sex", shape = "Study") + theme(legend.key.size = unit(0.8, "cm")) #result showed pc2 is sex, so need to adjust for sex as well

p_pca_sex

ggsave(file.path(output_folder, "pca_by_sex.png"), #save output
       p_pca_sex, width = 8, height = 6, dpi = 300)

# Participants whose sex doesn't match the majority of their PC2 group
# (possible sample mix-ups or recording errors worth checking)
sex_check %>%
  add_count(pc2_group, biological_sex) %>%
  group_by(pc2_group) %>%
  filter(n < max(n)) %>%
  ungroup() %>%
  select(participant_id, study_accession, PC2, pc2_group, biological_sex)

#one female that doesn't match PC2 group. don't know if we need to account for that?
