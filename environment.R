###############################################################################
# ENVIRONMENT SETUP 
###############################################################################

#####Set working directory##### 
# Dynamic user name
username <- Sys.getenv("HOME")

# Define path
input_path <- file.path(username, "Documents", "GitHub", "omics_m2phds", "data")

# Set wd
setwd(input_path)

#####Load packages##### 
pacman::p_load(tidyverse, here, edgeR)

#####Load datasets#####  
# Particpants
dem <- read.csv("participants.csv")
# Vaccination history 
vacc <- read.csv("participant_vaccineHistory.csv")
# Flow Cytometry 
flow <- read.csv("publicData_ex_vivo_flow.csv")
# RNA SEQ Datasets
rna_uga_2019 <- read.csv("publicData_rnaseq_2019_UGA.csv")
rna_uga_2020 <- read.csv("publicData_rnaseq_2020_UGA.csv")
rna_uga_2024 <- read.csv("publicData_rnaseq_2024_UGA.csv")
rna_SDY224_Bcells <- read_csv("publicData_rnaseq_SDY224_Bcells.csv")
rna_SDY2867 <- read_csv("publicData_rnaseq_SDY2867.csv")
rna_SDY2941 <- read_csv("publicData_rnaseq_SDY2941.csv")
rna_260821 <- read_csv("publicData_serology_260821.csv")

# Create output folder
dir.create(here("output"), showWarnings = FALSE)
output_folder <- here("output")


###############################################################################
# Set global theme 
###############################################################################

theme_set(theme_bw() +
            theme(title = element_text(size = 15, color = "grey30"),
                  axis.title.x = element_text(size = 15, color = "grey30"),
                  axis.title.y = element_text(size = 15, color = "grey30"),
                  axis.text.x = element_text(size = 15, color = "grey30"),
                  axis.text.y = element_text(size = 15, color = "grey30"),
                  panel.grid = element_blank(),
                  strip.text = element_text(size = 15, color = "grey30"),
                  strip.background = element_blank(), 
                  legend.position = "bottom",
                  legend.key.size = unit(3, "cm"),
                  legend.text = element_text(size = 15, color = "grey30"))) 

###############################################################################
# Filter and Merge Datasets
###############################################################################
### Take only monocytes at day 7
flow_mon <- flow %>%
  filter(name == "Classical_monocytes", timepoint == "7")

flow_mon_all_timepoints <- flow %>%
  filter(name == "Classical_monocytes")

### Find participants with monocyte data that have RNA seq data
# Define merge keys
keys <- c("participant_id", "study_accession", "subject")

# Make a list of unique IDs 
mon_ids <- flow_mon %>% distinct(across(all_of(keys)))

# Create a list of all RNA seq datasets

rna_list <- list(
  uga_2019     = rna_uga_2019,
  uga_2020     = rna_uga_2020,
  uga_2024     = rna_uga_2024,
  SDY224_Bcell = rna_SDY224_Bcells,
  SDY2867      = rna_SDY2867,
  SDY2941      = rna_SDY2941
)

# Filter RNA seq for day 0 (baseline) and join with monocyte data by keys
rna_mon <- rna_list %>%
  map(~ filter(.x, as.character(timepoint) == "0")) %>%
  map(~ semi_join(.x, mon_ids, by = keys)) %>%
  map(~ mutate(.x, comments = as.character(comments))) %>%
  bind_rows(.id = "rna_dataset")

flow_mon_rna_clean <- flow_mon %>%
  inner_join(rna_mon, by = keys,
             suffix = c("_flow", "_rna"),
             relationship = "many-to-many")

# Check the number of participants with both monocyte and RNA seq data
n_distinct(flow_mon_rna_clean$participant_id)

### Create a clean monocyte dataset with the study sample
flow_mon_clean <- flow_mon %>%
  semi_join(rna_mon, by = keys)

### Create a clean vaccine history dataset with a binary sum column (yes if particpant has received at least 1 vaccine)
vacc_clean <- vacc %>%
  mutate(vacc_binary = if_else(startsWith(received, "Yes"), 1, 0)) %>%
  group_by(subject) %>%
  summarise(
    vacc_any = as.integer(any(vacc_binary == 1, na.rm = TRUE)),
    .groups = "drop"
  )
  
### Merge demographic, vaccine, and monocyte datasets 
mon_dem_clean <- flow_mon_clean %>%
  left_join(dem, by = keys) %>%
  left_join(vacc_clean, by = "subject", relationship = "many-to-many")

###############################################################################
# Data Cleaning 
###############################################################################

### Standardized value column. 2024_UGA study was multiplied by 100 to convert from proportion to percentage.  
# Check that 2024 UGA values are a proportion 
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

# Clean data
mon_dem_clean <- mon_dem_clean %>% 
  mutate(value_clean = if_else(study_accession == "2024_UGA", value*100, value),
# Create binary column for vaccine history 
         vacc_binary = if_else(startsWith(received, "Yes"), 1, 0))

# Check harmonized results  
flow_mon_clean %>%
  group_by(study_accession) %>%
  summarise(min = min(value_clean), median = median(value_clean), max = max(value_clean))

### Filter for Genes commonly expressed 

# 1. Participant order: one row per participant in the study sample
sample_info <- flow_mon_clean %>%
  distinct(participant_id, study_accession) %>%
  arrange(participant_id)

nrow(sample_info)   # expect 77

# 2. Genes measured in BOTH studies
# A gene missing from one study would otherwise appear as 0 expression for all
# of that study's participants, creating a fake "study" signal.
common_genes <- rna_mon %>%
  filter(participant_id %in% sample_info$participant_id) %>%
  distinct(study_accession, ensembl_gene_id) %>%
  count(ensembl_gene_id) %>%
  filter(n == n_distinct(sample_info$study_accession)) %>%
  pull(ensembl_gene_id)

length(common_genes)   # report in methods

# 3. Raw count matrix: genes (rows) x participants (columns)
counts_mat <- rna_mon %>%
  filter(participant_id %in% sample_info$participant_id,
         ensembl_gene_id %in% common_genes) %>%
  select(participant_id, ensembl_gene_id, counts) %>%
  pivot_wider(names_from  = participant_id,
              values_from = counts,
              values_fn   = sum) %>%                
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

counts_mat <- counts_mat[, sample_info$participant_id]   # fixed participant order

dim(counts_mat)                 # genes x 77
sum(is.na(counts_mat))          # should be 0

# 4. Filter low-expressed genes
# filterByExpr keeps genes with enough reads in enough samples,
# taking group sizes (here: study) into account.
dge  <- DGEList(counts = counts_mat, group = sample_info$study_accession)
keep <- filterByExpr(dge, group = sample_info$study_accession)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

sum(keep)   # number of genes after filtering (report in methods)

###############################################################################
# TMM normalisation of day-0 RNA-SEQ (study sample only)
###############################################################################

# 1. TMM normalisation
# Computes a scaling factor per sample to correct for sequencing depth and
# composition differences between samples.
dge <- calcNormFactors(dge, method = "TMM")

summary(dge$samples$norm.factors)   # usually close to 1 (e.g. 0.7-1.3)

# 2. Log2 counts per million (TMM-adjusted)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

dim(logcpm)                                           # genes x 77
all(colnames(logcpm) == sample_info$participant_id)  # must be TRUE

### 7. Save for use in other scripts
# Saved in data/ (ignored by Git, so it is not pushed: confidential data)
saveRDS(dge,         here("data", "dge_tmm.rds"))     # for limma-voom (needs counts + norm factors)
saveRDS(logcpm,      here("data", "logcpm_tmm.rds"))  # for PCA and prediction
saveRDS(sample_info, here("data", "sample_info.rds"))

