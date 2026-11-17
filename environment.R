############Environment Setup############

#####Set working directory##### 
# Dynamic user name
username <- Sys.getenv("HOME")

# Define path
input_path <- file.path(username, "Documents", "GitHub", "omics_m2phds", "data")

# Set wd
setwd(input_path)

#####Load packages##### 
pacman::p_load(tidyverse)

#####Load datasets#####  
# Particpants
dem <- read.csv("participants.csv")
# Vaccination history 
dem_vacc <- read.csv("participant_vaccineHistory.csv")
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

#####Filter and Merge Datasets##### 
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

### Create a clean demographic dataset with the study sample
dem_clean <- dem %>%
  semi_join(flow_mon_clean, by = keys)

### Create a clean vaccine history dataset with the study sample
vacc_clean <- dem_vacc %>%
  semi_join(flow_mon_clean, by = "subject")

#####Data datasets#####
# Standardized value column. 2024_UGA study was multiplied by 100 to convert from proportion to percentage.  
flow_mon_clean <- flow_mon_clean %>% 
  mutate(value_clean = if_else(study_accession == "2024_UGA", value*100, value))
