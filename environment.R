###### Environment Setup 
### Set working directory
# Dynamic user name
username <- Sys.getenv("HOME")

# Define path
input_path <- file.path(username, "Documents", "GitHub", "omics_m2phds", "data")

# Set wd
setwd(input_path)

###### Load packages
pacman::p_load(tidyverse)

###### Load datasets 
# Particpants
ppt <- read.csv("participants.csv")
# Vaccination history 
ppt_vacc <- read.csv("participant_vaccineHistory.csv")
# Flow Cytometry 
flow <- read.csv("publicData_ex_vivo_flow.csv")
### RNA SEQ Datasets
rna_uga_2019 <- read.csv("publicData_rnaseq_2019_UGA.csv")
rna_uga_2020 <- read.csv("publicData_rnaseq_2020_UGA.csv")
rna_uga_2024 <- read.csv("publicData_rnaseq_2024_UGA.csv")
rna_SDY224_Bcells <- read_csv("publicData_rnaseq_SDY224_Bcells.csv")
rna_SDY2867 <- read_csv("publicData_rnaseq_SDY2867.csv")
rna_SDY2941 <- read_csv("publicData_rnaseq_SDY2941.csv")
rna_260821 <- read_csv("publicData_serology_260821.csv")

###### Filter and Merge Datasets

### Filter flow cytometry for classic monocytes 
### 3 studies have monocytes: rna_SDY224_Bcells, rna_SDY224, rna_SDY2941
flow_mon <- flow %>%
  filter(name == "Intermediate_monocytes")
