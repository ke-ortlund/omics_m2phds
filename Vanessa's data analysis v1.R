
###Cleaning

# Standardizing value column. 2024_UGA study was multiplied by 100 to convert from proportion to percentage.
flow_mon_rna_clean <- flow_mon_rna_clean %>% 
  mutate(value_clean = if_else(study_accession == "2024_UGA", value*100, value))

###Packages
#load packages
library(DESeq2)     
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)


###Aplying Variance Stabilization

#making a counts table
count_table <- flow_mon_rna_clean %>%
  distinct(participant_id, ensembl_gene_id, counts) %>%
  pivot_wider(names_from = participant_id, values_from = counts) %>%
  tibble::column_to_rownames("ensembl_gene_id") %>%
  as.matrix()
count_table <- count_table[complete.cases(count_table), ]

#making a table of information for participant
sample_info <- flow_mon_rna_clean %>%
  distinct(participant_id, rna_dataset, value_clean) %>%
  as.data.frame()

rownames(sample_info) <- sample_info$participant_id
sample_info <- sample_info[colnames(count_table), ]


#conducting vst
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData   = sample_info,
                              design    = ~ 1)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # remove genes with almost no reads
vsd <- vst(dds, blind = TRUE)

###PCA and Adjusting for Batch effects

#Run PCA
plotPCA(vsd, intgroup = "rna_dataset")   # colored by study
plotPCA(vsd, intgroup = "value_clean")   # colored by day 7 monocytes
saveRDS(vsd, "vsd.rds")                  # save so you don't have to redo it


#adjusting for batch effects and sex
sex_info <- dplyr::distinct(dem_clean, participant_id, biological_sex)
vsd$sex  <- factor(sex_info$biological_sex[match(colnames(vsd), sex_info$participant_id)])

vsd_adj <- vsd
assay(vsd_adj) <- removeBatchEffect(assay(vsd),
                                    batch  = vsd$rna_dataset,
                                    batch2 = vsd$sex)

plotPCA(vsd_adj, intgroup = "value_clean")   # colored by day 7 monocytes
plotPCA(vsd_adj, intgroup = "sex")           # check: should be one mixed cloud



###Gene Expression Anlaysis with voom-lima##

#install edge
library(edgeR)  

# Counts and normalization factors
d0 <- DGEList(count_table)
d0 <- calcNormFactors(d0)

# Filter low-expressed genes (automatic process , based library sizes and design)
keep <- filterByExpr(d0, mm)
d <- d0[keep, , keep.lib.sizes = FALSE]
y <- voom(d, mm, plot = TRUE)



#Sample info, in the same order as the count columns
ids <- colnames(d)

study <- factor(vsd$rna_dataset[match(ids, colnames(vsd))])
sex   <- factor(vsd$sex[match(ids, colnames(vsd))])
mono  <- vsd$value_clean[match(ids, colnames(vsd))]
age   <- as.numeric(dem_clean$age[match(ids, dem_clean$participant_id)])
race  <- factor(dem_clean$race[match(ids, dem_clean$participant_id)])


# Model 1: study + sex + monocytes (like the tutorial's batch + continuous examples)
mm1 <- model.matrix(~ study + sex + mono)

y1   <- voom(d, mm1, plot = TRUE)
fit1 <- lmFit(y1, mm1)
tmp1 <- contrasts.fit(fit1, coef = "mono")
tmp1 <- eBayes(tmp1)
res1 <- topTable(tmp1, sort.by = "P", n = Inf)

head(res1, 20)
length(which(res1$adj.P.Val < 0.05))
write.csv(res1, "limma_model1_study_sex.csv")

#Model 2: study + sex + age + monocytes
mm2 <- model.matrix(~ study + sex + age + mono)

y2   <- voom(d, mm2, plot = TRUE)
fit2 <- lmFit(y2, mm2)
tmp2 <- contrasts.fit(fit2, coef = "mono")
tmp2 <- eBayes(tmp2)
res2 <- topTable(tmp2, sort.by = "P", n = Inf)

head(res2, 20)
length(which(res2$adj.P.Val < 0.05))
write.csv(res2, "limma_model2_study_sex_age.csv")


#Model 3: study + sex + age + race + monocytes
mm3 <- model.matrix(~ study + sex + age + race + mono)
y3   <- voom(d, mm3, plot = TRUE)
fit3 <- lmFit(y3, mm3)
tmp3 <- contrasts.fit(fit3, coef = "mono")
tmp3 <- eBayes(tmp3)
res3 <- topTable(tmp3, sort.by = "P", n = Inf)

head(res3, 20)
length(which(res3$adj.P.Val < 0.05))
write.csv(res3, "limma_model3_study_sex_age_race.csv")

###Comparing 3 models
#Number of significant genes within each model

data.frame(
  model = c("1: study + sex", "2: + age", "3: + age + race"),
  n_significant = c(sum(res1$adj.P.Val < 0.05),
                    sum(res2$adj.P.Val < 0.05),
                    sum(res3$adj.P.Val < 0.05))
)

