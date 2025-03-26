library(DESeq2)
library(tidyverse)

# load file
bladder_data <- readRDS("./data/UROMOL_TaLG.teachingcohort.rds")

# separate expression data
expression_columns <- grep("^expr\\.", names(bladder_data), value = TRUE)
clinical_data <- bladder_data %>% select(-all_of(expression_columns))

# make deseq object
count_matrix <- bladder_data %>%
  select(all_of(expression_columns)) %>%
  t() %>%
  round() # rounding for deseq2 to work

coldata <- clinical_data %>%
  mutate(recurrence_rate = factor(recurrence_rate)) %>%
  column_to_rownames("patient_id")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ recurrence_rate
)

# running DEA
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)

# identify significant genes
significant_genes <- results %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames_to_column("gene") %>%
  pull(gene)

# create dataset with only differentially expressed genes
filtered_data <- bladder_data %>%
  select(
    all_of(names(clinical_data)),
    all_of(paste0("expr.", significant_genes))
  )

# save 
write_csv(filtered_data, "./data/UROMOL_TaLG.teachingcohort_filted.csv")
