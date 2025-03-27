library(DESeq2)
library(tidyverse)
library(dplyr)

# load file
bladder_data <- readRDS("/Users/emilyzhang/bladder_cancer_classifier/data/UROMOL_TaLG.teachingcohort.rds")

recurrence_1 <- bladder_data[bladder_data$Recurrence == 1, ]
recurrence_0 <- bladder_data[bladder_data$Recurrence == 0, ]

count_data <- cbind(recurrence_1$exprs[1:nrow(recurrence_0$exprs), ], recurrence_0$exprs)


dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)

# 5. Run analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Recurrence", "No_Recurrence"))

# 6. Get significant genes (padj < 0.05)
sig_genes <- rownames(res)[which(res$padj < 0.05)]

# 7. Create filtered dataset
filtered_exprs <- count_data[sig_genes, ]
clinical_data <- bladder_data %>% select(-exprs)
final_data <- cbind(clinical_data, t(filtered_exprs))


# 8. Save results
write_csv(final_data, "recurrence_filtered_data.csv")