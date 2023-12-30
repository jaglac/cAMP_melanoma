library(readr)
library(AnnotationHub)
library(tidyverse)
library(DESeq2)

################# set up 
anno.hub <- AnnotationHub()
##search for record: query(anno.hub, c("OrgDb","Homo sapiens"))
anno.hs <- anno.hub[["AH114084"]]
#example - df$symbol <- mapIds(anno.hs, keys=row.names(df),keytype="ENSEMBL", column="SYMBOL")

dataset_name <- "integrated_braf"  ### CHANGE THIS DEPENDING ON WHAT DATASET TO ANALYZE - for integrated set also change dds design matrix to ~source
######################## import processed data and annotation tables

metadata <- read.csv(file.path("tidy_input", dataset_name, "metadata.csv"), row.names = 1, header = TRUE)
counts <- read.csv(file.path("tidy_input", dataset_name, "roundedGeneCounts.csv"), row.names = 1, header = TRUE)
counts_dds <- counts
counts_dds$gene_ensembl <- NULL

#import processed annotated geneset dataframes
anno_files_list <- list.files(path="tidy_input/annotations", pattern = "csv")
tf_files_list <- list.files(path="tidy_input/annotations/individual_tfs")
aggr_files_list <- list.files(path="tidy_input/annotations/aggregated_tables")

anno <- list()
anno$tf_targets <- list()
anno$aggr <- list()

for (i in 1:length(anno_files_list)){
  name <- sub(".csv", "", anno_files_list[[i]])
  anno[[name]] <- read.csv(file.path("tidy_input/annotations", anno_files_list[[i]]), row.names = 1, header = TRUE)
}

for (i in 1:length(tf_files_list)){
  name <- sub("_targets.csv", "", tf_files_list[[i]])
  anno$tf_targets[[name]] <- read.csv(file.path("tidy_input/annotations/individual_tfs", tf_files_list[[i]]), row.names = 1, header = TRUE)
}

for (i in 1:length(aggr_files_list)){
  name <- sub(".csv", "", aggr_files_list[[i]])
  anno$aggr[[name]] <- read.csv(file.path("tidy_input/annotations/aggregated_tables", aggr_files_list[[i]]), row.names = 1, header = TRUE)
}

remove(anno_files_list, tf_files_list, i, name, aggr_files_list)

####################################### get log-space transformed expression values and z scores by gene (non-batch corrected)

df <- list()
df$dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_dds,
                                         colData = metadata,
                                         design = ~ source)
keep <- rowSums(counts(df$dds)) >= 50
df$dds <- df$dds[keep,]
remove(keep)

#without replicates for each cell line and using cell line as the model matrix design (which is necessary because that is what i want to compare without biasing based on previous subtype assignments by tsoi), must use blind=TRUE since it can't estimate dispersions by replicate group as it would if blind=FALSE
df$vst <- DESeq2::vst(df$dds, blind=TRUE) 

#build vst dataframe 
df$vst_df <- as.data.frame(assay(df$vst))


#get per-gene z scores
t <- as.data.frame(t(df$vst_df))
df$vst_df_z <- as.data.frame(t(sapply(t, function(t) (t-mean(t))/sd(t))))
rownames(df$vst_df_z) <- rownames(df$vst_df)
colnames(df$vst_df_z) <- colnames(df$vst_df)

remove(t)

df$vst_zscores_goi_all <- df$vst_df_z %>% filter(rownames(df$vst_df_z) %in% anno$aggr$aggr_all$gene_ensembl)
df$vst_zscores_goi_noTFtargets <- df$vst_df_z %>% filter(rownames(df$vst_df_z) %in% anno$aggr$aggr_without_tf_targets$gene_ensembl)

####################################### batch correct vst (log-space transformed) expression values and get z scores by gene

library(sva)

adjusted_counts <- as.matrix(counts_dds)
adjusted_counts <- ComBat_seq(adjusted_counts, batch=metadata$source, group=NULL)



df_corr <- list()

df_corr$dds <- DESeq2::DESeqDataSetFromMatrix(countData = adjusted_counts,
                                              colData = metadata,
                                              design = ~ source)

keep <- rowSums(counts(df_corr$dds)) >= 50
df_corr$dds <- df_corr$dds[keep,]
remove(keep)

#without replicates for each cell line and using cell line as the model matrix design (which is necessary because that is what i want to compare without biasing based on previous subtype assignments by tsoi), must use blind=TRUE since it can't estimate dispersions by replicate group as it would if blind=FALSE
df_corr$vst <- DESeq2::vst(df_corr$dds, blind=TRUE) 

#build vst dataframe 
df_corr$vst_df <- as.data.frame(assay(df_corr$vst))

#get per-gene z scores
t <- as.data.frame(t(df_corr$vst_df))
df_corr$vst_df_z <- as.data.frame(t(sapply(t, function(t) (t-mean(t))/sd(t))))
rownames(df_corr$vst_df_z) <- rownames(df_corr$vst_df)
colnames(df_corr$vst_df_z) <- colnames(df_corr$vst_df)

remove(t)

df_corr$vst_zscores_goi_all <- df_corr$vst_df_z %>% filter(rownames(df_corr$vst_df_z) %in% anno$aggr$aggr_all$gene_ensembl)
df_corr$vst_zscores_goi_noTFtargets <- df_corr$vst_df_z %>% filter(rownames(df_corr$vst_df_z) %in% anno$aggr$aggr_without_tf_targets$gene_ensembl)

###################################

dir.create(file.path("output/"))
dir.create(file.path("output", dataset_name))
dir.create(file.path("output", dataset_name, "non_batch_corrected"))
write.csv(df$vst_df, file = file.path("output", dataset_name, "non_batch_corrected" , "vst.csv"), quote = TRUE)
write.csv(df$vst_df_z, file = file.path("output", dataset_name, "non_batch_corrected", "vst_zscores.csv"), quote = TRUE)
write.csv(df$vst_zscores_goi_noTFtargets, file = file.path("output", dataset_name, "non_batch_corrected", "vst_zscores_goi_noTFtargets.csv"), quote = TRUE)
write.csv(df$vst_zscores_goi_all, file = file.path("output", dataset_name, "non_batch_corrected", "vst_zscores_goi_all.csv"), quote = TRUE)

write.csv(df_corr$vst_df, file = file.path("output", dataset_name, "vst.csv"), quote = TRUE)
write.csv(df_corr$vst_df_z, file = file.path("output", dataset_name, "vst_zscores.csv"), quote = TRUE)
write.csv(df_corr$vst_zscores_goi_noTFtargets, file = file.path("output", dataset_name, "vst_zscores_goi_noTFtargets.csv"), quote = TRUE)
write.csv(df_corr$vst_zscores_goi_all, file = file.path("output", dataset_name, "vst_zscores_goi_all.csv"), quote = TRUE)

rm(list = ls())
gc()    


