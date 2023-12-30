library(readr)
library(tidyverse)

################# set up 
anno.hub <- AnnotationHub::AnnotationHub()
##search for record: query(anno.hub, c("OrgDb","Homo sapiens"))
anno.hs <- anno.hub[["AH114084"]]
#example - df$symbol <- mapIds(anno.hs, keys=row.names(df),keytype="ENSEMBL", column="SYMBOL")

dataset_name <- "integrated_braf"
#########################################################
df <- list()
df$vst_zscores <- read.csv(file.path("output", dataset_name, "vst_zscores.csv"), row.names = 1, header = TRUE)
df$vst_zscores_raw <- df$vst_zscores
df$vst_zscores <- df$vst_zscores %>% add_column(gene_ensembl = rownames(df$vst_zscores), .before=1)
  
  
  
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

#get mean z score for genesets
df$bioplanet_gpcr_mean_z <- data.frame(sample = colnames(df$vst_zscores[,2:ncol(df$vst_zscores)]))

for (i in unique(anno$bioplanet_gpcr$bioplanet_id)) {
  list <- anno$bioplanet_gpcr %>% filter(bioplanet_id == i) 
  z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
  z_avg <- data.frame(name = colMeans(z[ ,2:ncol(z)]))
  df$bioplanet_gpcr_mean_z <- cbind(df$bioplanet_gpcr_mean_z, z_avg)
  remove(list,z,z_avg)
}
colnames(df$bioplanet_gpcr_mean_z) <- c("sample", unique(anno$bioplanet_gpcr$bioplanet_id))
t <- df$bioplanet_gpcr_mean_z %>% dplyr::select(!contains("sample")) %>% dplyr::select(sort(tidyselect::peek_vars()))
colnames(t) <- anno$bioplanet_key$new_col_name
df$bioplanet_gpcr_mean_z <- bind_cols(df$bioplanet_gpcr_mean_z %>% dplyr::select(sample),t)
remove(t)

#get mean z score for GO gpcr biological process genesets
df$go_all_mean_z <- data.frame(sample = colnames(df$vst_zscores[,2:ncol(df$vst_zscores)]))
for (i in unique(anno$go_all$go_id_safe)) {
  list <- anno$go_all %>% filter(go_id_safe == i) 
  z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
  z_avg <- data.frame(name = colMeans(z[ ,2:ncol(z)]))
  df$go_all_mean_z <- cbind(df$go_all_mean_z, z_avg)
  remove(list,z,z_avg)
}
go_key <- read.csv("tidy_input/annotations/go_key.csv", row.names = 1, header = TRUE)
colnames(df$go_all_mean_z) <- c("sample", anno$go_key$new_col_names)

######### individual transcription factor targets

tf_list <- names(anno$tf_targets)
df$tf_targets_mean_z <- data.frame(sample = colnames(df$vst_zscores[,2:ncol(df$vst_zscores)]))
rownames(df$tf_targets_mean_z) <- df$tf_targets_mean_z$sample

for (i in tf_list){
  z <- df$vst_zscores %>% filter(gene_ensembl %in% anno$tf_targets[[i]]$gene_ensembl)
  z_avg <- data.frame(colMeans(z[ ,2:ncol(z)]))
  df$tf_targets_mean_z <- bind_cols(df$tf_targets_mean_z, z_avg)
}
colnames(df$tf_targets_mean_z) <- c("sample", paste0(tf_list, "_targets"))
remove(i,tf_list, z, z_avg)

########## tsoi subtype signatures
df$tsoi_subtype_signatures_mean_z <- data.frame(sample = colnames(df$vst_zscores[,2:ncol(df$vst_zscores)]))
subtype_signature <- unique(anno$tsoi_subtype_signatures$subtype_signature)
for (i in unique(anno$tsoi_subtype_signatures$subtype_signature)) { 
  list <- anno$tsoi_subtype_signatures %>% filter(subtype_signature == i)
  z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
  z_avg <- data.frame(name = colMeans(z[ ,2:ncol(z)]))
  df$tsoi_subtype_signatures_mean_z <- cbind(df$tsoi_subtype_signatures_mean_z, z_avg)
  remove(list,z,z_avg)
}
colnames(df$tsoi_subtype_signatures_mean_z) <- c("sample", paste0(subtype_signature, "_signature"))
remove(i,subtype_signature)



######### known genes of interest (select rows from vst gene expression z score df)
### all genes, no summarizing to mean z
gene_z <- df$vst_zscores %>% filter(gene_ensembl %in% anno$known_genes_of_interest$gene_ensembl)

rownames(gene_z) <- mapIds(anno.hs, keys = gene_z$gene_ensembl, keytype = "ENSEMBL", column = "SYMBOL")

gene_z <- gene_z[ , 2:ncol(gene_z)] %>% t() %>% as.data.frame() %>% dplyr::select(order(colnames(.)))

###CREB family TFs
list <- anno$known_genes_of_interest %>% filter(characteristic == "CREB_family_tf") %>% dplyr::select(gene_ensembl)
creb_family_z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
creb_family_z_avg <- data.frame(creb_family_tfs = colMeans(creb_family_z[ ,2:ncol(creb_family_z)]))
remove(creb_family_z)

###CREB binding cofactors
list <- anno$known_genes_of_interest %>% filter(characteristic == "CREB_binding_cofactor") %>% dplyr::select(gene_ensembl)
creb_binding_z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
creb_binding_z_avg <- data.frame(creb_binding_cofactors = colMeans(creb_binding_z[ ,2:ncol(creb_binding_z)]))
remove(creb_binding_z)

###AP1 family TFs
list <- anno$known_genes_of_interest %>% filter(characteristic == "AP1_family_tf") %>% dplyr::select(gene_ensembl)
ap1_family_z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
ap1_family_z_avg <- data.frame(ap1_family_tfs = colMeans(ap1_family_z[ ,2:ncol(ap1_family_z)]))
remove(ap1_family_z)

###CREB activity markers
list <- anno$known_genes_of_interest %>% filter(grepl("marker_CREB_activity", characteristic)) %>% dplyr::select(gene_ensembl)
creb_marker_z <- df$vst_zscores %>% filter(gene_ensembl %in% list$gene_ensembl)
creb_marker_z_avg <- data.frame(creb_activity_markers = colMeans(creb_marker_z[ ,2:ncol(creb_marker_z)]))
remove(creb_marker_z)

df$known_genes_of_interest_z <- bind_cols(creb_family_z_avg, creb_binding_z_avg, creb_marker_z_avg, ap1_family_z_avg, gene_z)
remove(list,creb_family_z_avg, creb_binding_z_avg, creb_marker_z_avg, ap1_family_z_avg, gene_z)

################################ BIND ALL GENESET Z SCORES INTO ONE DF
df$bind_zscores <- bind_cols(df[c("tsoi_subtype_signatures_mean_z", "known_genes_of_interest_z", 
                                              "go_all_mean_z", "bioplanet_gpcr_mean_z", "tf_targets_mean_z")]) %>% 
                                dplyr::select(!starts_with("sample")) 

################################ 

dir.create(file.path("output", dataset_name))
dir.create(file.path("output", dataset_name, "mean_zscores"))

write.csv(df$bioplanet_gpcr_mean_z, file = file.path("output", dataset_name, "mean_zscores", "bioplanet_gpcr_z.csv"), quote = TRUE)
write.csv(df$tsoi_subtype_signatures_mean_z, file = file.path("output", dataset_name, "mean_zscores", "tsoi_subtype_sig_z.csv"), quote = TRUE)
write.csv(df$known_genes_of_interest_z, file = file.path("output", dataset_name, "mean_zscores", "known_goi_z.csv"), quote = TRUE)
write.csv(df$tf_targets_mean_z, file = file.path("output", dataset_name, "mean_zscores", "tf_targets_z.csv"), quote = TRUE)
write.csv(df$go_all_mean_z, file = file.path("output", dataset_name, "mean_zscores", "go_all_z.csv"), quote = TRUE)
write.csv(df$bind_zscores, file = file.path("output", dataset_name, "mean_zscores", "all_genesets_z.csv"), quote = TRUE)



rm(list = ls())
gc()    

