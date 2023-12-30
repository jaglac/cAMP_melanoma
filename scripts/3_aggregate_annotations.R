library(tidyverse)
library(readr)
library(AnnotationHub)
anno.hub <- AnnotationHub()
##search for record: query(anno.hub, c("OrgDb","Homo sapiens"))
anno.hs <- anno.hub[["AH114084"]]
#example - df$symbol <- mapIds(anno.hs, keys=row.names(df),keytype="ENSEMBL", column="SYMBOL")


anno <- list()
anno$tf_targets <- list()

#import processed annotated geneset dataframes
anno_files_list <- list.files(path="tidy_input/annotations", pattern = "csv")
tf_files_list <- list.files(path="tidy_input/annotations/individual_tfs")

for (i in 1:length(anno_files_list)){
  name <- sub(".csv", "", anno_files_list[[i]])
  anno[[name]] <- read.csv(file.path("tidy_input/annotations", anno_files_list[[i]]), row.names = 1, header = TRUE)
}

for (i in 1:length(tf_files_list)){
  name <- sub("_targets.csv", "", tf_files_list[[i]])
  anno$tf_targets[[name]] <- read.csv(file.path("tidy_input/annotations/individual_tfs", tf_files_list[[i]]), row.names = 1, header = TRUE)
}

remove(anno_files_list, tf_files_list, i, name)


anno$aggr <- list()
anno$aggr$go <- dplyr::full_join(aggregate(go_id ~ gene_ensembl, 
                                          data = anno$go_all, 
                                          paste, collapse = "; "),
                                  aggregate(go_category ~ gene_ensembl, 
                                            data = anno$go_all, 
                                            paste, collapse = "; "), 
                                  by = "gene_ensembl") %>%
            dplyr::full_join(aggregate(go_term_name ~ gene_ensembl, 
                             data = anno$go_all, 
                             paste, collapse = "; "), by = "gene_ensembl")

anno$aggr$bp <- merge(aggregate(bioplanet_id ~ gene_ensembl, 
                          data = anno$bioplanet_gpcr, 
                          paste, collapse = "; "),
                aggregate(bioplanet_name ~ gene_ensembl, 
                          data =anno$bioplanet_gpcr, 
                          paste, collapse = "; "), 
                by = "gene_ensembl")

anno$aggr$tf_targets <- dplyr::full_join(aggregate(tf_ensembl ~ gene_ensembl, 
                                               data = bind_rows(anno$tf_targets), 
                                               paste, collapse = "; "),
                                     aggregate(tf_symbol ~ gene_ensembl, 
                                               data = bind_rows(anno$tf_targets), 
                                               paste, collapse = "; "), 
                                     by = "gene_ensembl") %>% 
                  dplyr::full_join(aggregate(tf_entrez ~ gene_ensembl, 
                                     data = bind_rows(anno$tf_targets), 
                                     paste, collapse = "; "),by = "gene_ensembl") %>%
                  dplyr::full_join(aggregate(tf_name ~ gene_ensembl, 
                                   data = bind_rows(anno$tf_targets), 
                                   paste, collapse = "; "),by = "gene_ensembl")
colnames(anno$aggr$tf_targets) <- c("gene_ensembl", colnames(anno$aggr$tf_targets[,2:ncol(anno$aggr$tf_targets)]))

anno$aggr$no_tf_targets <- dplyr::full_join(anno$known_genes_of_interest[ , c("gene_ensembl", "characteristic")], 
                                  anno$tsoi_subtype_signatures[ , c("gene_ensembl", "subtype_signature")], by = "gene_ensembl") %>%
                                dplyr::full_join(anno$aggr$go, by = "gene_ensembl") %>% 
                                dplyr::full_join(anno$aggr$bp, by = "gene_ensembl") 
anno$aggr$no_tf_targets <- anno$aggr$no_tf_targets %>% 
                              add_column(gene_symbol = 
                                mapIds(anno.hs, keys=anno$aggr$no_tf_targets$gene_ensembl, keytype="ENSEMBL", column="SYMBOL"), 
                                .after = "gene_ensembl")
anno$aggr$no_tf_targets <- anno$aggr$no_tf_targets %>% 
                              add_column(gene_name = 
                                           mapIds(anno.hs, keys=anno$aggr$no_tf_targets$gene_ensembl, keytype="ENSEMBL", column="GENENAME"), 
                                         .after = "gene_symbol")




anno$aggr$all <- dplyr::full_join(anno$known_genes_of_interest[ , c("gene_ensembl", "characteristic")], 
                              anno$tsoi_subtype_signatures[ , c("gene_ensembl", "subtype_signature")], by = "gene_ensembl") %>%
               dplyr::full_join(anno$aggr$go, by = "gene_ensembl") %>% 
               dplyr::full_join(anno$aggr$bp, by = "gene_ensembl") %>%
               dplyr::full_join(anno$aggr$tf_targets, by = "gene_ensembl")

anno$aggr$all <- anno$aggr$all %>% 
  add_column(gene_symbol = 
               mapIds(anno.hs, keys=anno$aggr$all$gene_ensembl, keytype="ENSEMBL", column="SYMBOL"), 
             .after = "gene_ensembl")
anno$aggr$all <- anno$aggr$all %>% 
  add_column(gene_name = 
               mapIds(anno.hs, keys=anno$aggr$all$gene_ensembl, keytype="ENSEMBL", column="GENENAME"), 
             .after = "gene_symbol")

################
dir.create(file.path("tidy_input/annotations/aggregated_tables"))

write.csv(anno$aggr$go, file = "tidy_input/annotations/aggregated_tables/aggr_go_terms.csv")
write.csv(anno$aggr$bp, file = "tidy_input/annotations/aggregated_tables/aggr_bioplanet.csv")
write.csv(anno$aggr$no_tf_targets, file = "tidy_input/annotations/aggregated_tables/aggr_without_tf_targets.csv")
write.csv(anno$aggr$all, file = "tidy_input/annotations/aggregated_tables/aggr_all.csv")


rm(list = ls())
gc()    
