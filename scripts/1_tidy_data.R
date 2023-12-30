library(BiocManager)
library(readr)
library(AnnotationHub)
library(tidyverse)


anno.hub <- AnnotationHub()
##search for record: query(anno.hub, c("OrgDb","Homo sapiens"))
anno.hs <- anno.hub[["AH114084"]]
#example - df$symbol <- mapIds(anno.hs, keys=row.names(df),keytype="ENSEMBL", column="SYMBOL")



############### CLEAN COUNTS DATA (TSOI 53 CELL LINES) 

tsoi_all_metadata <- read_table("input_raw/tsoi_basal/grein_GSE80829_filtered_metadata_clean.txt")
table_s1 <- read_table("input_raw/tsoi_supplemental_tables/table_s1_cellLineInfo_clean.txt")
tsoi_all_metadata <- merge(tsoi_all_metadata, table_s1, by = "cell_line", all = TRUE)
tsoi_all_metadata <- tsoi_all_metadata %>% arrange(sample)
rownames(tsoi_all_metadata) <- tsoi_all_metadata$sample
remove(table_s1)

tsoi_all_counts <- read_csv("input_raw/tsoi_basal/grein_GSE80829_GeneLevel_Raw_data.csv")
tsoi_all_counts <- tsoi_all_counts %>% column_to_rownames(var="...1")
tsoi_all_counts <- tsoi_all_counts[,2:ncol(tsoi_all_counts)] %>% dplyr::select(sort(tidyselect::peek_vars()))
tsoi_all_counts <- round(tsoi_all_counts, digits = 0)
tsoi_all_counts <- tsoi_all_counts %>% add_column(gene_ensembl = rownames(tsoi_all_counts), .before = 1) 

### FILTER TO ONLY BRAF-MUTANT CELL LINES
tsoi_braf_metadata <- tsoi_all_metadata %>% filter(mut_status == "BRAF")
tsoi_braf_counts <- tsoi_all_counts %>% dplyr::select(any_of(rownames(tsoi_braf_metadata)))


### tsoi drug treatment timecourse
tsoi_dt_metadata <- read_csv("input_raw/tsoi_drug_timecourse/grein_GSE110054_filtered_metadata_clean.csv")
tsoi_dt_counts <- read_csv("input_raw/tsoi_drug_timecourse/grein_GSE110054_GeneLevel_Raw_data.csv")
tsoi_dt_counts <- tsoi_dt_counts %>% column_to_rownames(var="...1")
tsoi_dt_counts <- tsoi_dt_counts[,2:ncol(tsoi_dt_counts)] %>% dplyr::select(sort(tidyselect::peek_vars()))
tsoi_dt_counts <- round(tsoi_dt_counts, digits = 0)
tsoi_dt_counts <- tsoi_dt_counts %>% add_column(gene_ensembl = rownames(tsoi_dt_counts), .before = 1) 


### ccle data
ccle_braf_metadata <- read_delim("input_raw/ccle_cell_lines/BRAF_Mutant_Melanoma_Cell_Lines.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
colnames(ccle_braf_metadata) <- c("cell_line", "ccle_id", "depmap_broad_id", "histology", "braf_mutation")
ccle_braf_metadata$safe_id <- str_replace(ccle_braf_metadata$depmap_broad_id, "-", "")


ctrp_comb_ccle_diffstate <- read_delim("input_raw/ccle_cell_lines/CTRP combined melanoma CCLE cell line predictions Tsoi Graeber et al.txt", delim = "\t", escape_double = FALSE, col_types = cols(...5 = col_skip()), 
                                       trim_ws = TRUE) %>% dplyr::mutate(Final = tolower(Final),
                                                                         CCLE = tolower(CCLE),
                                                                         GDSC = tolower(GDSC))
colnames(ctrp_comb_ccle_diffstate) <- c("cell_line", "gdsc_diffstate", "ccle_diffstate", "comb_diffstate")
ctrp_comb_ccle_diffstate$comb_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$comb_diffstate, 
                                                  pattern = fixed(" ("), replacement = "_")
ctrp_comb_ccle_diffstate$comb_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$comb_diffstate, 
                                                  pattern = fixed(")"), replacement = "")
ctrp_comb_ccle_diffstate$comb_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$comb_diffstate, 
                                                  pattern = fixed("-"), replacement = "_")
ctrp_comb_ccle_diffstate$comb_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$comb_diffstate, 
                                                  pattern = fixed(" "), replacement = "_")
ctrp_comb_ccle_diffstate$ccle_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$ccle_diffstate, 
                                                           pattern = fixed("-"), replacement = "_")
ctrp_comb_ccle_diffstate$ccle_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$ccle_diffstate, 
                                                           pattern = fixed(" "), replacement = "_")
ctrp_comb_ccle_diffstate$gdsc_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$gdsc_diffstate, 
                                                           pattern = fixed("-"), replacement = "_")
ctrp_comb_ccle_diffstate$gdsc_diffstate <- str_replace_all(ctrp_comb_ccle_diffstate$gdsc_diffstate, 
                                                           pattern = fixed(" "), replacement = "_")


ccle_braf_metadata <- left_join(ccle_braf_metadata, ctrp_comb_ccle_diffstate)
ccle_braf_metadata <- ccle_braf_metadata %>% arrange(depmap_broad_id)



ccle_counts <- read_csv("input_raw/ccle_cell_lines/CCLE_expression_proteincoding_genes_expected_count.csv", 
                        trim_ws = FALSE)
colnames(ccle_counts) <- c("sample", colnames(ccle_counts[,2:ncol(ccle_counts)]))
ccle_braf_counts <- ccle_counts %>% filter(sample %in% ccle_braf_metadata$depmap_broad_id)
ccle_braf_counts <- t(ccle_braf_counts)
colnames(ccle_braf_counts) <- ccle_braf_counts[1,]
ccle_braf_counts <- ccle_braf_counts[2:nrow(ccle_braf_counts),] 
ccle_braf_counts <- as.data.frame(ccle_braf_counts) %>% dplyr::select(sort(tidyselect::peek_vars()))
colnames(ccle_braf_counts) <- ccle_braf_metadata$safe_id
ccle_braf_counts[] <- sapply(ccle_braf_counts, as.numeric)
ccle_braf_counts[] <- round(ccle_braf_counts, digits=0)
ccle_braf_counts <- ccle_braf_counts %>% rownames_to_column(var="gene_info")
ccle_braf_counts <- ccle_braf_counts %>% separate(gene_info, into = c("gene_symbol", "gene_entrez"), sep = " ", remove=FALSE) %>% 
                                         mutate(gene_entrez = str_replace_all(string = gene_entrez, pattern = "[\\(\\)]", replacement = "")) 
ccle_braf_counts <- ccle_braf_counts %>% add_column(gene_ensembl = mapIds(anno.hs, keys=ccle_braf_counts$gene_entrez, 
                                                                      keytype="ENTREZID", column="ENSEMBL"), .after = "gene_entrez") 
ccle_braf_counts <- ccle_braf_counts %>% filter(!is.na(gene_ensembl) & !duplicated(.[["gene_ensembl"]]))
rownames(ccle_braf_counts) <- ccle_braf_counts$gene_ensembl
ccle_braf_counts <- ccle_braf_counts[, 4:ncol(ccle_braf_counts)]
  
########################
#integrate metadata and counts tables 

integrate_ccle <- ccle_braf_metadata %>% mutate(mut_status = "BRAF",
                                                    sample = safe_id,
                                                    subtype = NA, 
                                                    treatment = NA, #not actually dmso just need to match controls from tsoi_dt
                                                    dose_nM = NA,
                                                    time_days = NA, 
                                                    source = "ccle") %>% 
                                          dplyr::select(sample, cell_line, mut_status, subtype, treatment, dose_nM, time_days, source, contains("diffstate"))
integrate_tsoi <- tsoi_all_metadata %>% mutate(treatment = NA, #not actually dmso just need to match controls from tsoi_dt
                                               dose_nM = NA,
                                               time_days = NA, 
                                               source = "tsoi_basal") %>%
                                        dplyr::select(sample, cell_line, mut_status, subtype, treatment, dose_nM, time_days, source)
integrate_tsoi_dt <- tsoi_dt_metadata %>% mutate(dose_nM = as.character(dose_nM),
                                                 time_days = as.character(time_days),
                                                 source = "tsoi_dt") %>% 
                                        dplyr::select(sample, cell_line, mut_status, subtype, treatment, dose_nM, time_days, source)
integrate_metadata <- bind_rows(integrate_ccle, integrate_tsoi, integrate_tsoi_dt)
rownames(integrate_metadata) <- NULL


integrate_counts <- inner_join(ccle_braf_counts, tsoi_all_counts) %>% inner_join(tsoi_dt_counts)
integrate_counts <- integrate_counts %>% mutate(sum = rowSums(integrate_counts[,2:ncol(integrate_counts)]),
                                 mean = rowMeans(integrate_counts[,2:ncol(integrate_counts)])) %>%
                          filter(sum > 100 & mean > 5) %>% dplyr::select(-c(sum, mean))
integrate_counts <- integrate_counts %>% column_to_rownames(var="gene_ensembl")
integrate_counts <- integrate_counts %>% add_column(gene_ensembl = rownames(integrate_counts), .before=1)

integrate_metadata_braf <- integrate_metadata %>% filter(mut_status == "BRAF" & source != "tsoi_dt") %>% 
                                                  dplyr::select(sample, cell_line, subtype, source, contains("diffstate")) %>% 
                                                  mutate(tsoi_diffstate = subtype, .keep="unused")  %>% replace_na(list(gdsc_diffstate="na_nd",
                                                                                                                        ccle_diffstate="na_nd",
                                                                                                                        comb_diffstate="na_nd",
                                                                                                                        tsoi_diffstate="na_nd"))
integrate_counts_braf <- integrate_counts%>% dplyr::select(any_of(c("gene_ensembl", integrate_metadata_braf$sample)))





dir.create(file.path("tidy_input/"))
dir.create(file.path("tidy_input/tsoi_basal"))
write.csv(tsoi_all_metadata, file = "tidy_input/tsoi_basal/metadata.csv", quote=TRUE)
write.csv(tsoi_all_counts, file = "tidy_input/tsoi_basal/roundedGeneCounts.csv", quote=TRUE)

dir.create(file.path("tidy_input/tsoi_basal_braf"))
write.csv(tsoi_braf_metadata, file = "tidy_input/tsoi_basal_braf/metadata.csv", quote=TRUE)
write.csv(tsoi_braf_counts, file = "tidy_input/tsoi_basal_braf/roundedGeneCounts.csv", quote=TRUE)

dir.create(file.path("tidy_input/tsoi_dt"))
write.csv(tsoi_dt_metadata, file = "tidy_input/tsoi_dt/metadata.csv", quote=TRUE)
write.csv(tsoi_dt_counts, file = "tidy_input/tsoi_dt/roundedGeneCounts.csv", quote=TRUE)

dir.create(file.path("tidy_input/ccle"))
write.csv(ccle_braf_counts, file = "tidy_input/ccle/roundedGeneCounts.csv", quote=TRUE)
write.csv(ccle_braf_metadata, file = "tidy_input/ccle/metadata.csv", quote=TRUE)

dir.create(file.path("tidy_input/integrated"))
write.csv(integrate_metadata, file = "tidy_input/integrated/metadata.csv", quote=TRUE)
write.csv(integrate_counts, file = "tidy_input/integrated/roundedGeneCounts.csv", quote=TRUE)

dir.create(file.path("tidy_input/integrated_braf"))
write.csv(integrate_metadata_braf, file = "tidy_input/integrated_braf/metadata.csv", quote=TRUE)
write.csv(integrate_counts_braf, file = "tidy_input/integrated_braf/roundedGeneCounts.csv", quote=TRUE)


rm(list = ls())
gc()


