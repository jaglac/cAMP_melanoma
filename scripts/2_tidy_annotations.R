library(readr)
library(tidyverse)
library(AnnotationHub)
library(GO.db)

############### CLEAN AND BUILD ANNOTATION TABLES/GENESETS 

anno.hub <- AnnotationHub()
##search for record: query(anno.hub, c("OrgDb","Homo sapiens"))
anno.hs <- anno.hub[["AH114084"]]
#example - df$symbol <- mapIds(anno.hs, keys=row.names(df),keytype="gene_ensembl", column="SYMBOL")



anno <- list()

### TSOI DIFFERENTIALLY EXPRESSED GENE SIGNATURES FOR DIFFERENTIATION STATE SUBTYPES
anno$tsoi_subtype_signatures <- read_csv("input_raw/tsoi_supplemental_tables/table_s3_subtypeSignatures_diffUpGenes.csv")
colnames(anno$tsoi_subtype_signatures) <- c("gene_symbol", "subtype_signature", "gene_name")
anno$tsoi_subtype_signatures$subtype_signature <- tolower(anno$tsoi_subtype_signatures$subtype_signature)
anno$tsoi_subtype_signatures$gene_ensembl <- mapIds(anno.hs, keys=anno$tsoi_subtype_signatures$gene_symbol, keytype="ALIAS", column="ENSEMBL")
anno$tsoi_subtype_signatures <- anno$tsoi_subtype_signatures %>% 
  mutate(subtype_signature = recode(subtype_signature, 
                                     "undifferentiated-neural crest-like" = "undifferentiated_neural_crest_like",
                                    "neural crest-like" = "neural_crest_like", 
                                    "neural crest-like-transitory" = "neural_crest_like_transitory"))




### KNOWN GENES OF INTEREST
marker_genes <- data.frame(gene_symbol = c("MITF", "SOX10", "AXL", "NGFR", 
                                           "MKI67", "RPS6KA1", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KB1", 
                                           "PCK1", "NR4A1", "NR4A2", "NR4A3"), 
                                characteristic = c("marker_melanocytic", "marker_differentiated", "marker_undifferentated", "marker_neural_crest_like",
                                                   "marker_proliferation", "marker_proliferation", "marker_proliferation", 
                                                   "marker_proliferation", "marker_proliferation", "marker_proliferation",
                                                   "marker_CREB_activity", "marker_CREB_activity", 
                                                   "marker_CREB_activity", "marker_CREB_activity"))
ap1_tf_genes <- data.frame(gene_symbol = c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),
                                characteristic = "AP1_family_tf")
creb_tf_genes <- data.frame(gene_symbol = c("CREB1", "CREB3", "CREB5", "CREBL2", "CREB3L1", "CREB3L2", "CREB3L3", "CREB3L4",
                                                 "CREM", "CREBZF", "CREBRF", "ATF1", "ATF2", "ATF4", "ATF6B"),
                                 characteristic = "CREB_family_tf")
creb_binding_genes <- data.frame(gene_symbol = c("CREBBP", "CRTC1", "CRTC2", "CRTC3", "EP300", "EID1", "PPARGC1A"),
                                 characteristic = "CREB_binding_cofactor")
anno$genes_of_interest <- bind_rows(marker_genes, ap1_tf_genes, creb_tf_genes, creb_binding_genes)
anno$genes_of_interest$gene_ensembl <- mapIds(anno.hs, keys=anno$genes_of_interest$gene_symbol, keytype="ALIAS", column="ENSEMBL")
anno$genes_of_interest$gene_name <- mapIds(anno.hs, keys=anno$genes_of_interest$gene_symbol, keytype="ALIAS", column="GENENAME")
remove(marker_genes, ap1_tf_genes, creb_tf_genes, creb_binding_genes)


### downloaded genesets from bioplanet website
anno$bioplanet <- read_delim("input_raw/annotations/bioplanet_pathway_genesets.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(GENE_ID = col_character()), 
                             trim_ws = TRUE)
colnames(anno$bioplanet) <- c("bioplanet_id", "bioplanet_name", tolower(colnames(anno$bioplanet[ , 3:4])))
anno$bioplanet$gene_ensembl <- mapIds(anno.hs, keys=anno$bioplanet$gene_id, keytype="ENTREZID", column="ENSEMBL")
anno$bioplanet$gene_name <- mapIds(anno.hs, keys=anno$bioplanet$gene_id, keytype="ENTREZID", column="GENENAME")
anno$bioplanet_gpcr <- anno$bioplanet %>% dplyr::filter(grepl(c("GPCR|G-protein|G alpha|G beta|cAMP|PKA|CREB|arrestin|adenylate cyclase"), bioplanet_name))
anno$bioplanet_gpcr_list <- unique(anno$bioplanet_gpcr[,1:2])
anno$bioplanet_key <- anno$bioplanet_gpcr %>% dplyr::select(contains("bioplanet")) %>% distinct() %>% mutate(new_col_name = paste0(bioplanet_name, " (", bioplanet_id, ")")) %>% mutate(new_col_name = str_replace(new_col_name, "G-protein coupled receptors", "GPCRs")) %>% mutate(new_col_name = str_replace(new_col_name, "bioplanet_", "bp_")) %>% arrange(bioplanet_id)


### GO:BIOLOGICAL PROCESS GPCR SIGNALING & SUB-SETS
anno$go_gpcr <- read_delim("input_raw/annotations/amiGO_GO0007186_GPCRsignaling.txt", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)[, c(1,2,4)]
colnames(anno$go_gpcr) <- c("uniprot", "protein_name", "go_id")
anno$go_gpcr <- distinct(anno$go_gpcr) %>% 
  tidyr::separate(uniprot, into = c(NA, "uniprot"), sep = ":") %>% 
  tidyr::separate(go_id, into = c(NA, "go_id_number"), sep = ":", remove=FALSE) %>% 
  group_by(go_id) %>% filter(n() >= 5)
anno$go_gpcr$go_id_safe <- as.character(paste0("GO_", anno$go_gpcr$go_id_number))
anno$go_gpcr$go_term_name <- mapIds(GO.db, keys=anno$go_gpcr$go_id, column="TERM", keytype="GOID")
anno$go_gpcr$gene_ensembl <- mapIds(anno.hs, keys=anno$go_gpcr$uniprot, keytype="UNIPROT", column="ENSEMBL")
anno$go_gpcr$gene_symbol <- mapIds(anno.hs, keys=anno$go_gpcr$uniprot, keytype="UNIPROT", column="SYMBOL")
anno$go_gpcr$gene_name <- mapIds(anno.hs, keys=anno$go_gpcr$uniprot, keytype="UNIPROT", column="GENENAME")
anno$go_gpcr$go_category <- "gpcr_signaling"
anno$go_gpcr <- anno$go_gpcr %>% filter( !is.na(gene_ensembl))
### GO:BIOLOGICAL PROCESS CELL CYCLE REGULATION (POSITIVE)
anno$go_posCellCycle <- read_delim("input_raw/annotations/amiGO_posRegOfCellCycle.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   col_names = FALSE, trim_ws = TRUE)
colnames(anno$go_posCellCycle) <- c("uniprot", "protein_name", "go_id")
anno$go_posCellCycle <- distinct(anno$go_posCellCycle) %>% 
  tidyr::separate(uniprot, into = c(NA, "uniprot"), sep = ":") %>% 
  tidyr::separate(go_id, into = c(NA, "go_id_number"), sep = ":", remove=FALSE) %>% 
  group_by(go_id) %>% filter(n() >= 5)
anno$go_posCellCycle$go_id_safe <- as.character(paste0("GO_", anno$go_posCellCycle$go_id_number))
anno$go_posCellCycle$go_term_name <- mapIds(GO.db, keys=anno$go_posCellCycle$go_id, column="TERM", keytype="GOID")
anno$go_posCellCycle$gene_ensembl <- mapIds(anno.hs, keys=anno$go_posCellCycle$uniprot, keytype="UNIPROT", column="ENSEMBL")
anno$go_posCellCycle$gene_symbol <- mapIds(anno.hs, keys=anno$go_posCellCycle$uniprot, keytype="UNIPROT", column="SYMBOL")
anno$go_posCellCycle$gene_name <- mapIds(anno.hs, keys=anno$go_posCellCycle$uniprot, keytype="UNIPROT", column="GENENAME")
anno$go_posCellCycle$go_category <- "positive_regulation_cell_cycle"
anno$go_posCellCycle <- anno$go_posCellCycle %>% filter(!is.na(gene_ensembl))
### GO:BIOLOGICAL PROCESS CELL CYCLE REGULATION (NEGATIVE)
anno$go_negCellCycle <- read_delim("input_raw/annotations/amiGO_negRegOfCellCycle.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   col_names = FALSE, trim_ws = TRUE)
colnames(anno$go_negCellCycle) <- c("uniprot", "protein_name", "go_id")
anno$go_negCellCycle <- distinct(anno$go_negCellCycle) %>% 
  tidyr::separate(uniprot, into = c(NA, "uniprot"), sep = ":") %>% 
  tidyr::separate(go_id, into = c(NA, "go_id_number"), sep = ":", remove=FALSE) %>% 
  group_by(go_id) %>% filter(n() >= 5)
anno$go_negCellCycle$go_id_safe <- as.character(paste0("GO_", anno$go_negCellCycle$go_id_number))
anno$go_negCellCycle$go_term_name <- mapIds(GO.db, keys=anno$go_negCellCycle$go_id, column="TERM", keytype="GOID")
anno$go_negCellCycle$gene_ensembl <- mapIds(anno.hs, keys=anno$go_negCellCycle$uniprot, keytype="UNIPROT", column="ENSEMBL")
anno$go_negCellCycle$gene_symbol <- mapIds(anno.hs, keys=anno$go_negCellCycle$uniprot, keytype="UNIPROT", column="SYMBOL")
anno$go_negCellCycle$gene_name <- mapIds(anno.hs, keys=anno$go_negCellCycle$uniprot, keytype="UNIPROT", column="GENENAME")
anno$go_negCellCycle$go_category <- "negative_regulation_cell_cycle"
anno$go_negCellCycle <- anno$go_negCellCycle %>% filter(!is.na(gene_ensembl)) %>% filter(go_id_safe != "GO_0090267")
### create lookup table key for go terms 
anno$go_key <- distinct(bind_rows(anno$go_gpcr[, c(3:5,10,6)], anno$go_posCellCycle[, c(3:5,10,6)], anno$go_negCellCycle[, c(3:5,10,6)]))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("G protein-coupled receptor", "GPCR", go_term_name))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("receptor signaling pathway", "GPCR signaling", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("signaling pathway", "signaling", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("G protein-coupled ", "", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("adenylate cyclase", "AC", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("phospholipase C", "PLC", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("protein kinase C", "PKC", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("positive regulation", "pos. reg.", short_desc))
anno$go_key <- anno$go_key %>% mutate(short_desc = gsub("negative regulation", "neg. reg.", short_desc))
anno$go_key$new_col_names <- paste0(anno$go_key$short_desc, " (", anno$go_key$go_id_safe, ")")

anno$go_all <- bind_rows(anno$go_gpcr, anno$go_posCellCycle, anno$go_negCellCycle)


################## TRANSCRIPTION FACTOR TARGETS 
anno$ChEA_tf_targets <- read_delim("input_raw/annotations/ChEA_harmonizome_gene_attribute_edges.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
anno$ChEA_tf_targets <- anno$ChEA_tf_targets[2:nrow(anno$ChEA_tf_targets), c(1,3,4,6)]
colnames(anno$ChEA_tf_targets) <- c("gene_symbol", "gene_entrez", "tf_symbol", "tf_entrez")
anno$ChEA_tf_targets$gene_ensembl <- mapIds(anno.hs, keys=anno$ChEA_tf_targets$gene_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$ChEA_tf_targets$gene_name <- mapIds(anno.hs, keys=anno$ChEA_tf_targets$gene_entrez, keytype="ENTREZID", column="GENENAME")
anno$ChEA_tf_targets$tf_ensembl <- mapIds(anno.hs, keys=anno$ChEA_tf_targets$tf_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$ChEA_tf_targets$tf_name <- mapIds(anno.hs, keys=anno$ChEA_tf_targets$tf_entrez, keytype="ENTREZID", column="GENENAME")
anno$ChEA_tf_targets$source <- "ChEA_harmonizome"

       

## encode transcription factor targets
anno$encode_tf_targets <- read_delim("input_raw/annotations/ENCODE_harmonizome_gene_attribute_edges.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
anno$encode_tf_targets <- anno$encode_tf_targets[2:nrow(anno$encode_tf_targets), c(1,3,4,6)]
colnames(anno$encode_tf_targets) <- c("gene_symbol", "gene_entrez", "tf_symbol", "tf_entrez")
anno$encode_tf_targets$gene_ensembl <- mapIds(anno.hs, keys=anno$encode_tf_targets$gene_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$encode_tf_targets$gene_name <- mapIds(anno.hs, keys=anno$encode_tf_targets$gene_entrez, keytype="ENTREZID", column="GENENAME")
anno$encode_tf_targets$tf_ensembl <- mapIds(anno.hs, keys=anno$encode_tf_targets$tf_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$encode_tf_targets$tf_name <- mapIds(anno.hs, keys=anno$encode_tf_targets$tf_entrez, keytype="ENTREZID", column="GENENAME")
anno$encode_tf_targets$source <- "encode_harmonizome"


## jaspar transcription factor targets
anno$jaspar_tf_targets <- read_delim("input_raw/annotations/JASPAR_harmonizome_gene_attribute_edges.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
anno$jaspar_tf_targets <- anno$jaspar_tf_targets[2:nrow(anno$jaspar_tf_targets), c(1,3,4,6)]
colnames(anno$jaspar_tf_targets) <- c("gene_symbol", "gene_entrez", "tf_symbol", "tf_entrez")
anno$jaspar_tf_targets$gene_ensembl <- mapIds(anno.hs, keys=anno$jaspar_tf_targets$gene_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$jaspar_tf_targets$gene_name <- mapIds(anno.hs, keys=anno$jaspar_tf_targets$gene_entrez, keytype="ENTREZID", column="GENENAME")
anno$jaspar_tf_targets$tf_ensembl <- mapIds(anno.hs, keys=anno$jaspar_tf_targets$tf_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$jaspar_tf_targets$tf_name <- mapIds(anno.hs, keys=anno$jaspar_tf_targets$tf_entrez, keytype="ENTREZID", column="GENENAME")
anno$jaspar_tf_targets$source <- "jaspar_harmonizome"


## transfac transcription factor targets
anno$transfac_tf_targets <- read_delim("input_raw/annotations/TRANSFAC_harmonizome_gene_attribute_edges.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
anno$transfac_tf_targets <- anno$transfac_tf_targets[2:nrow(anno$transfac_tf_targets), c(1,3,4,6)]
colnames(anno$transfac_tf_targets) <- c("gene_symbol", "gene_entrez", "tf_symbol", "tf_entrez")
anno$transfac_tf_targets$gene_ensembl <- mapIds(anno.hs, keys=anno$transfac_tf_targets$gene_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$transfac_tf_targets$gene_name <- mapIds(anno.hs, keys=anno$transfac_tf_targets$gene_entrez, keytype="ENTREZID", column="GENENAME")
anno$transfac_tf_targets$tf_ensembl <- mapIds(anno.hs, keys=anno$transfac_tf_targets$tf_entrez, keytype="ENTREZID", column="ENSEMBL")
anno$transfac_tf_targets$tf_name <- mapIds(anno.hs, keys=anno$transfac_tf_targets$tf_entrez, keytype="ENTREZID", column="GENENAME")
anno$transfac_tf_targets$source <- "transfac_harmonizome"



anno$union_tf_targets <- bind_rows(anno$ChEA_tf_targets, anno$encode_tf_targets, anno$jaspar_tf_targets, anno$transfac_tf_targets)
anno$union_tf_lists <- distinct(anno$union_tf_targets[, c("tf_symbol", "tf_name", "tf_ensembl", "source")])



#### tf's of possible interest in each dataset: 
#ChEA: ATF3, CREB1, CREM, JUN, KDM5A, KDM5B, KDM6A, MITF, NR4A2
#encode: ATF1, ATF2, ATF3, CREB1, CREBBP, EP300, FOS, FOSL1, FOSL2, JUN, JUND, KDM1A, KDM4A, KDM5A, KDM5B, PPARGC1A
#jaspar: CREB1, FOS, JUN, JUND, NR4A2
#transfac: ATF1, ATF2, ATF3, ATF4, ATF6, CREB1, JUN, NF1

######################### make consensus lists of transcription factor targets/extract data for tfs only present in 1 dataset
## creb1 consensus (target must be present in 3 of 4 databases)
anno$creb1_targets <- list()
anno$creb1_targets$chea <- anno$ChEA_tf_targets %>% filter(tf_symbol == 'CREB1')
anno$creb1_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'CREB1')
anno$creb1_targets$jaspar <- anno$jaspar_tf_targets %>% filter(tf_symbol == 'CREB1')
anno$creb1_targets$transfac <- anno$transfac_tf_targets %>% filter(tf_symbol == 'CREB1')
anno$creb1_targets$lists <- list()
anno$creb1_targets$lists$chea <- anno$creb1_targets$chea$gene_symbol
anno$creb1_targets$lists$encode <- anno$creb1_targets$encode$gene_symbol
anno$creb1_targets$lists$jaspar <- anno$creb1_targets$jaspar$gene_symbol
anno$creb1_targets$lists$transfac <- anno$creb1_targets$transfac$gene_symbol
anno$creb1_targets$mat <- ComplexHeatmap::make_comb_mat(anno$creb1_targets$lists)
anno$creb1_targets$list <- unique(c(ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "1111"), 
                                            ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "1110"), 
                                            ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "1101"),  
                                            ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "1101"), 
                                            ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "0111"),
                                            ComplexHeatmap::extract_comb(anno$creb1_targets$mat, "1100")))
anno$creb1_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% 
                  filter(tf_symbol == 'CREB1', gene_symbol %in% anno$creb1_targets$list) %>% distinct()


## jun consensus (target must be present in 2 of 4 databases)
anno$jun_targets <- list()
anno$jun_targets$chea <- anno$ChEA_tf_targets %>% filter(tf_symbol == 'JUN')
anno$jun_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'JUN')
anno$jun_targets$jaspar <- anno$jaspar_tf_targets %>% filter(tf_symbol == 'JUN')
anno$jun_targets$transfac <- anno$transfac_tf_targets %>% filter(tf_symbol == 'JUN')
anno$jun_targets$lists <- list()
anno$jun_targets$lists$chea <- anno$jun_targets$chea$gene_symbol
anno$jun_targets$lists$encode <- anno$jun_targets$encode$gene_symbol
anno$jun_targets$lists$jaspar <- anno$jun_targets$jaspar$gene_symbol
anno$jun_targets$lists$transfac <- anno$jun_targets$transfac$gene_symbol
anno$jun_targets$mat <- ComplexHeatmap::make_comb_mat(anno$jun_targets$lists)
anno$jun_targets$list <- unique(c(ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1111"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1110"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1101"),  
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1101"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "0111"),
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1100"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1010"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "1001"), 
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "0110"),
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "0101"),
                                              ComplexHeatmap::extract_comb(anno$jun_targets$mat, "0011")))
anno$jun_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>%
                filter(tf_symbol == 'JUN', gene_symbol %in% anno$jun_targets$list) %>% distinct()


## jund consensus (target must be present in 2 of 2 databases)
anno$jund_targets <- list()
anno$jund_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'JUND')
anno$jund_targets$jaspar <- anno$jaspar_tf_targets %>% filter(tf_symbol == 'JUND')
anno$jund_targets$lists$encode <- anno$jund_targets$encode$gene_symbol
anno$jund_targets$lists$jaspar <- anno$jund_targets$jaspar$gene_symbol
anno$jund_targets$mat <- ComplexHeatmap::make_comb_mat(anno$jund_targets$lists)
anno$jund_targets$list <- unique(c(ComplexHeatmap::extract_comb(anno$jund_targets$mat, "11")))
anno$jund_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% 
                  filter(tf_symbol == 'JUND', gene_symbol %in% anno$jund_targets$list) %>% distinct()


## nr4a2 consensus (target must be present in 1 of 2 databases)
anno$nr4a2_targets <- list()
anno$nr4a2_targets$chea <- anno$ChEA_tf_targets %>% filter(tf_symbol == 'NR4A2')
anno$nr4a2_targets$jaspar <- anno$jaspar_tf_targets %>% filter(tf_symbol == 'NR4A2')
anno$nr4a2_targets$list <- unique(c(anno$nr4a2_targets$chea$gene_symbol,
                                        anno$nr4a2_targets$jaspar$gene_symbol))
anno$nr4a2_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% 
                  filter(tf_symbol == 'NR4A2', gene_symbol %in% anno$nr4a2_targets$list) %>% distinct()


## fos consensus (target must be present in 2 of 2 databases)
anno$fos_targets <- list()
anno$fos_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'FOS')
anno$fos_targets$jaspar <- anno$jaspar_tf_targets %>% filter(tf_symbol == 'FOS')
anno$fos_targets$lists$encode <- anno$fos_targets$encode$gene_symbol
anno$fos_targets$lists$jaspar <- anno$fos_targets$jaspar$gene_symbol
anno$fos_targets$mat <- ComplexHeatmap::make_comb_mat(anno$fos_targets$lists)
anno$fos_targets$list <- unique(c(ComplexHeatmap::extract_comb(anno$fos_targets$mat, "11")))
anno$fos_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% 
                filter(tf_symbol == 'FOS', gene_symbol %in% anno$fos_targets$list) %>% distinct()



## get df for transcription factors only present in one dataset
anno$crem_targets <- list()
    anno$crem_targets$chea <- anno$ChEA_tf_targets %>% filter(tf_symbol == 'CREM')
    anno$crem_targets$list <- anno$crem_targets$chea$gene_symbol
    anno$crem_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'CREM')

    
anno$mitf_targets <- list()
    anno$mitf_targets$chea <- anno$ChEA_tf_targets %>% filter(tf_symbol == 'MITF')
    anno$mitf_targets$list <- anno$mitf_targets$chea$gene_symbol
    anno$mitf_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'MITF')

    
anno$crebbp_targets <- list()
    anno$crebbp_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'CREBBP')
    anno$crebbp_targets$list <- anno$crebbp_targets$encode$gene_symbol
    anno$crebbp_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'CREBBP')

    
anno$ep300_targets <- list()
    anno$ep300_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'EP300')
    anno$ep300_targets$list <- anno$ep300_targets$encode$gene_symbol
    anno$ep300_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'EP300')

    
anno$fosl1_targets <- list()
    anno$fosl1_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'FOSL1')
    anno$fosl1_targets$list <- anno$fosl1_targets$encode$gene_symbol
    anno$fosl1_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'FOSL1')    

    
anno$fosl2_targets <- list()
    anno$fosl2_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'FOSL2')
    anno$fosl2_targets$list <- anno$fosl2_targets$encode$gene_symbol
    anno$fosl2_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'FOSL2')        

    
anno$ppargc1a_targets <- list()
    anno$ppargc1a_targets$encode <- anno$encode_tf_targets %>% filter(tf_symbol == 'PPARGC1A')
    anno$ppargc1a_targets$list <- anno$ppargc1a_targets$encode$gene_symbol
    anno$ppargc1a_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'PPARGC1A')   

    
anno$nf1_targets <- list()
    anno$nf1_targets$transfac <- anno$transfac_tf_targets %>% filter(tf_symbol == 'NF1')
    anno$nf1_targets$list <- anno$nf1_targets$transfac$gene_symbol
    anno$nf1_targets$df <- anno$union_tf_targets[ , 1:c(ncol(anno$union_tf_targets)-1)] %>% filter(tf_symbol == 'NF1')   

    
############################
    dir.create(file.path("tidy_input/annotations"))
    dir.create(file.path("tidy_input/annotations/reference_tables"))
    dir.create(file.path("tidy_input/annotations/individual_tfs"))
    
    write.csv(anno$tsoi_subtype_signatures, file = "tidy_input/annotations/tsoi_subtype_signatures.csv", quote = TRUE)
    write.csv(anno$genes_of_interest, file = "tidy_input/annotations/known_genes_of_interest.csv", quote = TRUE)
    write.csv(anno$bioplanet, file = "tidy_input/annotations/reference_tables/bioplanet.csv", quote = TRUE)
    write.csv(anno$bioplanet_gpcr, file = "tidy_input/annotations/bioplanet_gpcr.csv", quote = TRUE)
    write.csv(anno$bioplanet_key, file = "tidy_input/annotations/bioplanet_key.csv", quote = TRUE)
    
    ### write go term related annotation tables to csv files
    write.csv(anno$go_gpcr, file = "tidy_input/annotations/reference_tables/go_gpcr.csv", quote=TRUE)
    write.csv(anno$go_posCellCycle, file = "tidy_input/annotations/reference_tables/go_posCellCycle.csv", quote=TRUE)
    write.csv(anno$go_negCellCycle, file = "tidy_input/annotations/reference_tables/go_negCellCycle.csv", quote=TRUE)
    write.csv(anno$go_key, file = "tidy_input/annotations/go_key.csv", quote=TRUE)
    write.csv(anno$go_all, file = "tidy_input/annotations/go_all.csv", quote=TRUE)
    
    #write unfiltered tf_target dataframes to csv files
    write.csv(anno$ChEA_tf_targets, file = "tidy_input/annotations/reference_tables/ChEA_tf_targets.csv", quote=TRUE)
    write.csv(anno$encode_tf_targets, file = "tidy_input/annotations/reference_tables/encode_tf_targets.csv", quote=TRUE)
    write.csv(anno$jaspar_tf_targets, file = "tidy_input/annotations/reference_tables/jaspar_tf_targets.csv", quote=TRUE)
    write.csv(anno$transfac_tf_targets, file = "tidy_input/annotations/reference_tables/transfac_tf_targets.csv", quote=TRUE)
    write.csv(anno$union_tf_targets, file = "tidy_input/annotations/reference_tables/union_harmonizome_tf_targets.csv", quote=TRUE)
    write.csv(anno$union_tf_lists, file = "tidy_input/annotations/reference_tables/union_harmonizome_tf_list.csv", quote=TRUE)
    
    # write specific tf targets to csv
    write.csv(anno$creb1_targets$df, file = "tidy_input/annotations/individual_tfs/creb1_targets.csv", quote=TRUE)
    write.csv(anno$jun_targets$df, file = "tidy_input/annotations/individual_tfs/jun_targets.csv", quote=TRUE)
    write.csv(anno$jund_targets$df, file = "tidy_input/annotations/individual_tfs/jund_targets.csv", quote=TRUE)
    write.csv(anno$nr4a2_targets$df, file = "tidy_input/annotations/individual_tfs/nr4a2_targets.csv", quote=TRUE)
    write.csv(anno$fos_targets$df, file = "tidy_input/annotations/individual_tfs/fos_targets.csv", quote=TRUE)
    write.csv(anno$crem_targets$df, file = "tidy_input/annotations/individual_tfs/crem_targets.csv", quote=TRUE)
    write.csv(anno$mitf_targets$df, file = "tidy_input/annotations/individual_tfs/mitf_targets.csv", quote=TRUE)
    write.csv(anno$crebbp_targets$df, file = "tidy_input/annotations/individual_tfs/crebbp_targets.csv", quote=TRUE)
    write.csv(anno$ep300_targets$df, file = "tidy_input/annotations/individual_tfs/ep300_targets.csv", quote=TRUE)
    write.csv(anno$fosl1_targets$df, file = "tidy_input/annotations/individual_tfs/fosl1_targets.csv", quote=TRUE)
    write.csv(anno$ppargc1a_targets$df, file = "tidy_input/annotations/individual_tfs/ppargc1a_targets.csv", quote=TRUE)
    write.csv(anno$ppargc1a_targets$df, file = "tidy_input/annotations/individual_tfs/ppargc1a_targets.csv", quote=TRUE)
    write.csv(anno$nf1_targets$df, file = "tidy_input/annotations/individual_tfs/nf1_targets.csv", quote=TRUE)
    
################    
    
    rm(list = ls())
    gc()    
    
