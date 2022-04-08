##### LOAD PACKAGES #####
libraries <- c("DESeq2", "pheatmap", "vsn", "hexbin", "cowplot", "apeglm", "tximport", "gprofiler2",
               "viridis", "PoiClaClu", "genefilter", "vidger", "reshape2", "ggplot2", "stringr", 
               "gplots", "sva", "tidyverse", "biomaRt", "RColorBrewer", "gprofiler2", "DEGreport",
               "tximport", "enrichplot", "DOSE", "org.Hs.eg.db", "plotly", "clusterProfiler", "gprofiler2",
               "ape", "motmot", "geiger", "phytools", "OUwie", "ggtree", "casper", "taxizedb", "simplifyEnrichment",
               "ggsci", "Biostrings", "rtracklayer", "GenomicFeatures", "rstatix", "ggpubr", "BSgenome", "zFPKM", "ComplexHeatmap")
lapply(libraries, library, character.only = TRUE)

##### FUNCTION #######
#Function for extracting gene IDs where expression values in oocytes are >= 2
expr.extr <- function(variables){
  mat_ids <- c() #Initialize empty result vector
  
  for(i in 1:length(variables)){
    temp <- get(variables[i])$abundance
    colnames(temp) <- str_remove_all(colnames(temp), "_[0-9]+$")
    
    #Using zFPKM to subset for expressed genes in egg stages
    temp_expr <- rownames(temp[rowMeans(as.matrix(temp)[, grepl("stage1", colnames(temp)), drop = F]) >= 2, ])
    
    #Subsetting TPM table for expressed genes
    temp <- subset(temp, rownames(temp) %in% temp_expr)
    temp <- rownames(temp)
    
    mat_ids <- c(mat_ids, temp)
  }
  
  return(mat_ids)
}

source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### LOAD SALMON COUNTS #####
files_ext <- list.files(path = "/home/ferenkagan/Documents/H_dujardini",  pattern = "*salmon")
files <- file.path("/home/ferenkagan/Documents/H_dujardini/", files_ext, "quant.sf")
names(files) <- c(paste(rep("stage1", 4), 1:4, sep = "_"),
                  paste(rep("stage2", 3), 1:3, sep = "_"),
                  paste(rep("stage3", 3), 1:3, sep = "_"),
                  paste(rep("stage4", 3), 1:3, sep = "_"),
                  paste(rep("stage5", 3), 1:3, sep = "_"),
                  paste(rep("stage6", 3), 1:3, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("/home/ferenkagan/Documents/H_dujardini//tx2gene.tsv", sep = "\t")
txi_hduj <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

##### CREATE ANNOTATION TABLE #####
#Load UniProt annotation crossreferences
crosref_names <- read.csv("~/Documents/H_dujardini/hduj_annot.csv", header = F) #name crossreference
crosref <- read.csv("~/Documents/H_dujardini/hduj_annot_id.csv", header = F) #ID crosreference
crosref$V2 <- str_remove_all(str_remove_all(str_remove_all(crosref$V2, "sp|"), "^\\|"), "\\|.+$")
uniprot_entrez <- read.table("~/Documents/uniprot-entrez.tsv", header = F, sep = "\t") #uniprot ID to entrez 
crosref$V3 <- uniprot_entrez$V2[match(crosref$V2, uniprot_entrez$V1)]
crosref$V4 <- crosref_names$V2[match(crosref$V1, crosref_names$V1)]

#Results from reciprocal best hit blast against human proteome
ensembl <- read.csv("~/Documents/H_dujardini/rbh_hduj_stripped.tbl", header = F)
ensembl$V1 <- str_remove_all(ensembl$V1, "t[0-9]+$")
ensembl$V2 <- str_remove_all(ensembl$V2, "\\.[0-9]+$")

#Getting Entrez IDs for Enesmbl IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#View(listFilters(mart))
entrez_tbl <- getBM(attributes = c("ensembl_peptide_id",  "entrezgene_id", "external_gene_name"),  
                    values = ensembl$V2, bmHeader = T, mart = mart)
ensembl$V3 <- entrez_tbl$`NCBI gene (formerly Entrezgene) ID`[match(ensembl$V2, entrez_tbl$`Protein stable ID`)]

#Filling some missing entrez IDs using gprofiler
gConv <- gconvert(ensembl$V2, filter_na = T,
                  target = "ENTREZGENE_ACC")[, c(2, 4)]
gConv <- as.data.frame(apply(gConv, 2, function(x){str_replace_all(x, "nan", NA_character_)}))
ensembl <- unique(merge(ensembl, gConv, by.x = "V2", by.y = "input", all = T))
ensembl <-tibble(ensembl) %>%                                            
  mutate(EntrezID = coalesce(as.integer(V3), as.integer(target)))

#Assembling final crosreference table
crosref$V5 <- ensembl$V3[match(crosref$V1, ensembl$V1)]
crosref <- tibble(crosref) %>%                                            
  mutate(EntrezIDComb = coalesce(as.integer(V5), as.integer(V3)))

#Reading GO annotations
go_annot <- read_tsv("~/Documents/Gene_expr_evol/GO/Annotations/hduj_go.tsv") %>% filter(ARGOT_PPV >= 0.5)
go_annot$qpid <- str_remove_all(go_annot$qpid, "t[0-9]+$")
go_annot$goid <- paste0("GO:", go_annot$goid)
go_annot <- go_annot %>% dplyr::select(qpid, goid, desc)

##### BUILD DESEQ OBJECT ######
coldata_hduj <- data.frame(row.names=colnames(txi_hduj$abundance), 
                           Stages_hduj = factor(str_remove_all(colnames(txi_hduj$abundance), "_[0-9]+")))
dds_hduj <- DESeqDataSetFromTximport(txi_hduj, colData = coldata_hduj, design = ~Stages_hduj)

##### LOW COUNT FILTER #####
dds_hduj <- estimateSizeFactors(dds_hduj)
keep <- rowSums(counts(dds_hduj, normalized = T) >= 10) >= 2
dds_hduj <- dds_hduj[keep,]

##### VARIANCE STABILIZING #####

dds_rld <- rlogTransformation(dds_hduj)
dds_vst <- vst(dds_hduj)

(p <- plotPCA(dds_vst, intgroup=c("Stages_hduj")))
meanSdPlot(assay(dds_vst))

##### BATCH CORRECTION #####

dds_hduj <- DESeq(dds_hduj)

dat <- counts(dds_hduj, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~Stages_hduj, colData(dds_hduj))
mod0 <- model.matrix(~1, colData(dds_hduj))
svseq <- svaseq(dat, mod, mod0)

#5 surrogate variable were found, including this in the modeling
ddssva_hduj <- dds_hduj
SV1 <- svseq$sv[, 1]
SV2 <- svseq$sv[, 2]
SV3 <- svseq$sv[, 3]
SV4 <- svseq$sv[, 4]
SV5 <- svseq$sv[, 5]
colData(ddssva_hduj) <- cbind(colData(ddssva_hduj), SV1)
colData(ddssva_hduj) <- cbind(colData(ddssva_hduj), SV2)
colData(ddssva_hduj) <- cbind(colData(ddssva_hduj), SV3)
colData(ddssva_hduj) <- cbind(colData(ddssva_hduj), SV4)
colData(ddssva_hduj) <- cbind(colData(ddssva_hduj), SV5)
design(ddssva_hduj) <- ~Stages_hduj + SV1 + SV2 + SV3 + SV4 +SV5
ddssva_hduj <- DESeq(ddssva_hduj)

#Excluding gene without convergence
ddssva_hduj <- ddssva_hduj[which(mcols(ddssva_hduj)$betaConv),]

###### EXTRACT RESULTS #####
#Contrasts for MZT and extracting results
contrast_hduj_1 <- c("Stages_hduj", "stage1", "stage2")
res_1 <- results(ddssva_hduj, contrast_hduj_1, alpha = 0.05)
summary(res_1)

#Extracting gene IDs where oocyte stages have expression
maternal <- expr.extr(ls(pattern = "txi_....$"))

#Extracting gene IDs which belong to degraded subset
mzt_deg <- rownames(subset(res_1, res_1$log2FoldChange <= -2 & res_1$padj <= 0.05))

#What proportion of the transcriptome is expressed?
length(maternal) / dim(txi_hduj$abundance)[1]

#Mean expressions of maternal genes
temp <- collapseReplicate(subset(txi_hduj$abundance, rownames(txi_hduj$abundance) %in% maternal)[, 1:4])
head(temp[order(temp, decreasing = T)], 10)

#What are the annotations for top genes?
crosref[crosref$V1 %in% names(head(temp[order(temp, decreasing = T)], 70)), ]
subset(pfam, pfam$query_name %in% names(head(temp[order(temp, decreasing = T)], 20))) #pfam variable is read in later

##### UTR ANALYSIS #####
#This section was used for analyzing UTR motif enrichment and motif occurence results in combination with
#expression data

#Read in enrichment results from the two separate databases
crosref <- read.csv("~/Documents/H_dujardini/hduj_annot.csv", header = F)
att_crosref <- read.table("~/Documents/Gene_expr_evol/Gene_lengths/ATtRACT_db.txt", sep = "\t", header = T)
sea_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_dge_ray2013.out",
                          sep = "\t", header = T)
sea_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_dge_attract.out",
                          sep = "\t", header = T)
sea_dge_att$ALT_ID <- att_crosref$Gene_name[match(sea_dge_att$ID, att_crosref$Matrix_id)]
sea_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_mat_ray2013.out",
                          sep = "\t", header = T)
sea_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_mat_attract.out",
                          sep = "\t", header = T)
sea_mat_att$ALT_ID <- att_crosref$Gene_name[match(sea_mat_att$ID, att_crosref$Matrix_id)]

#Combine the two database based results
sea_dge <- rbind(sea_dge_att, sea_dge_ray)
sea_mat <- rbind(sea_mat_att, sea_mat_ray)

#Read in motif occurence results for the two databases for the degraded subset of genes
fimo_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_dge_ray2013.fimo",
                           sep = "\t", header = T)
fimo_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_dge_attract.fimo",
                           sep = "\t", header = T)

#Combine the two results with the results of simple enrichment results
fimo_dge <- rbind(fimo_dge_att, fimo_dge_ray)
fimo_dge <- subset(fimo_dge, fimo_dge$motif_id %in% sea_dge$ID)

#Format result table
fimo_dge$motif_alt_id <- sea_dge$ALT_ID[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sequence_alt_name <- crosref$V2[match(fimo_dge$sequence_name, crosref$V1)]
fimo_dge$enr_val <- sea_dge$ENR_RATIO[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sea.q.value <- sea_dge$QVALUE[match(fimo_dge$motif_id, sea_dge$ID)]

#Read in motif occurence results for the two databases for the non-degraded subset of genes
fimo_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_mat_ray2013.fimo",
                           sep = "\t", header = T)
fimo_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Hduj_mat_attract.fimo",
                           sep = "\t", header = T)

#Combine the two results with the results of simple enrichment results
fimo_mat <- rbind(fimo_mat_att, fimo_mat_ray)
fimo_mat <- subset(fimo_mat, fimo_mat$motif_id %in% sea_mat$ID)

#Format result table
fimo_mat$motif_alt_id <- sea_mat$ALT_ID[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sequence_alt_name <- crosref$V2[match(fimo_mat$sequence_name, crosref$V1)]
fimo_mat$enr_val <- sea_mat$ENR_RATIO[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sea.q.value <- sea_mat$QVALUE[match(fimo_mat$motif_id, sea_mat$ID)]

#Inspect motif occurence for each transcript
fimo_dge <- tibble(fimo_dge) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)
fimo_mat <- tibble(fimo_mat) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)

temp <- t(scale(t(collapseReplicate(assay(dds_rld)))))
rownames(temp) <- rownames(assay(dds_rld))

#Subsetting expression data for plotting target gene expressions
temp <- subset(temp, rownames(temp) %in% subset(fimo_mat, fimo_mat$motif_alt_id %in% c("AGO1", "AGO2"))$sequence_name)
colnames(temp) <- c("oocyte", "segmentation", "stage 17", "stage18", "hatching", "juvenile")

#svg("~/Desktop/Figures/fimo_heatmap_Hduj3.svg")
#png("~/Desktop/Figures/fimo_heatmap_Hduj3.png",
#    width=750, height=1000)
print(ComplexHeatmap::Heatmap(matrix = temp,
                              column_title = "AGO1/2 motif containing genes",
                              show_row_names = F,
                              show_row_dend = T,
                              row_dend_reorder = F,
                              cluster_columns = F,
                              col = rev(brewer.pal(9,"RdBu")),
                              column_names_rot = 45, column_names_gp = gpar(fontface = "italic", fontsize = 10),
                              name = "Gene Z-score"))
#dev.off()

#GO enrichment of motif target genes
#Preparing data
ids <- subset(fimo_dge, fimo_dge$q.value <= 0.05)
variables <- unique(ids$motif_alt_id)
utr_list <- vector(mode = "list", length = length(variables))
resrbind <- as.data.frame(res_1)

#Looping through RBP target data and performing enrichment
for(k in 1: length(variables)){
  temp_name <- variables[k] 
  temp_ids <- subset(ids, ids$motif_alt_id == variables[k]) #Enrichment table subsetting for target motif
  
  temp_res <- resrbind[order(abs(resrbind$log2FoldChange), decreasing = T),]
  temp_res <- subset(temp_res, rownames(temp_res) %in% unique(temp_ids$sequence_name)) #Subsetting expression data according to motif targets
  
  #Perform enrichment analysis
  temp_ego<- enricher(rownames(temp_res), TERM2GENE = go_annot[,c(2,1)], TERM2NAME = go_annot[, c(2,3)],
                      minGSSize = 5)
  temp_ego@result <- temp_ego@result[temp_ego@result$qvalue <= 0.05, ]
  
  #If no enrichment found to avoid error message just skip
  if(is.null(temp_ego@result)){next}
  
  utr_list[[k]] <- temp_ego@result
  names(utr_list)[[k]] <- temp_name
  
}


##### RBP #####
#This section was used to analyze RNA binding domain (RBD) abundances within maternal transcriptome

#Read in data
pfam <- rhmmer::read_tblout("/home/ferenkagan/Documents/H_dujardini/hduj_v_pfam.tblout")
#Subset for RNA binding the whole table, this is used as background universe for enrichment
pfam$query_name <- str_remove_all(pfam$query_name, "t[0-9]+$")
pfam_rbp <- pfam[str_detect(pfam$description, "RNA"), ]
pfam_rbp <- pfam_rbp[str_detect(pfam_rbp$description, "binding"), ]

#Subset the HMMER annotation with hduj results
maternal_domains <- subset(pfam, pfam$query_name %in% maternal)
mzt_dge_domains <- subset(pfam, pfam$query_name %in% mzt_deg)
mzt_zga_domains <- subset(pfam, pfam$query_name %in% mzt_zga)

#Subsetting with grep with searchwords of RNA and binding
maternal_domains_rbp <- maternal_domains[str_detect(maternal_domains$description, "RNA"), ]
mzt_dge_domains_rbp <- mzt_dge_domains[str_detect(mzt_dge_domains$description, "RNA"), ]
mzt_zga_domains_rbp <- mzt_zga_domains[str_detect(mzt_zga_domains$description, "RNA"), ]
reference_domains_rbp <- pfam[str_detect(pfam$description, "RNA"), ]
 
maternal_domains_rbp <- maternal_domains_rbp[str_detect(maternal_domains_rbp$description, "binding"), ]
mzt_dge_domains_rbp <- mzt_dge_domains_rbp[str_detect(mzt_dge_domains_rbp$description, "binding"), ]
mzt_zga_domains_rbp <- mzt_zga_domains_rbp[str_detect(mzt_zga_domains_rbp$description, "binding"), ]
reference_domains_rbp <- pfam[str_detect(pfam$description, "binding"), ]

#Proportions of RBP at different groupings
length(unique(maternal_domains_rbp$query_name)) / length(unique(maternal_domains$query_name))
length(unique(mzt_dge_domains_rbp$query_name)) / length(unique(mzt_dge_domains$query_name))
length(unique(mzt_zga_domains_rbp$query_name)) / length(unique(mzt_zga_domains$query_name))
length(unique(pfam_rbp$query_name)) / length(unique(pfam$query_name))

#For DGE genes hypergeometric test
phyper(length(unique(mzt_dge_domains_rbp$query_name)), 
       length(unique(mzt_dge_domains$query_name)), 
       length(unique(pfam$query_name)) - length(unique(mzt_dge_domains$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)

#For maternal genes hypergeometric test
phyper(length(unique(maternal_domains_rbp$query_name)), 
       length(unique(maternal_domains$query_name)), 
       length(unique(pfam$query_name)) - length(unique(maternal_domains$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)

#Control sequences from later development
phyper(length(unique(reference_domains_rbp$query_name)), 
       length(unique(contr_1_domains$query_name)), 
       length(unique(pfam$query_name)) - length(unique(contr_1_domains$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)
phyper(length(unique(contr_2_domains_rbp$query_name)), 
       length(unique(contr_2_domains$query_name)), 
       length(unique(pfam$query_name)) - length(unique(contr_2_domains$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)

##### MFUZZ CLUSTERING #####
#This section was used for fuzzy clustering of gene expression data and GO enrichment on resultant clusters
library(Mfuzz)

#Collapsing replicates
cts <- DESeq2::collapseReplicates(ddssva_hduj, groupby = coldata_hduj$Stages_hduj)
cts <- assay(rlog(cts))

#Adding NA to 0 counts for subsequent filtering
cts[cts == 0] <- NA

#Missing values
cts <- filter.NA(ExpressionSet(cts))
cts <- fill.NA(cts)

#Filtering
cts <- filter.std(cts, min.std = 0)

#Standardise
cts <- standardise(cts)

#Estimating clustering parameters
m <- mestimate(cts)
c <- Dmin(cts, m = m, crange = seq(6,40,2), repeats = 5)

#Fuzzy clustering, centers determined by Dmin plot inspection
cts_clst <- mfuzz(cts, centers = 12, m = m)

stages <- c("oocyte", "segmentation", "stage 17", "stage18", "hatchling", "juvenile")

#Export plot for each cluster
for(i in 1:length(unique(cts_clst$cluster))){
  svg(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_hduj", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_hduj", ".png"))),
  #  width=550, height=550)
  print(mfuzz.plot2(cts, cts_clst, min.mem = 0.5, colo = 'fancy', time.labels = stages,
              centre = T, single = i, x11 = F))
  print(mtext(paste0("Cluster size: ", sum(cts_clst$cluster %in% i))))
  dev.off()
}

#Get cluster cores for enrichment analysis of those
cores <- acore(cts, cts_clst)

#Initializing variables
variables <- unique(cts_clst$cluster)
go_list <- vector(mode = "list", length = length(variables))
resrbind <- as.data.frame(res_1)
resrbind <- resrbind[unique(rownames(resrbind)), ]

#Looping through clusters for enrichment
for(k in 1:length(variables)){
  temp_clust <- cores[[k]]$NAME
  temp_name <- paste0("Cluster_", k)
  
  #Perform enrichment analysis
  temp_ego<- enricher(temp_clust, TERM2GENE = go_annot[,c(2,1)], TERM2NAME = go_annot[, c(2,3)])
  temp_ego@result <- temp_ego@result[temp_ego@result$qvalue <= 0.05, ]
  
  #If no enrichment found to avoid error message just skip
  if(is.null(temp_ego@result)){next}
  
  go_list[[k]] <- temp_ego@result
  names(go_list)[[k]] <- temp_name
  
}


##### BUILDING EXPRESSION PATHWAY HEATMAPS ######
#This section was used to query KEGG pathways within Hduj transcriptome
library(KEGGREST)

rld <- collapseReplicate(as.data.frame(assay(dds_rld)))

#These pathways are queried, for further information visit: https://www.genome.jp/kegg/pathway.html
keggid <- c("04010", "04014", "04310", "04330", "04340", "04350", "04390", "04630", "04064",
            "04668", "04068", "04020", "04070", "04072", "04024", "04151", "04514", "04512",
            "04210", "04360")

#Initializing variables
expr_list <- vector(mode = "list", length = length(keggid))
missing <- c()

for(i in 1:length(keggid)){
  #Getting pathway name
  temp_name <- keggFind(database = "pathway", query = keggid[i])
  
  #Getting genes involved in pathway
  temp_id <- keggGet(paste("hsa", keggid[i], sep =  ""))[[1]]$GENE
  temp_id_df <- data.frame(Entrez = temp_id[seq(2,length(temp_id),2)-1],
                           Name = str_remove_all(temp_id[seq(2,length(temp_id),2)], ";.*"))
  
  #What is the missing proportion of genes
  missing[i] <- 1 - (length(crosref$V3[crosref$V3 %in% temp_id_df$Entrez]) / 
                       length(temp_id_df$Entrez))
  names(missing)[i] <- temp_name
  
  #Build expression table
  temp_id <- rownames(rld)
  temp_id <- crosref$V3[match(temp_id, crosref$V1)]
  mat_temp <- as.data.frame(rld[temp_id %in% temp_id_df$Entrez, ])
  
  temp_id_gn <- crosref$V3[match(rownames(mat_temp), crosref$V1)]
  temp_id_gn <- temp_id_df$Name[match(temp_id_gn, temp_id_df$Entrez)]
  rownames(mat_temp) <- paste(rownames(mat_temp), temp_id_gn, sep = "_")
  
  #Assign exprs table to list
  expr_list[[i]] <- mat_temp
  names(expr_list)[i] <- temp_name
}

#Plotting
anno <- data.frame(Stages = colnames(expr_list[[1]]))

for(i in 1:length(expr_list)){
  temp <-  apply(expr_list[[i]], 2, function(x){as.double(x)})
  temp <- t(scale(t(temp)))
  
  rownames(temp) <- rownames(expr_list[[i]])
  colnames(temp) <- c("oocyte", "segmentation", "stage 17", "stage18", "hatching", "juvenile")
  
  ha <- rowAnnotation(foo = anno_mark(at = match(unique(str_remove_all(rownames(temp), "Hduj[0-9]+_")), str_remove_all(rownames(temp), "Hduj[0-9]+_")), 
                                      labels = rownames(temp[match(unique(str_remove_all(rownames(temp), "Hduj[0-9]+_")), str_remove_all(rownames(temp), "Hduj[0-9]+_")), ])))
  
  
  #svg(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("_hduj", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("_hduj", ".png"))),
  #    width=750, height=1000)
   print(ComplexHeatmap::Heatmap(matrix = temp,
                                column_title = names(expr_list)[i],
                                show_row_names = F,
                                right_annotation = ha,
                                show_row_dend = T,
                                row_dend_reorder = F,
                                cluster_columns = F,
                                col = rev(brewer.pal(9,"RdBu")),
                                column_names_rot = 45, column_names_gp = gpar(fontface = "italic", fontsize = 10),
                                name = "Gene Z-score",
                                rect_gp = gpar(col = "white", lwd = 0.5)))
  #dev.off()
  
}

###### PLOT INDIVIDUAL GENES #####
#This section was used for plotting gene expression values across development
library(wesanderson)

#rlog transformed data
rld <- as.data.frame(assay(dds_rld))

#Target genes
targets <- c("Hduj000470", "Hduj005359", "Hduj012802", "Hduj002128", "Hduj003414")

#Loopin through each target gene
for(i in 1:length(targets)){
  #Could be missing from table
  if(!(targets[i] %in% rownames(rld))){next}
  
  #Subsetting for gene and formatting names
  expr <- as.data.frame(t(subset(rld, rownames(rld) %in% targets[i])))
  colnames(expr) <- "expr"
  expr$replicate <- str_remove_all(rownames(expr), "stage[0-9]+_")
  expr$staging <- factor(c(rep("oocyte", 4), rep("segmentation", 3), rep("stage 17", 3), rep("stage 18", 3), rep("hatching", 3), rep("juvenile", 3)),
                         levels = c("oocyte", "segmentation", "stage 17", "stage 18", "hatching", "juvenile"))
  
  #svg(paste0("~/Desktop/Figures/", paste0(crosref$V4[crosref$V1 == targets[i]], paste0("_hduj", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(crosref$V4[crosref$V1 == targets[i]], paste0("_hduj", ".png"))),
  #    width=550, height=550)
  print(ggplot(expr, aes(x = staging, y = expr, color = replicate, group = replicate)) +
          geom_point(size = 6) +
          #annotate("text", x=5, y=10, label= crosref$V4[crosref$V1 == targets[i]], size = 9) +
          theme_bw() +
          ggtitle(crosref$V4[crosref$V1 %in% targets[i]]) +
          scale_color_manual(values=wes_palettes$Rushmore[seq(1, length(wes_palettes), by = 2)]) +
          xlab("") +
          ylab("Regularized log transformed expression") +
          theme(text = element_text(size = 18)))
  #dev.off()
  
}


