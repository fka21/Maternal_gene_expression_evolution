##### LOAD PACKAGES #####
libraries <- c("DESeq2", "pheatmap", "vsn", "hexbin", "cowplot", "apeglm", "tximport", "gprofiler2",
               "viridis", "PoiClaClu", "genefilter", "vidger", "reshape2", "ggplot2", "stringr", 
               "gplots", "sva", "tidyverse", "biomaRt", "RColorBrewer", "gprofiler2", "DEGreport",
               "tximport", "enrichplot", "DOSE", "org.Hs.eg.db", "plotly", "clusterProfiler", "gprofiler2",
               "ape", "motmot", "geiger", "phytools", "OUwie", "ggtree", "casper", "taxizedb", "simplifyEnrichment",
               "ggsci", "Biostrings", "rtracklayer", "GenomicFeatures", "rstatix", "ggpubr", "BSgenome", "zFPKM")
lapply(libraries, library, character.only = TRUE)

##### FUNCTION #######

#Function for extracting gene IDs where expression values in oocytes are >= 2
expr.extr <- function(variables){
  mat_ids <- c() #Initialize empty result vector
  
  for(i in 1:length(variables)){
    temp <- get(variables[i])$abundance
    
    #Using cutoff to subset for expressed genes in egg stages
    temp_expr <- rownames(temp[rowMeans(as.matrix(temp)[, grepl("S01", colnames(temp)), drop = F]) >= 2, ])
    
    #Subsetting TPM table for expressed genes
    temp <- subset(temp, rownames(temp) %in% temp_expr)
    temp <- rownames(temp)
    
    mat_ids <- c(mat_ids, temp)
  }
  
  return(mat_ids)
}

collapseReplicate <- function(df){
  temp <- as.data.frame(t(df))
  
  temp$ID <- str_remove_all(rownames(temp), "_F[0-9]+")
  
  temp2 <- tibble(temp) %>% dplyr::group_by(ID) %>%
    dplyr::summarize(across(1:dim(.)[2]-1, mean, na.rm = T))
  
  temp <- as.data.frame(temp2)[, -1]
  rownames(temp) <- temp2$ID; temp <- t(temp)
  return(temp)
}


##### LOAD SALMON COUNTS #####
files_ext <- list.files(path = "/home/ferenkagan/Documents/T_transversa/quantification",  pattern = "*salmon")
files <- file.path("/home/ferenkagan/Documents/T_transversa/quantification", files_ext, "quant.sf")
names(files) <- paste(str_extract_all(files_ext, "S[0-9]+", simplify = T),
      str_extract_all(files_ext, "F[0-9]", simplify = T), sep = "_")
all(file.exists(files))
tx2gene <- read.table("/home/ferenkagan/Documents/T_transversa/quantification/tx2gene.tsv", sep = "\t")
txi_ttra <- tximport(files, type = "salmon", tx2gene = tx2gene)

##### CREATE ANNOTATION TABLE #####
#Load UniProt annotation crossreferences
crosref_names <- read.csv("~/Documents/T_transversa/ttra_annot.csv", header = F) #name crossreference
crosref <- read.csv("~/Documents/T_transversa//ttra_annot_id.csv", header = F) #ID crosreference
crosref$V2 <- str_remove_all(str_remove_all(str_remove_all(crosref$V2, "sp|"), "^\\|"), "\\|.+$")
uniprot_entrez <- read.table("~/Documents/uniprot-entrez.tsv", header = F, sep = "\t") #uniprot ID to entrez 
crosref$V3 <- uniprot_entrez$V2[match(crosref$V2, uniprot_entrez$V1)]
crosref$V4 <- crosref_names$V2[match(crosref$V1, crosref_names$V1)]

stages <- read.csv("/home/ferenkagan/Documents/T_transversa/quantification/stages.csv", header = F)

#Results from reciprocal best hit blast against human proteome
ensembl <- read.csv("~/Documents/T_transversa/rbh_ttra_stripped.tbl", header = F)
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
go_annot <- read_tsv("~/Documents/Gene_expr_evol/GO/Annotations/ttra_go.tsv") %>% filter(ARGOT_PPV >= 0.5)
go_annot$qpid <- str_remove_all(go_annot$qpid, "t[0-9]+$")
go_annot$goid <- paste0("GO:", go_annot$goid)
go_annot <- go_annot %>% dplyr::select(qpid, goid, desc)

##### BUILD DESEQ OBJECT ######
coldata_ttra <- data.frame(row.names=colnames(txi_ttra$abundance), 
                           Stages_ttra = factor(str_remove_all(colnames(txi_ttra$abundance), "_F[0-9]+")))
dds_ttra <- DESeqDataSetFromTximport(txi_ttra, colData = coldata_ttra, design = ~Stages_ttra)

##### LOW COUNT FILTER #####
dds_ttra <- estimateSizeFactors(dds_ttra)
keep <- rowSums(counts(dds_ttra, normalized = T) >= 10) >= 2 #Atleast 2 columns have more than 10 counts
dds_ttra <- dds_ttra[keep,]

##### VARIANCE STABILIZING #####

dds_rld <- rlogTransformation(dds_ttra)
dds_vst <- vst(dds_ttra)

(p <- plotPCA(dds_vst, intgroup=c("Stages_ttra")))
meanSdPlot(assay(dds_vst))

##### BATCH CORRECTION #####
#Following SVA vignette steps
dds_ttra <- DESeq(dds_ttra)

dat <- counts(dds_ttra, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~Stages_ttra, colData(dds_ttra))
mod0 <- model.matrix(~1, colData(dds_ttra))
svseq <- svaseq(dat, mod, mod0)

#3 surrogate variable were found, including this in the modeling
ddssva_ttra <- dds_ttra
SV1 <- svseq$sv[, 1]
SV2 <- svseq$sv[, 2]
SV3 <- svseq$sv[, 3]
colData(ddssva_ttra) <- cbind(colData(ddssva_ttra), SV1)
colData(ddssva_ttra) <- cbind(colData(ddssva_ttra), SV2)
colData(ddssva_ttra) <- cbind(colData(ddssva_ttra), SV3)

design(ddssva_ttra) <- ~Stages_ttra + SV1 + SV2 + SV3
ddssva_ttra <- DESeq(ddssva_ttra)

###### EXTRACT RESULTS #####

#Contrasts for MZT
contrast_ttra_1 <- c("Stages_ttra", "S02", "S01")
contrast_ttra_2 <- c("Stages_ttra", "S03", "S01")
contrast_ttra_3 <- c("Stages_ttra", "S03", "S02")

#Extracting results
res_1 <- results(ddssva_ttra, contrast_ttra_1, alpha = 0.05)
res_2 <- results(ddssva_ttra, contrast_ttra_2, alpha = 0.05)
res_3 <- results(ddssva_ttra, contrast_ttra_3, alpha = 0.05)

#Inspecting results
summary(res_1)
summary(res_2)
summary(res_3)

#Extracting gene IDs where oocyte stages have expression
maternal <- expr.extr(ls(pattern = "txi_....$"))

#Extracting gene IDs which belong to degraded subset
mzt_deg <- unique(c(rownames(subset(res_1, res_1$log2FoldChange <= -2 & res_1$padj <= 0.05)),
                    rownames(subset(res_2, res_2$log2FoldChange <= -2 & res_2$padj <= 0.05)),
                    rownames(subset(res_3, res_3$log2FoldChange <= -2 & res_3$padj <= 0.05))))

#Exporting result tables
write.table(data.frame(X1 = maternal), "~/Documents/T_transversa/Results/maternal.list",
            quote = F, row.names = F, col.names = F)
write.table(data.frame(X1 = mzt_deg), "~/Documents/T_transversa/Results/degr.list",
            quote = F, row.names = F, col.names = F)

#What proportion of the transcriptome is expressed?
length(maternal) / dim(txi_ttra$abundance)[1]

##### RBP #####
#This section was used to analyze RNA binding domain (RBD) abundances within maternal transcriptome

#Read in data
pfam <- rhmmer::read_tblout("/home/ferenkagan/Documents/T_transversa/ttra_v_pfam.tblout")
#Subset for RNA binding the whole table, this is used as background universe for enrichment
pfam$query_name <- str_remove_all(pfam$query_name, "t[0-9]+$")
pfam_rbp <- pfam[str_detect(pfam$description, "RNA"), ]
pfam_rbp <- pfam_rbp[str_detect(pfam_rbp$description, "binding"), ]

#Subset the HMMER annotation wth ttra results
maternal_domains <- subset(pfam, pfam$query_name %in% maternal)
mzt_dge_domains <- subset(pfam, pfam$query_name %in% mzt_deg)
mzt_zga_domains <- subset(pfam, pfam$query_name %in% mzt_zga)
contr_1_domains <- subset(pfam, pfam$query_name %in% rownames(res_5))
contr_2_domains <- subset(pfam, pfam$query_name %in% rownames(res_6))

#Subsetting with grep with search words of RNA and binding
maternal_domains_rbp <- maternal_domains[str_detect(maternal_domains$description, "RNA"), ]
mzt_dge_domains_rbp <- mzt_dge_domains[str_detect(mzt_dge_domains$description, "RNA"), ]
mzt_zga_domains_rbp <- mzt_zga_domains[str_detect(mzt_zga_domains$description, "RNA"), ]
contr_1_domains_rbp <- contr_1_domains[str_detect(contr_1_domains$description, "RNA"), ]
contr_2_domains_rbp <- contr_2_domains[str_detect(contr_2_domains$description, "RNA"), ]

maternal_domains_rbp <- maternal_domains_rbp[str_detect(maternal_domains_rbp$description, "binding"), ]
mzt_dge_domains_rbp <- mzt_dge_domains_rbp[str_detect(mzt_dge_domains_rbp$description, "binding"), ]
mzt_zga_domains_rbp <- mzt_zga_domains_rbp[str_detect(mzt_zga_domains_rbp$description, "binding"), ]
contr_1_domains_rbp <- contr_1_domains_rbp[str_detect(contr_1_domains_rbp$description, "binding"), ]
contr_2_domains_rbp <- contr_2_domains_rbp[str_detect(contr_2_domains_rbp$description, "binding"), ]

#Proportions of RBP at different groupings
length(unique(maternal_domains_rbp$query_name)) / length(unique(maternal_domains$query_name))
length(unique(pfam_rbp$query_name)) / length(unique(pfam$query_name))
length(unique(mzt_dge_domains_rbp$query_name)) / length(unique(mzt_dge_domains$query_name))
length(unique(mzt_zga_domains_rbp$query_name)) / length(unique(mzt_zga_domains$query_name))
length(unique(contr_1_domains_rbp$query_name)) / length(unique(contr_1_domains$query_name))
length(unique(contr_2_domains_rbp$query_name)) / length(unique(contr_2_domains$query_name))

#For DGE genes hypergeometric test
phyper(length(unique(mzt_dge_domains_rbp$query_name)), 
       length(unique(mzt_dge_domains$query_name)), 
       length(unique(pfam$query_name)) - length(unique(mzt_dge_domains$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = T)

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

##### UTR ANALYSIS #####
#This section was used for analyzing UTR motif enrichment and motif occurence results in combination with
#expression data

#Read in enrichment results from the two separate databases
att_crosref <- read.table("~/Documents/Gene_expr_evol/Gene_lengths/ATtRACT_db.txt", sep = "\t", header = T)
sea_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_dge_ray2013.out",
                          sep = "\t", header = T)
sea_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_dge_attract.out",
                          sep = "\t", header = T)
sea_dge_att$ALT_ID <- att_crosref$Gene_name[match(sea_dge_att$ID, att_crosref$Matrix_id)]
sea_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_mat_ray2013.out",
                          sep = "\t", header = T)
sea_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_mat_attract.out",
                          sep = "\t", header = T)
sea_mat_att$ALT_ID <- att_crosref$Gene_name[match(sea_mat_att$ID, att_crosref$Matrix_id)]

#Combine the two database based results
sea_dge <- rbind(sea_dge_att, sea_dge_ray)
sea_mat <- rbind(sea_mat_att, sea_mat_ray)

#Inspect results
tibble(sea_dge) %>% dplyr::arrange(QVALUE, desc(ENR_RATIO))
tibble(sea_mat) %>% dplyr::arrange(QVALUE, desc(ENR_RATIO))

#Import motif occurence results for degraded subset and create a combination table of SEA and FIMO
fimo_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_dge_ray2013.fimo",
                           sep = "\t", header = T)
fimo_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_dge_attract.fimo",
                           sep = "\t", header = T)
fimo_dge <- rbind(fimo_dge_att, fimo_dge_ray)
fimo_dge <- subset(fimo_dge, fimo_dge$motif_id %in% sea_dge$ID)
fimo_dge$motif_alt_id <- sea_dge$ALT_ID[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sequence_alt_name <- crosref$V4[match(fimo_dge$sequence_name, crosref$V1)]
fimo_dge$enr_val <- sea_dge$ENR_RATIO[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sea.q.value <- sea_dge$QVALUE[match(fimo_dge$motif_id, sea_dge$ID)]

#Import motif occurence results for maternal subset and create a combination table of SEA and FIMO
fimo_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_mat_ray2013.fimo",
                           sep = "\t", header = T)
fimo_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Ttra_mat_attract.fimo",
                           sep = "\t", header = T)

#Combine result tables
fimo_mat <- rbind(fimo_mat_att, fimo_mat_ray)
fimo_mat <- subset(fimo_mat, fimo_mat$motif_id %in% sea_mat$ID)

#Format result table
fimo_mat$motif_alt_id <- sea_mat$ALT_ID[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sequence_alt_name <- crosref$V4[match(fimo_mat$sequence_name, crosref$V1)]
fimo_mat$enr_val <- sea_mat$ENR_RATIO[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sea.q.value <- sea_mat$QVALUE[match(fimo_mat$motif_id, sea_mat$ID)]

fimo_dge <- tibble(fimo_dge) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)
fimo_mat <- tibble(fimo_mat) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)

#Inspect motif occurence for each transcript
View(tibble(fimo_mat) %>% filter(q.value <= 0.05) %>% group_by(sequence_name, motif_alt_id, sequence_alt_name) %>% dplyr::summarise(n()))
View(tibble(fimo_dge) %>% filter(q.value <= 0.05) %>% group_by(sequence_name, motif_alt_id, sequence_alt_name) %>% dplyr::summarise(n()))

#Inspect significant motif enrichments
unique((tibble(fimo_mat) %>% filter(q.value <= 0.05))$motif_alt_id)

#Subsetting expression data for plotting target gene expressions
temp <- t(scale(t(collapseReplicate(assay(dds_vst)))))

ids <- subset(fimo_mat, fimo_mat$q.value <= 0.05)
ids <- subset(ids, grepl("ELAV", ids$motif_alt_id))

rownames(temp) <- rownames(assay(dds_rld))
temp_subs <- subset(temp, rownames(temp) %in% unique(ids$sequence_name))
#temp_subs <- subset(temp_subs, rownames(temp_subs) %in% mzt_deg)
colnames(temp_subs) <- stages$V2

#svg("~/Desktop/Figures/fimo_heatmap_ttra.svg")
#png("~/Desktop/Figures/fimo_heatmap_ttra.png",
#    width=750, height=1000)
print(ComplexHeatmap::Heatmap(matrix = temp_subs,
                              column_title = "ELAV motif containing genes",
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
ids <- subset(fimo_mat, fimo_mat$q.value <= 0.05)
variables <- unique(ids$motif_alt_id)
utr_list <- vector(mode = "list", length = length(variables))
resrbind <- rbind(as.data.frame(res_1), as.data.frame(res_2))
resrbind <- rbind(resrbind, as.data.frame(res_3))
resrbind <- resrbind[unique(rownames(resrbind)), ]

#Looping through RBP target data and performing enrichment
for(k in 1: length(variables)){
  temp_name <- variables[k]
  temp_ids <- subset(ids, ids$motif_alt_id == variables[k])  #Enrichment table subsetting for target motif
  
  temp_res <- subset(resrbind, rownames(resrbind) %in% unique(temp_ids$sequence_name)) #Subsetting expression data according to motif targets
  temp_res <- temp_res[order(temp_res$log2FoldChange, decreasing = T),]

  temp_ego<- enricher(rownames(temp_res), TERM2GENE = go_annot[,c(2,1)], TERM2NAME = go_annot[, c(2,3)])
  temp_ego@result <- temp_ego@result[temp_ego@result$qvalue <= 0.05, ]
  
  #If no enrichment found to avoid error message just skip
  if(is.null(temp_ego@result)){next}
  
  utr_list[[k]] <- temp_ego@result
  names(utr_list)[[k]] <- temp_name
  
}


##### MFUZZ CLUSTERING #####
#This section was used for fuzzy clustering of gene expression data and GO enrichment on resultant clusters
library(Mfuzz)

#Collapsing replicates
cts <- DESeq2::collapseReplicates(ddssva_ttra, groupby = coldata_ttra$Stages_ttra)
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

#Export plot for each cluster
for(i in 1:length(unique(cts_clst$cluster))){
  svg(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_ttra", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_ttra", ".png"))),
  #  width=550, height=550)
  mfuzz.plot2(cts, cts_clst, min.mem = 0.5, colo = 'fancy', time.labels = stages$V2, 
              centre = T, single = i, x11 = F)
  mtext(paste0("Cluster size: ", sum(cts_clst$cluster %in% i)))
   dev.off()
}

#Get cluster cores for enrichment analysis of those
cores <- acore(cts, cts_clst)

#Initializing variables
variables <- unique(cts_clst$cluster)
go_list <- vector(mode = "list", length = length(variables))
resrbind <- rbind(as.data.frame(res_1), as.data.frame(res_2))
resrbind <- rbind(resrbind, as.data.frame(res_3))
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
#This section was used to query KEGG pathways within Ttra transcriptome
library(KEGGREST)

rld <- collapseReplicate(as.data.frame(assay(dds_rld)))

#These pathways are queried, for further information visit: https://www.genome.jp/kegg/pathway.html
keggid <- c("04010", "04014", "04310", "04330", "04340", "04350", "04390", "04630", "04064",
            "04668", "04068", "04020", "04070", "04072", "04024", "04151", "04514", "04512",
            "04210", "04360")

#Initializing variables
expr_list <- vector(mode = "list", length = length(keggid))
missing <- c()

#Looping through KEGG pathways
for(i in 1:length(keggid)){
  #Getting pathway name
  temp_name <- keggFind(database = "pathway", query = keggid[i])
  
  #Getting genes involved in pathway
  temp_id <- keggGet(paste("hsa", keggid[i], sep =  ""))[[1]]$GENE
  temp_id_df <- data.frame(Entrez = temp_id[seq(2,length(temp_id),2)-1],
                           Name = str_remove_all(temp_id[seq(2,length(temp_id),2)], ";.*"))
  
  #What is the missing proportion of genes
  missing[i] <- 1 - (dim(crosref[crosref$EntrezIDComb %in% temp_id_df$Entrez, ])[1] / 
                       length(temp_id_df$Entrez))
  names(missing)[i] <- temp_name
  
  #Build expression table
  temp_id <- rownames(rld)
  temp_id <- crosref$EntrezIDComb[match(temp_id, crosref$V1)]
  mat_temp <- as.data.frame(rld[temp_id %in% temp_id_df$Entrez, ])
  
  #Attach name information to ID
  temp_id_gn <- crosref$EntrezIDComb[match(rownames(mat_temp), crosref$V1)]
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
  colnames(temp) <- stages$V2
  
  ha <- rowAnnotation(foo = anno_mark(at = match(unique(str_remove_all(rownames(temp), "Ttra[0-9]+_")), str_remove_all(rownames(temp), "Ttra[0-9]+_")), 
                                      labels = rownames(temp[match(unique(str_remove_all(rownames(temp), "Ttra[0-9]+_")), str_remove_all(rownames(temp), "Ttra[0-9]+_")), ])))
  
  
  #svg(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("_ttra", ".svg"))))
  png(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("_ttra", ".png"))),
      width=750, height=1000)
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
  dev.off()
  
}

###### PLOT INDIVIDUAL GENES #####
#This section was used for plotting gene expression values across development
library(wesanderson)

#rlog transformed data
rld <- as.data.frame(assay(dds_rld))

#Target genes
targets <- c("Ttra000521", "Ttra001758", "Ttra000600", "Ttra017865",
             "Ttra011062","Ttra010320", "Ttra013696",
             "Ttra000136", "Ttra000079", "Ttra008558", "Ttra015656",
             "Ttra000102", "Ttra003215", "Ttra002432", "Ttra001066",
             "Ttra017865", "Ttra000600", "Ttra001827")
stages_ttra <- data.frame(X1 = stages$V2,
                          X2 = c(paste0(rep("stage", 14), 1:14)))

#Loopin through each target gene
for(i in 1:length(targets)){
  #Could be missing from table
  if(!(targets[i] %in% rownames(rld))){next}
  
  #Subsetting for gene and formatting names
  expr <- as.data.frame(t(subset(rld, rownames(rld) %in% targets[i])))
  colnames(expr) <- "expr"
  expr$replicate <- str_remove_all(rownames(expr), "S[0-9]+_F")
  expr$staging <- factor(stages$V2[match(str_remove_all(rownames(expr), "_F[0-9]+"), stages$V1)],
                         levels = stages$V2)
  
  #svg(paste0("~/Desktop/Figures/", paste0(crosref$V4[crosref$V1 == targets[i]], paste0("_ttra", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(crosref$V4[crosref$V1 == targets[i]], paste0("_ttra", ".png"))),
  #    width=550, height=550)
  print(ggplot(expr, aes(x = staging, y = expr, color = replicate, group = replicate)) +
          geom_point(size = 6) +
          theme_bw() +
          ggtitle(crosref_names$V2[crosref_names$V1 %in% targets[i]]) +
          scale_color_manual(values=wes_palettes$Rushmore[seq(1, length(wes_palettes), by = 2)]) +
          xlab("") +
          ylab("Regularized log transformed expression") +
          theme(text = element_text(size = 18),
                axis.text.x = element_text(angle = 45, hjust = 1)))
  #dev.off()
  
}
