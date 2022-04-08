#Script for Differential gene expression analysis
#18.05.2021
#Ferenc.Kagan@uib.no

##### LIBRARIES ######
libraries <- c("DESeq2", "pheatmap", "vsn", "hexbin", "cowplot", "apeglm", "tximport", "gprofiler2",
               "viridis", "PoiClaClu", "genefilter", "vidger", "reshape2", "ggplot2", "stringr", 
               "gplots", "sva", "tidyverse", "biomaRt", "RColorBrewer", "gprofiler2", "DEGreport",
               "tximport", "enrichplot", "DOSE", "org.Hs.eg.db", "plotly", "clusterProfiler", "gprofiler2",
               "ape", "motmot", "geiger", "phytools", "OUwie", "ggtree", "casper", "taxizedb", "simplifyEnrichment",
               "ggsci", "Biostrings", "rtracklayer", "GenomicFeatures", "rstatix", "ggpubr", "BSgenome", "zFPKM")
lapply(libraries, library, character.only = TRUE)

##### FUNCTION #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

##### CROSREFERENCE TABLE FORMATTING #####
#Crosreference tables containing 1-to-1 orthologs and genome ID crosreferences within P.caudatus
crosref <- read.table("/home/ferenkagan/Documents/Pcau_project/quantification/Genome_v2/one_to_one_crosref.txt")
crosref2 <- read.table("/home/ferenkagan/Documents/Pcau_project/quantification/Genome_v2/crossref_pcau.txt")

#Building crossreference table of diverse IDs
crosref2$V1 <- str_remove_all(crosref2$V1, "\\.[0-9]+")
crosref$V1 <- str_remove_all(crosref$V1, "\\.[0-9]+")
crosref$V3 <- crosref2$V1[match(crosref$V2, crosref2$V2)]

#Getting Entrez IDs for Enesmbl IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#View(listFilters(mart))
entrez_tbl <- getBM(attributes = c("ensembl_gene_id",  "entrezgene_id", "external_gene_name"),  
                    values = str_remove_all(crosref$V1, "\\.[0-9]+$"), bmHeader = T, mart = mart)

#Attaching Entrez IDs to crosreference table
crosref$V4 <- entrez_tbl$`NCBI gene (formerly Entrezgene) ID`[match(crosref$V1, entrez_tbl$`Gene stable ID`)]
crosref$V5 <- entrez_tbl$`Gene name`[match(crosref$V1, entrez_tbl$`Gene stable ID`)]

#Merging two crosreference tables
crosref <- merge(crosref2, crosref, by.x = "V1", by.y = "V3", all.x = T)
colnames(crosref) <- c("GeneID", "GeneName", "EnsemblID", "EntrezName", "EntrezID", "ExternalName")
crosref$GeneName_edit <- str_remove_all(crosref$GeneName, "-[0-9]+") #Removing pattern as it seems to affect later conversions

#Using gene symbols from annotation and gconvert getting more entrez IDs
gConv <- gconvert(crosref$GeneName_edit, filter_na = T,
                  target = "ENTREZGENE_ACC")[, c(2, 4)]
gConv <- as.data.frame(apply(gConv, 2, function(x){str_replace_all(x, "nan", NA_character_)}))

#Combine the orthology inference derived IDs and annotation derived IDs into one table
crosref <- unique(merge(crosref, gConv, by.x = "GeneName_edit", by.y = "input", all = T))
crosref <- tibble(crosref) %>%                                            #Combining two columns as 
  mutate(EntrezIDComb = coalesce(EntrezID, as.integer(target))) %>%       #they have some IDs which are
  dplyr::select(GeneName, GeneName_edit, GeneID, EntrezIDComb, ExternalName, EnsemblID)  #missing from the other one  

crosref <- readRDS("/home/ferenkagan/Documents/Pcau_project/quantification/Genome_v2/Intermediate_files/combined_crosref.RDS")

##### LOAD SALMON COUNTS #####
files_ext <- list.files(path = "~/Documents/Pcau_project/quantification/Genome_v2",  pattern = "*salmon")
files <- file.path("~/Documents/Pcau_project/quantification/Genome_v2/", files_ext, "quant.sf")
names(files) <- c(paste("stage1", 1:2, sep = "_"), paste("stage2", 1:2, sep = "_"), paste("stage3", 1:2, sep = "_"),
                  paste("stage4", 1:2, sep = "_"), paste("stage5", 1:2, sep = "_"))
all(file.exists(files))
tx2gene <- read.table("~/Documents/Pcau_project/quantification/Genome_v2/tx2gene.tsv", header = F)
txi_pcau <- tximport(files, type = "salmon", tx2gene = tx2gene, infRepStat = rowMeans2)

##### BUILD DESEQ OBJECT ######
coldata_pcau <- data.frame(row.names=colnames(txi_pcau$abundance), 
                           Stages_pcau = factor(str_remove_all(colnames(txi_pcau$abundance), "_[0-9]+")))
dds <- DESeqDataSetFromTximport(txi_pcau,  colData = coldata_pcau, design = ~Stages_pcau)

##### LOW COUNT FILTER #####
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = T) >= 10) >= 2
dds <- dds[keep,]

##### VARIANCE STABILIZING #####
vst_transf <- vst(dds) 
rld_transf <- rlogTransformation(dds)

(p <- plotPCA(rld_transf, intgroup=c("Stages_pcau")) +
  theme_bw() +
  theme(text = element_text(size = 18)))
meanSdPlot(assay(rld_transf))

##### BATCH CORRECTION #####
dds <- DESeq(dds)

dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~Stages_pcau, colData(dds))
mod0 <- model.matrix(~1, colData(dds))
svseq <- svaseq(dat, mod, mod0)

#2 surrogate variable were found, including this in the modeling
ddssva_pcau <- dds
colData(ddssva_pcau) <- cbind(colData(ddssva_pcau), svseq$sv[,1], svseq$sv[,2])
colnames(colData(ddssva_pcau))[2:3] <- c("SV1", "SV2")
design(ddssva_pcau) <- ~Stages_pcau + SV1 + SV2
ddssva_pcau <- DESeq(ddssva_pcau)  

###### EXTRACT RESULTS #####
#Contrasts for MZT and extracting results
contrast_pcau_1 = c("Stages_pcau", "stage2", "stage1")
contrast_pcau_2 = c("Stages_pcau", "stage3", "stage2")

res_pcau_1 <- results(ddssva_pcau, contrast_pcau_1, alpha = 0.05)
res_pcau_2 <- results(ddssva_pcau, contrast_pcau_2, alpha = 0.05)

res_0v1 <- subset(res_pcau_1, padj <= 0.05 & abs(log2FoldChange) >= 2) 
res_1v2 <- subset(res_pcau_2, padj <= 0.05 & abs(log2FoldChange) >= 2) 

#Adding annotation information to result table
res_0v1$Names <- crosref$EntrezIDComb[match(rownames(res_0v1), crosref$GeneID)]
res_1v2$Names <- crosref$EntrezIDComb[match(rownames(res_1v2), crosref$GeneID)]
res_0v1$GeneNames <- crosref$GeneName_edit[match(rownames(res_0v1), crosref$GeneID)]
res_1v2$GeneNames <- crosref$GeneName_edit[match(rownames(res_1v2), crosref$GeneID)]

#Subsetting for only degraded subset of genes
res_0v1_down <- subset(res_0v1, res_0v1$log2FoldChange <= -2)
res_1v2_down <- subset(res_1v2, res_1v2$log2FoldChange <= -2)

down <- as.data.frame(rbind(res_0v1_down, res_1v2_down))
mzt_deg <- unique(rownames(down))

#Defining non-degraded subset of the maternal transcriptome
maternal <- txi_pcau$abundance[, 1:2]
maternal <- rownames(txi_pcau$abundance[rowMeans(maternal) >= 2, ])

##### UTR ANALYSIS #####
#This section was used for analyzing UTR motif enrichment and motif occurence results in combination with
#expression data

#Read in enrichment results from the two separate databases
att_crosref <- read.table("~/Documents/Gene_expr_evol/Gene_lengths/ATtRACT_db.txt", sep = "\t", header = T)
sea_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_dge_ray2013.out",
                          sep = "\t", header = T)
sea_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_dge_attract.out",
                          sep = "\t", header = T)
sea_dge_att$ALT_ID <- att_crosref$Gene_name[match(sea_dge_att$ID, att_crosref$Matrix_id)]
sea_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_mat_ray2013.out",
                          sep = "\t", header = T)
sea_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_mat_attract.out",
                          sep = "\t", header = T)
sea_mat_att$ALT_ID <- att_crosref$Gene_name[match(sea_mat_att$ID, att_crosref$Matrix_id)]

#Combine the two database based results
sea_dge <- rbind(sea_dge_att, sea_dge_ray)
sea_mat <- rbind(sea_mat_att, sea_mat_ray)

#Read in motif occurence results for the two databases for the degraded subset of genes
fimo_dge_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_dge_ray2013.fimo",
                          sep = "\t", header = T)
fimo_dge_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_dge_attract.fimo",
                          sep = "\t", header = T)

#Combine the two results with the results of simple enrichment results
fimo_dge <- rbind(fimo_dge_att, fimo_dge_ray)
fimo_dge <- subset(fimo_dge, fimo_dge$motif_id %in% sea_dge$ID)

#Format result table
fimo_dge$motif_alt_id <- sea_dge$ALT_ID[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sequence_alt_name <- crosref$GeneName_edit[match(fimo_dge$sequence_name, crosref$GeneID)]
fimo_dge$enr_val <- sea_dge$ENR_RATIO[match(fimo_dge$motif_id, sea_dge$ID)]
fimo_dge$sea.q.value <- sea_dge$QVALUE[match(fimo_dge$motif_id, sea_dge$ID)]

#Read in motif occurence results for the two databases for the non-degraded subset of genes
fimo_mat_ray <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_mat_ray2013.fimo",
                           sep = "\t", header = T)
fimo_mat_att <- read.table("/home/ferenkagan/Documents/Gene_expr_evol/Gene_lengths/Results/Pcau_mat_attract.fimo",
                           sep = "\t", header = T)

#Combine the two results with the results of simple enrichment results
fimo_mat <- rbind(fimo_mat_att, fimo_mat_ray)
fimo_mat <- subset(fimo_mat, fimo_mat$motif_id %in% sea_mat$ID)

#Format result table
fimo_mat$motif_alt_id <- sea_mat$ALT_ID[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sequence_alt_name <- crosref$GeneName_edit[match(fimo_mat$sequence_name, crosref$GeneID)]
fimo_mat$enr_val <- sea_mat$ENR_RATIO[match(fimo_mat$motif_id, sea_mat$ID)]
fimo_mat$sea.q.value <- sea_mat$QVALUE[match(fimo_mat$motif_id, sea_mat$ID)]

#Inspect motif occurence for each transcript
fimo_dge <- tibble(fimo_dge) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)
fimo_mat <- tibble(fimo_mat) %>% select(-strand, -score, -start, -stop, -p.value, -matched_sequence)

temp <- t(scale(t(collapseReplicate(assay(rld_transf)))))

ids <- subset(fimo_mat, fimo_mat$q.value <= 0.05)
ids <- subset(ids, ids$motif_alt_id %in% c("CELF1"))

#Subsetting expression data for plotting target gene expressions
rownames(temp) <- rownames(assay(rld_transf))
temp_subs <- subset(temp, rownames(temp) %in% unique(ids$sequence_name))
#temp_subs <- subset(temp_subs, rownames(temp_subs) %in% mzt_deg)
colnames(temp_subs) <- c("oocyte", "cleavage", "gastrula", "introvertula", "larva")

#svg("~/Desktop/Figures/fimo_heatmap_pcau.svg")
#png("~/Desktop/Figures/fimo_heatmap_pcau.png",
#    width=750, height=1000)
print(ComplexHeatmap::Heatmap(matrix = temp_subs,
                              show_row_names = F,
                              show_row_dend = F,
                              row_dend_reorder = F,
                              cluster_columns = F,
                              col = rev(brewer.pal(9,"RdBu")),
                              column_names_rot = 45, column_names_gp = gpar(fontface = "italic", fontsize = 22),
                              name = "Gene Z-score"))
#dev.off()

#GO enrichment of motif target genes
#Preparing data
ids <- subset(fimo_dge, fimo_dge$q.value <= 0.05)
variables <- unique(ids$motif_alt_id)
utr_list <- vector(mode = "list", length = length(variables))
resrbind <- rbind(as.data.frame(res_0v1), as.data.frame(res_1v2))
resrbind <- rbind(resrbind, as.data.frame(res_2v3))
resrbind <- rbind(resrbind, as.data.frame(res_3v4))
resrbind <- resrbind[unique(rownames(resrbind)), ]


#Looping through RBP target data and performing enrichment
for(k in 1: length(variables)){
  temp_name <- variables[k]
  temp_ids <- subset(ids, ids$motif_alt_id == variables[k]) #Enrichment table subsetting for target motif
  
  temp_res <- resrbind[order(abs(resrbind$log2FoldChange), decreasing = T),]
  temp_res <- subset(temp_res, rownames(temp_res) %in% unique(temp_ids$sequence_name)) #Subsetting expression data according to motif targets

  #Perform enrichment analysis
  temp_ego <- enrichGO(gene = temp_res$Names, OrgDb = org.Hs.eg.db, 
                       ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
  temp_ego@result <- temp_ego@result[temp_ego@result$qvalue <= temp_ego@qvalueCutoff, ]
  
  #If no enrichment found to avoid error message just skip
  if(is.null(temp_ego@result)){next}

  utr_list[[k]] <- temp_ego@result
  names(utr_list)[[k]] <- temp_name
  
}

#ELAV protein expression plots
#Fetch ELAV expressions
asd <- as.data.frame(subset(assay(rld_transf), rownames(assay(rld_transf)) %in% c("PCA14917", "PCA15855", "PCA27694")))
asd$gene <- rownames(asd)

#Format table
asd_long <- gather(as.data.frame(asd), stage, expr, stage1_1:stage5_2)
asd_long$alpha <- asd_long$gene == "PCA27694" #used for transparency 
asd_long$stage <- str_remove_all(asd_long$stage, "_[0-9]+$")

#svg("~/Desktop/Figures/pcau_elav_genes.svg")
print(ggplot(asd_long, aes(x = stage, y = expr, color = gene, alpha = alpha)) +
  geom_point(size = 5) + theme_bw() + 
  scale_alpha_discrete(range = c(0.3, 1)) + 
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#33A02C")) +
  theme(text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  xlab(NULL) +
  ylab("Regularized log transformed expression") +
  scale_x_discrete(labels = c("oocyte", "cleavage", "gastrula", "introvertula", "larva")))
#dev.off()

#ELAVL4 motif occurences
fimo_mat_means <- fimo_mat %>% 
  group_by(motif_alt_id, sequence_name) %>% 
  dplyr::summarize(n = n()) %>% 
  filter(motif_alt_id %in% c("ELAVL4"))

#svg("~/Desktop/Figures/pcau_elav_motifnr.svg")
print(ggplot(fimo_mat_means, aes(x = n, fill = motif_alt_id)) + 
  geom_bar(alpha = 0.8) +
  theme_bw() +
  xlab(NULL) +
  theme(text = element_text(size = 22), legend.position = 'none') +
  scale_fill_manual(values = "#33A02C"))
#dev.off()

##### RBP #####
#This section was used to analyze RNA binding domain (RBD) abundances within maternal transcriptome

#Read in data
pfam <- rhmmer::read_tblout("/home/ferenkagan/Documents/Pcau_project/quantification/Genome_v2/PriCau_pfam.tblout")
pfam$query_name <- str_remove_all(pfam$query_name, "\\.[0-9]+")

#Subset for RNA binding the whole table, this is used as background universe for enrichment
pfam_rbp <- pfam[str_detect(pfam$description, "RNA"), ]
pfam_rbp <- pfam_rbp[str_detect(pfam_rbp$description, "binding"), ]
  
#Subset the HMMER annotation
maternal_domains <- subset(pfam, pfam$query_name %in% maternal)
mzt_dge_domains <- subset(pfam, pfam$query_name %in% mzt_deg)
mzt_zga_domains <- subset(pfam, pfam$query_name %in% mzt_zga)

#Subset for uncharacterized genes
unchar_maternal <- subset(maternal_domains, maternal_domains$query_name %in% crosref$GeneID[grep("Unchar", crosref$GeneName)])
unchar_mzt_dge <- subset(mzt_dge_domains, mzt_dge_domains$query_name %in% crosref$GeneID[grep("Unchar", crosref$GeneName)])

#Subsetting with grep with searchwords of RNA and binding
maternal_domains_rbp <- maternal_domains[str_detect(maternal_domains$description, "RNA"), ]
mzt_dge_domains_rbp <- mzt_dge_domains[str_detect(mzt_dge_domains$description, "RNA"), ]
contr_1_domains_rbp <- contr_1_domains[str_detect(contr_1_domains$description, "RNA"), ]
contr_2_domains_rbp <- contr_2_domains[str_detect(contr_2_domains$description, "RNA"), ]

maternal_domains_rbp <- maternal_domains_rbp[str_detect(maternal_domains_rbp$description, "binding"), ]
mzt_dge_domains_rbp <- mzt_dge_domains_rbp[str_detect(mzt_dge_domains_rbp$description, "binding"), ]
unchar_maternal_rbp <- unchar_maternal_rbp[str_detect(unchar_maternal_rbp$description, "binding"), ]
unchar_mzt_dge_rbp <- unchar_mzt_dge_rbp[str_detect(unchar_mzt_dge_rbp$description, "binding"), ]

#Proportions of RBP at different groupings
length(unique(maternal_domains_rbp$query_name)) / length(unique(maternal_domains$query_name))
length(unique(pfam_rbp$query_name)) / length(unique(pfam$query_name))
length(unique(mzt_dge_domains_rbp$query_name)) / length(unique(mzt_dge_domains$query_name))
length(unique(unchar_maternal_rbp$query_name)) / length(unique(unchar_maternal$query_name))
length(unique(unchar_mzt_dge_rbp$query_name)) / length(unique(unchar_mzt_dge$query_name))

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
phyper(length(unique(unchar_maternal_rbp$query_name)), 
       length(unique(unchar_maternal$query_name)), 
       length(unique(pfam$query_name)) - length(unique(unchar_maternal$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)
phyper(length(unique(unchar_mzt_dge_rbp$query_name)), 
       length(unique(unchar_mzt_dge$query_name)), 
       length(unique(pfam$query_name)) - length(unique(unchar_mzt_dge$query_name)),
       length(unique(pfam_rbp$query_name)), lower.tail = F)


#Plotting expression distribution for RBD containing genes
library(wesanderson)

df <- data.frame(mean = c(rowMeans(txi_pcau$abundance[,1:2]), 
                    rowMeans(txi_pcau$abundance[rownames(txi_pcau$abundance) %in% mzt_dge_domains_rbp$query_name, 1:2]),
                    rowMeans(txi_pcau$abundance[rownames(txi_pcau$abundance) %in% maternal_domains_rbp$query_name, 1:2])),
           category = c(rep("reference", length(rowMeans(txi_pcau$abundance[,1:2]))),
                        rep("degraded", length(rowMeans(txi_pcau$abundance[rownames(txi_pcau$abundance) %in% mzt_dge_domains_rbp$query_name, 1:2]))),
                        rep("maternal", length(rowMeans(txi_pcau$abundance[rownames(txi_pcau$abundance) %in% maternal_domains_rbp$query_name, 1:2])))
                                                ))
#svg("~/Desktop/Figures/pcau_rbd_expression.svg")
print(ggplot(df, aes(x = log(mean), fill = category)) + 
  geom_density(alpha = 0.8) +
  theme_bw() +
  scale_fill_manual(values = wes_palettes$GrandBudapest1) + 
  theme(text = element_text(size = 22),
        legend.position = 'none'))
#dev.off()

ggplot(df, aes(x = stage, y = expr, fill = gene)) +
  geom_dotplot()
 
##### MFUZZ CLUSTERING #####
#This section was used for fuzzy clustering of gene expression data and GO enrichment on resultant clusters
library(Mfuzz)

#Collapsing replicates
cts <- DESeq2::collapseReplicates(ddssva_pcau, groupby = coldata_pcau$Stages_pcau)
cts <- counts(cts, normalized = T)

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
c <- Dmin(cts, m = m, crange = seq(4,50,2), repeats = 5)

#Fuzzy clustering, centers determined by Dmin plot inspection
cts_clst <- mfuzz(cts, centers = 10, m = m)

stages <- c("oocyte", "cleavage", "gastrula", "introvertula", "larva")

#Export plot for each cluster
for(i in 1:length(unique(cts_clst$cluster))){
  svg(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_pcau", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(paste0("cluster", i), paste0("_pcau", ".png"))),
  #  width=550, height=550)
  mfuzz.plot2(cts, cts_clst, min.mem = 0.5, colo = 'fancy', time.labels = c("oocyte", "cleavage", "gastrula", "introvertula", "larva"),
              centre = T, single = i, x11 = F)
  mtext(paste0("Cluster size: ", sum(cts_clst$cluster %in% i)))
  dev.off()
}

#Get cluster cores for enrichment analysis of those
cores <- acore(cts, cts_clst)

#Initialize variables
variables <- unique(cts_clst$cluster)
go_list <- vector(mode = "list", length = length(variables))
resrbind <- rbind(as.data.frame(res_0v1), as.data.frame(res_1v2))
resrbind <- rbind(resrbind, as.data.frame(res_2v3))
resrbind <- rbind(resrbind, as.data.frame(res_3v4))
resrbind <- resrbind[unique(rownames(resrbind)), ]

for(k in 1: length(variables)){
  temp_clust <- cores[[k]]$NAME
  temp_name <- paste0("Cluster_", k)
  
  #Subeset for entrezIDs  
  temp_id <- crosref$EntrezIDComb[crosref$GeneID %in% temp_clust] 
   
  #Perform unviersal enrichment analysis
  temp_ego <- enrichGO(gene = temp_id, 
                       universe = crosref$EntrezIDComb[crosref$GeneID %in% rownames(txi_pcau$abundance)],
                       OrgDb = org.Hs.eg.db, 
                       ont = "ALL", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
  temp_ego@result <- temp_ego@result[temp_ego@result$qvalue <= temp_ego@qvalueCutoff, ]
    
  #If no enrichment found to avoid error message just skip
  if(is.null(temp_ego@result)){next}
  
  go_list[[k]] <- temp_ego@result
  names(go_list)[[k]] <- temp_name

}


##### BUILDING EXPRESSION TABLES ######
#This section was used to query KEGG pathways within Pcau transcriptome
library(KEGGREST)

rld <- collapseReplicate(as.data.frame(assay(rld_transf)))

#These pathways are queried, for further information visit: https://www.genome.jp/kegg/pathway.html
keggid <- c("04010", "04014", "04310", "04330", "04340", "04350", "04390", "04630", "04064",
            "04668", "04068", "04020", "04070", "04072", "04024", "04151", "04514", "04512",
            "04210", "04360")

#Initialize variables
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
  missing[i] <- 1 - (length(crosref$EntrezIDComb[crosref$EntrezIDComb %in% temp_id_df$Entrez]) / 
                       length(temp_id_df$Entrez))
  names(missing)[i] <- temp_name
  
  #Build expression table
  temp_id <- rownames(rld)
  temp_id <- crosref$EntrezIDComb[match(temp_id, crosref$GeneID)]
  mat_temp <- as.data.frame(rld[temp_id %in% temp_id_df$Entrez, ])
  temp_row <- paste(rownames(mat_temp), 
                    temp_id_df$Name[match(crosref$EntrezIDComb[match(rownames(mat_temp), crosref$GeneID)], temp_id_df$Entrez)],
                    sep = "_")
  
  rownames(mat_temp) <- temp_row
  
  
  #Assign exprs table to list
  expr_list[[i]] <- mat_temp
  names(expr_list)[i] <- temp_name
}

#Plotting
anno <- data.frame(Stages = colnames(expr_list[[1]]))
#rownames(anno) <- colnames(mat_temp)[-length(colnames(mat_temp))]

for(i in 1:length(expr_list)){
  temp <-  apply(expr_list[[i]], 2, function(x){as.double(x)})
  temp <- t(scale(t(temp)))
  
  rownames(temp) <- rownames(expr_list[[i]])
  colnames(temp) <- c("oocyte", "cleavage", "gastrula", "introvertula", "larva")
  
  ha <- rowAnnotation(foo = anno_mark(at = match(unique(str_remove_all(rownames(temp), "PCA[0-9]+_")), str_remove_all(rownames(temp), "PCA[0-9]+_")), 
                                      labels = rownames(temp[match(unique(str_remove_all(rownames(temp), "PCA[0-9]+_")), str_remove_all(rownames(temp), "PCA[0-9]+_")), ])))
  
  
  #svg(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("pcau", ".svg"))))
  png(paste0("~/Desktop/Figures/", paste0(names(expr_list)[i], paste0("pcau", ".png"))),
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

##### INDIVIDUAL GENE PLOTS #####
#This section was used for plotting gene expression values across development

#rlog transformed data
rld <- as.data.frame(assay(rld_transf))

#Target genes
targets <- c("PCA25350", "PCA16978", "PCA23158", "PCA00491", 
             "PCA01390", "PCA27210", "PCA27215","PCA30505",
             "PCA16393", "PCA12942", "PCA08236", "PCA09112",
             "PCA24948", "PCA24830", "PCA26146", "PCA26273",
             "PCA26447", "PCA03101", "PCA11916", "PCA00491")

#Loopin through each target gene
for(i in 1:length(targets)){
  #Could be missing from table
  if(!(targets[i] %in% rownames(rld))){next}
  
  #Subsetting for gene and formatting names
  expr <- as.data.frame(t(subset(rld, rownames(rld) %in% targets[i])))
  colnames(expr) <- "expr"
  expr$replicate <- str_remove_all(rownames(expr), "stage[0-9]+_")
  expr$staging <- factor(c(rep("oocyte", 2), rep("cleavage", 2), rep("gastrula", 2), rep("introvertula", 2), rep("larva", 2)),
                         levels = c("oocyte", "cleavage", "gastrula", "introvertula", "larva"))
  
#  svg(paste0("~/Desktop/Figures/", paste0(crosref$GeneName_edit[crosref$GeneID == targets[i]], paste0("_pcau", ".svg"))))
  #png(paste0("~/Desktop/Figures/", paste0(crosref$GeneName_edit[crosref$GeneID == targets[i]], paste0("_pcau", ".png"))),
  #  width=550, height=550)
  print(ggplot(expr, aes(x = staging, y = expr, color = replicate, group = replicate)) +
    geom_point(size = 6) +
    ggtitle(crosref$GeneName_edit[crosref$GeneID == targets[i]]) +
    theme_bw() +
    scale_color_manual(values=wes_palettes$Rushmore[seq(1, length(wes_palettes), by = 2)]) +
    xlab("") +
    ylab("Regularized log transformed expression") +
    theme(text = element_text(size = 18)))
 # dev.off()
  
}
