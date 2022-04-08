###### LIBRARIES #####
library(ggstatsplot)
library(Biostrings)
library(rtracklayer)
library(ggpubr)
library(GenomicFeatures)
library(rstatix)
library(ggbiplot)
library(tidyverse)
library(sjPlot)


##### READ IN #####
#Excludingsome species which do not contain informations of UTRs
#Using GTF annotations to extract feature lengths
files <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Datasets/Gene_models", pattern = "*\\.g.f*", full.names = T)

#Need crossreference table for gene  and transcript IDs
tx2genes_paths <- list.files(path = "/Users/ferenckagan/Documents/Bioinformatic_analysis/Datasets/Gene_models", pattern = "*tsv", full.names = T)

#Read in transcript to gene crossreferences in a loop
tx2genes <- vector(mode = "list", length = length(files))
for(i in 1:length(tx2genes_paths)){
  tx2genes[[i]] <- read.table(tx2genes_paths[i], header = F)
  names(tx2genes)[i] <- str_remove_all(basename(tx2genes_paths)[i], "_tx2gene\\.tsv$")
}

#Adjusting IDs
for(i in 1:length(tx2genes)){
  if(names(tx2genes)[i] %in% c("caenorhabditis_elegans", "priapulus_caudatus")){
    temp <- tx2genes[[i]]
    tx2genes[[i]] <- temp
  } else{
  temp <- tx2genes[[i]]
  temp[,1] <- str_remove_all(temp[,1], "\\.[0-9]+$")
  tx2genes[[i]] <- temp}
  }
            

#Import maternal and degraded IDs
maternal_cleared_dge <- readRDS( "/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/downregulated_IDs_NFE-F.RDS")
maternal_broad <- readRDS("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_expr_evolution/Quantification/Intermediate_files/maternal_IDs.RDS")

###### EXTRACT FEATURES ######

#Initizialize variables
df_lengths_mat_dge <- vector(mode = "list", length = length(files))
df_lengths_mat_tpm <- vector(mode = "list", length = length(files))
df_lengths_zyg <- vector(mode = "list", length = length(files))
df_final <- vector(mode = "list", length = length(files))
df_lengths <- vector(mode = "list")

#Looping through each GTF file
for(i in 1:length(files)){
  format_type <- str_extract(files[i], "g.f")
  
  #Build TxDb object with corresponding transcript to gene crosreference table
  txdb <- makeTxDbFromGFF(file = files[i], format = format_type)
  tx2gene <- tx2genes[names(tx2genes) %in% str_to_lower(str_extract_all(basename(files[i]), "^[A-z]+_[a-z]+"))]
  
  
  #Extracting information from gtf annotation files for features
  threeUTRs <- threeUTRsByTranscript(txdb, use.names=TRUE)
  
  #3`UTRs are a filtering criteria
  if(is_empty(threeUTRs)){next}
  
  cdslength <- cdsBy(txdb, use.names = TRUE)
  length_threeUTRs   <- width(ranges(threeUTRs))
  length_cds <- width(ranges(cdslength))
  length_fiveUTRs <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  length_fiveUTRs   <- width(ranges(length_fiveUTRs))
  introns <- intronsByTranscript(txdb, use.names = TRUE)
  length_introns <- width(ranges(introns))
  
  #Extract lengths and calculate means or sums or them
  three_utr_lengths       <- as.data.frame(length_threeUTRs)
  three_utr_lengths        <- three_utr_lengths %>% group_by(group, group_name) %>% dplyr::summarise(round(mean(value)))
  three_utr_lengths        <- unique(three_utr_lengths[,c("group_name", "round(mean(value))")])
  colnames(three_utr_lengths) <- c("Transcript", "Mean 3' UTR Length")
  cds_lengths        <- as.data.frame(length_cds)
  cds_lengths        <- cds_lengths %>% group_by(group, group_name) %>% dplyr::mutate(length = sum(value), count = n())
  cds_lengths        <- unique(cds_lengths[,c("group_name", "length", "count")])
  colnames(cds_lengths) <- c("Transcript", "Total CDS Length", "Exon Count")
  five_utr_lengths        <- as.data.frame(length_fiveUTRs)
  five_utr_lengths        <- five_utr_lengths %>% group_by(group, group_name) %>% dplyr::summarise(round(mean(value)))
  five_utr_lengths        <- unique(five_utr_lengths[,c("group_name", "round(mean(value))")])
  colnames(five_utr_lengths) <- c("Transcript", "Mean 5' UTR Length")
  intron_lengths        <- as.data.frame(length_introns)
  intron_lengths        <- intron_lengths %>% group_by(group, group_name) %>% dplyr::mutate(length = sum(value), count = n())
  intron_lengths        <- unique(intron_lengths[,c("group_name", "length", "count")])
  colnames(intron_lengths) <- c("Transcript", "Total Intron Length", "Intron Count")
  cds_lengths_avg        <- as.data.frame(length_cds)
  cds_lengths_avg        <- cds_lengths_avg %>% group_by(group, group_name) %>% dplyr::summarise(round(mean(value)))
  cds_lengths_avg        <- unique(cds_lengths_avg[,c("group_name", "round(mean(value))")])
  colnames(cds_lengths_avg) <- c("Transcript", "Mean Exon Length")
  intron_lengths_avg       <- as.data.frame(length_introns)
  intron_lengths_avg        <- intron_lengths_avg %>% group_by(group, group_name) %>% dplyr::summarise(round(mean(value)))
  intron_lengths_avg        <- unique(intron_lengths_avg[,c("group_name", "round(mean(value))")])
  colnames(intron_lengths_avg) <- c("Transcript", "Mean Intron Length")
  
  #Merging lengths by transcripts 
  df_lengths[[i]] <- merge(cds_lengths, cds_lengths_avg, by = "Transcript") %>% merge(three_utr_lengths, by = "Transcript") %>%
    merge(five_utr_lengths, by = "Transcript") %>% 
    merge(intron_lengths, by = "Transcript") %>% merge(intron_lengths_avg, by = "Transcript") %>%
    relocate(`Exon Count`, `Intron Count`, .after = `Total Intron Length`)
  
  #Calculating proportion of given features compared to the sum of all known features and adding it to the table
  df_lengths[[i]]$Transcript <- tx2gene[[1]][match(df_lengths[[i]][,1], tx2gene[[1]][,1]), 2]
  df_lengths[[i]]$`Proportion CDS Length` <- df_lengths[[i]]$`Total CDS Length`/rowSums(df_lengths[[i]][,c(2,4,5,6)])
  df_lengths[[i]]$`Proportion Intron Length` <- df_lengths[[i]]$`Total Intron Length`/rowSums(df_lengths[[i]][,c(2,4,5,6)])
  df_lengths[[i]]$`Proportion 3' UTR Length` <- df_lengths[[i]]$`Mean 3' UTR Length`/rowSums(df_lengths[[i]][,c(2,4,5,6)])
  df_lengths[[i]]$`Proportion 5' UTR Length` <- df_lengths[[i]]$`Mean 5' UTR Length`/rowSums(df_lengths[[i]][,c(2,4,5,6)])
  
  #Subsetting for maternal genes
  df_lengths_mat_dge[[i]] <- subset(df_lengths[[i]], df_lengths[[i]][,1] %in% maternal_cleared_dge)
  df_lengths_mat_tpm[[i]] <- subset(df_lengths[[i]], df_lengths[[i]][,1] %in% maternal_broad)
  
  #Excluding overlaps of different categories
  df_lengths[[i]] <- subset(df_lengths[[i]], 
                              !(df_lengths[[i]][,1] %in% df_lengths_mat_dge[[i]][,1]) & 
                              !(df_lengths[[i]][,1] %in% df_lengths_mat_tpm[[i]][,1]))
  df_lengths_mat_tpm[[i]] <- subset(df_lengths_mat_tpm[[i]], 
                                      !(df_lengths_mat_tpm[[i]][,1] %in% df_lengths[[i]][,1]) & 
                                      !(df_lengths_mat_tpm[[i]][,1] %in% df_lengths_mat_dge[[i]][,1]))
  df_lengths_mat_dge[[i]] <- subset(df_lengths_mat_dge[[i]], 
                                      !(df_lengths_mat_dge[[i]][,1] %in% df_lengths[[i]][,1]) & 
                                      !(df_lengths_mat_dge[[i]][,1] %in% df_lengths_mat_tpm[[i]][,1]))

  
  #Addign conditions for later plot
  df_lengths[[i]]$'cond' <- "Reference"
  df_lengths_mat_dge[[i]]$'cond' <- "Degraded"
  df_lengths_mat_tpm[[i]]$'cond' <- "Persistently expressed"

  #Add results to final list of tables with species name attached
  df_final[[i]] <- rbind(df_lengths[[i]], df_lengths_mat_dge[[i]])
  df_final[[i]] <- rbind(df_final[[i]], df_lengths_mat_tpm[[i]])
  df_final[[i]]$'species' <- names(tx2gene)
  

}

#Removing enviorment variables to save up space
rm(list = ls(pattern = "tx2*"))
rm(list = ls(pattern = "*intron*"))
rm(list = ls(pattern = "*three*"))
rm(list = ls(pattern = "*five*"))
rm(list = ls(pattern = "*cds*"))

#Combine length data from list
df <- do.call(rbind.data.frame,df_final)
df <- reshape2::melt(df)

#Subset for absolute lengths
df_lengths_abs_res <- df %>% group_by(cond) %>% 
  filter(variable %in% c("Total CDS Length", "Total Intron Length", "Mean 3' UTR Length", "Mean 5' UTR Length", "Mean Exon Length", "Mean Intron Length")) %>%
  ungroup()
#Subset for feature occurence numbers
df_counts_res <- df %>% group_by(cond) %>% filter(variable == "Exon Count" | variable == "Intron Count") %>%
  ungroup() 
#Subset for proportions
df_lengths_prop_res <- df %>% group_by(cond) %>% 
  filter(variable %in% c("Proportion CDS Length", "Proportion Intron Length", "Proportion 3' UTR Length", "Proportion 5' UTR Length")) %>%
  ungroup() 

#Summarizing medians of all variables
df_summary <- df %>%
  select(cond,species,variable,value) %>% 
  group_by(cond, variable, species) %>% 
  dplyr::summarise(median = median(value)) %>% 
  pivot_wider(names_from = cond, values_from = median)
df_summary$species <- str_to_sentence(df_summary$species) #Adjusting to match tip labels downstream

#Export to a table
tab_df(df_summary,
       alternate.rows = T,
       file="~/Documents/sjt_test.doc")


##### PLOT WITH TREE #####

#Read in tree
tree <- read.tree('~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Orthofinder/congr_sp_tree.dated.tre')
species_tree$tip.label[species_tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

#Pruning species tree for species hacing feature annotations
tree <- keep.tip(tree, df_summary$species)

#Scaling tree height for better visualization
scale<-0.8*max(nodeHeights(tree))

#For every variable plot barplot with tree
for(i in 1:length(unique(df_summary$variable))){
  temp <- subset(df_summary, df_summary$variable %in% unique(df_summary$variable)[i])
  temp <- as.data.frame(temp)[, -1]
  rownames(temp) <- temp$species
  temp <- temp[, -1]
  colnames(temp) <- c("Downregulated", "Maternal", "Reference")
  
  plot_name <- str_replace_all(paste0(unique(df_summary$variable)[i], ".svg"), " ", "_")
  
  svg(file=paste("~/Documents/Bioinformatic_analysis/Gene_expr_evolution/Gene_lengths/Results", plot_name, sep = "/"))
  print(plotTree.barplot(tree, as.matrix(temp) , width = 2 , args.barplot=list(beside=TRUE, 
                                                            legend.text=TRUE, 
                                                            col=c("#003f5c", "#bc5090", "#31a354"),
                                                            xlab=unique(df_summary$variable)[i],
                                                            args.legend=list(x = 'top', bty = "n", ncol = 3, cex = 1.3,
                                                                             inset = c(0, -0.04)))))
  dev.off()
  
}

###### SAVE DATA ######
#save.image("/Users/ferenckagan/Documents/Bioinformatic_analysis/Quantification/Intermediate_files/stat_analysis.RData")
#load("/Users/ferenckagan/Documents/Bioinformatic_analysis/Quantification/Intermediate_files/stat_analysis.RData")

