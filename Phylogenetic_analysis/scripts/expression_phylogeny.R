##### LIBRARIES ######
library(phangorn)
library(tidyverse)

##### READ IN DATA #####
table1 <- read.table("~/Downloads/binary_OG_DF_norm-final.tsv") #expression values transformed to binary data (expression or no expression)
table2 <- read.table("~/Downloads/discr_OG_DF_norm-final.tsv") #expression values transformed to 5 category data based on cutoff values
categories <- read.table("~/Documents/Gene_expr_evol/Models/OG_categories.tsv")

##### DISCRETE DATA TABLE ######
#Discretized table formatting and exporting as alignment
table2[is.na(table2)] <- "?" #Missing values formatting
asd <- phyDat(table2, type = "USER", levels = c(0:9, "?"))
write.phyDat(asd, "~/Documents/Gene_expr_evol/Models/discr_mat.phy", format = "phylip")

##### BINARY DATA TABLES #####
#Data tables need pruning before use, doing that in a loop in this section
#Initialize variables
table_res_mat <- data.frame()
table_res_deg <- data.frame()

#Looping through orthogroups
for(i in 1:dim(table1)[1]){
  print(paste0("Working on maternal: ", i))
  
  #Maternal gene pruning
  temp_cat <- categories[categories$Category == "M", ]
  temp <- table1[i,][rownames(table1[i,]) %in% temp_cat$ID]
  
  #Adding binary data if present or not in category
  temp[names(temp) %in% temp_cat$species[temp_cat$ID %in% rownames(temp)]] <- 1
  temp[!(names(temp) %in% temp_cat$species[temp_cat$ID %in% rownames(temp)])] <- 0
  
  #Assigning to initialized variable
  table_res_mat <- rbind(table_res_mat, temp)
  
  print(paste0("Working on degraded: ", i))
  
  #Degraded gene pruning
  temp_cat <- categories[categories$Category == "D", ]
  temp <- table1[i,][rownames(table1[i,]) %in% temp_cat$ID]
  
  #Adding binary data if present or not in category
  temp[names(temp) %in% temp_cat$species[temp_cat$ID %in% rownames(temp)]] <- 1
  temp[!(names(temp) %in% temp_cat$species[temp_cat$ID %in% rownames(temp)])] <- 0
  
  #Assigning to initialized variable
  table_res_deg <- rbind(table_res_deg, temp)
  
}

#Transforming tables to alignment structure and exporting
table_res_deg <- phyDat(table_res_deg, type = "USER", levels = c(0:9, "?"))
write.phyDat(table_res_deg, "~/Documents/Gene_expr_evol/Models/binary_deg.phy", format = "phylip")

#Transforming tables to alignment structure and exporting
table_res_mat <- phyDat(table_res_mat, type = "USER", levels = c(0:9, "?"))
write.phyDat(table_res_mat, "~/Documents/Gene_expr_evol/Models/binary_mat.phy", format = "phylip")

##### PLOT OF IQTREE2 #####
#This section was used after running the IQ-TREE2 on each alignment

#Read in generated trees and species tree as reference point
bin_deg <- read.tree("~/Documents/Gene_expr_evol/Models/binary_deg.phy.contree")
bin_mat <- read.tree("~/Documents/Gene_expr_evol/Models/binary_mat.phy.contree")
dis_mat <- read.tree("~/Documents/Gene_expr_evol/Models/discr_mat.phy.contree")
spe_tre <- read.tree("~/tree_calibration/congr_sp_tree.dated.tre")
spe_tre$tip.label[spe_tre$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

#Create association matrixes used for cophylogenetic plots
assoc_deg <- matrix(rep(bin_deg$tip.label, 2), ncol = 2, nrow = length(bin_deg$tip.label))
assoc_mat <- matrix(rep(bin_mat$tip.label, 2), ncol = 2, nrow = length(bin_deg$tip.label))
assoc_mat_dis <- matrix(rep(dis_mat$tip.label, 2), ncol = 2, nrow = length(bin_deg$tip.label))
assoc1 <- cophylo(bin_deg, spe_tre, assoc_deg, rotate = T)
assoc2 <- cophylo(bin_mat, spe_tre, assoc_mat, rotate = T)
assoc3 <- cophylo(dis_mat, spe_tre, assoc_mat_dis, rotate = T)

#RF distances
treedist(spe_tre, bin_deg)
treedist(spe_tre, bin_mat)
treedist(spe_tre, dis_mat)

wRF.dist(spe_tre, dis_mat)
wRF.dist(spe_tre, bin_mat)
wRF.dist(spe_tre, bin_deg)

#Plotting co-phylo
tiff("~/Desktop/Figures/sp-tree_bin-deg_cophylo.tiff", 
     units="in", width=8.5, height=6.5, res=300)
svg("~/Desktop/Figures/sp-tree_bin-deg_cophylo.svg")
plot(assoc1,assoc.type="curved",link.lwd=1,link.lty="solid",
     link.col=make.transparent("blue",0.1),fsize=0.8)
nodelabels.cophylo(node=as.numeric(bin_deg$node.label[which(as.numeric(bin_deg$node.label) < 80)]),
                   which = "left", frame = "circle", bg = "firebrick", cex = 0.5)
dev.off()

tiff("~/Desktop/Figures/sp-tree_bin-mat_cophylo.tiff", 
     units="in", width=8.5, height=6.5, res=300)
svg("~/Desktop/Figures/sp-tree_bin-mat_cophylo.svg")
plot(assoc2,assoc.type="curved",link.lwd=1,link.lty="solid",
     link.col=make.transparent("blue",0.1),fsize=0.8)
nodelabels.cophylo(node=as.numeric(bin_mat$node.label[which(as.numeric(bin_mat$node.label) < 80)]),
                   which = "left", frame = "circle", bg = "firebrick", cex = 0.5)
dev.off()

tiff("~/Desktop/Figures/sp-tree_dis-mat_cophylo.tiff", 
     units="in", width=8.5, height=6.5, res=300)
svg("~/Desktop/Figures/sp-tree_dis-mat_cophylo.svg")
plot(assoc3,assoc.type="curved",link.lwd=1,link.lty="solid",
     link.col=make.transparent("blue",0.1),fsize=0.8)
nodelabels.cophylo(node=as.numeric(dis_mat$node.label[which(as.numeric(dis_mat$node.label) < 80)]),
                   which = "left", frame = "circle", bg = "firebrick", cex = 0.5)
dev.off()


#Plotting trees
svg("~/Desktop/Figures/sp-tree_bin-deg.svg")
print(ggtree(bin_deg) +
  geom_tiplab(size = 3) +                
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 80),
              size = 5, 
              shape = 21,
              fill = "forestgreen") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 80 & as.numeric(label) >= 50),
              size = 5, 
              shape = 21,
              fill = "gold") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 50),
              size = 5, 
              shape = 21,
              fill = "firebrick2") +
  theme(text = element_text(size = 15)) +
  theme_tree2() +
  xlim(c(0, .5)))
dev.off()


svg("~/Desktop/Figures/sp-tree_bin-mat.svg")
print(ggtree(bin_mat) +
  geom_tiplab(size = 3) +                
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 80),
              size = 5, 
              shape = 21,
              fill = "forestgreen") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 80 & as.numeric(label) >= 50),
              size = 5, 
              shape = 21,
              fill = "gold") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 50),
              size = 5, 
              shape = 21,
              fill = "firebrick2") +
    theme_tree2() +
  theme(text = element_text(size = 15)) +
  xlim(c(0, .6)))
dev.off()

svg("~/Desktop/Figures/sp-tree_dis-mat.svg")
print(ggtree(dis_mat) +
  geom_tiplab(size = 3) +                
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 80),
              size = 5, 
              shape = 21,
              fill = "forestgreen") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 80 & as.numeric(label) >= 50),
              size = 5, 
              shape = 21,
              fill = "gold") +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 50),
              size = 5, 
              shape = 21,
              fill = "firebrick2") +
    theme_tree2() +
  theme(text = element_text(size = 15)) +
  xlim(c(0, 3)))
dev.off()
       
print(ggtree(spe_tre) +
        geom_tiplab(size = 3) +                
        geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 80),
                    size = 5, 
                    shape = 21,
                    fill = "forestgreen") +
        geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 80 & as.numeric(label) >= 50),
                    size = 5, 
                    shape = 21,
                    fill = "gold") +
        geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 50),
                    size = 5, 
                    shape = 21,
                    fill = "firebrick2") +
        theme(text = element_text(size = 15)) +
        theme_tree2() +
        xlim(0, 1100))

      