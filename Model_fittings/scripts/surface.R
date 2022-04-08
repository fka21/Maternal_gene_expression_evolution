#Performing gene expression evoltuion analysis using the EVE model
#Ferenc Kagan
#19.12.2020
#Libraries
library(evemodel)
library(tidyverse)
library(phytools)
library(wesanderson)
library(geiger)
library(surface)

##### FUNCTION #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")
##### SAVE AND LOAD #####

#Read in expression table already pruned for maternal genes without degradation
df <- readRDS("~/Documents/Gene_expr_evol/Models/Results/mat_subset_expr_tbl.RDS")
#Read in ortohoup result table filtered for OU models and support > 0.5
ous <- read.table("~/Documents/Gene_expr_evol/Models/ou_models.tsv", sep = "\t", header = T, quote = "", row.names = 1)

#Read in species tree and adjust species name
species_tree <- read.tree("~/tree_calibration/congr_sp_tree.dated.tre")
species_tree$tip.label[species_tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

##### Building matrix and grouping variable #####
#Subset expression table for OU models
exprsMat <- as.matrix(subset(df, rownames(df) %in% rownames(ous)[ous$Category == "Maternal"]))

#Collapsing replicates
colExprsMat <- collapseReplicate(exprsMat)
colExprsMat[is.nan(colExprsMat)] <- NA

##### SRUFACE #####

#Initialize variables
surface_res <- vector(mode = "list", length =dim(colExprsMat)[1])
path_pic <- "~/Documents/Gene_expr_evol/Models/surface/"

#Looping through orthogroups
for(i in 1:length(t(colExprsMat))){
  print(i)
  #Preparing data
  sortedDat <- motmot::sortTraitData(phy = species_tree, 
                                     y = t(colExprsMat[i, , drop = F]), 
                                     log.trait = F,
                                     pass.ultrametric = T)
  
  #Naming missing node names
  test_tree <- nameNodes(sortedDat$phy)
  
  #Fitting surface
  test <- runSurface(test_tree, as.data.frame(sortedDat$trait))
  test_sum <- surfaceSummary(fwd = test$fwd, bwd = test$bwd)
  
  #Fitting single regime OU model
  test_ou <- startingModel(convertTreeData(sortedDat$phy, as.data.frame(sortedDat$trait))[[1]],
                           convertTreeData(sortedDat$phy, as.data.frame(sortedDat$trait))[[2]])
  
  #Calculating probability of single OU to explain data
  #Will use alpha of 0.05 for significance
  prob <- exp((test$bwd[[length(test$bwd)]]$aic - test_ou[[1]]$aic)/2)
  
  #If multiple regimes are better fit then save results
  if(prob <= 0.05){
    surface_res[[i]] <- list(test_sum, prob, colnames(t(colExprsMat[i, , drop = F])))
    
    #Plot result
    svg(paste0(path_pic, paste0(colnames(t(colExprsMat[i, , drop = F])), "_tree.svg")), height = 8)
    print(surfaceTreePlot(test_tree, test$bwd[[length(test$bwd)]], labelshifts = T,
                          edge.width = 2, label.offset = 4))
    dev.off()
    
    #Plot result
    svg(paste0(path_pic, paste0(colnames(t(colExprsMat[i, , drop = F])), ".svg")))
    print(surfaceTraitPlot(as.data.frame(sortedDat$trait), test$bwd[[length(test$bwd)]], 
                           convcol = T))
    dev.off()
    
  } else {
    next
  }
}

#Filter our empty list values
names(surface_res) <- seq_along(surface_res)
surface_res <- Filter(Negate(is.null), surface_res)

#Line used to inspect individual orthogroups
surface_res[sapply(surface_res, function(x) x[[3]]  == "OG0001726")]

#Codes used for extracting ditributions of regimes
asd <- data.frame(X1 = c(sapply(surface_res, function(x) x[[1]]$n_regimes[2]), 
                  rep(1, (3513 - length(surface_res))))) 

#Plot regimes
ggplot(asd, aes(x = as.character(X1)))+
  geom_bar(width = 0.8, fill = "black", alpha = 0.85) +
  theme_bw() +
  xlab("Regimes") +
  ylab("Frequency") +
  theme(text = element_text(size = 18))
#ggsave("~/Desktop/Figures/surface_#_regimes.tiff")

asd <- data.frame(X1 = sapply(surface_res, function(x) x[[1]]$n_regimes[1])) 
ggplot(asd, aes(x = X1))+
  geom_bar(width = 0.8, fill = "black", alpha = 0.85) +
  theme_bw() +
  xlab("Regime shifts") +
  ylab("Frequency") +
  theme(text = element_text(size = 18))
#ggsave("~/Desktop/Figures/surface_#_regime_shifts.svg")

#save.image('~/Documents/Gene_expr_evol/Intermediate_files/surface.RData')
#load("~/Documents/Gene_expr_evol/Intermediate_files/surface.RData")


##### NESTED MODEL #####
#These lines were used to inspect the presence of a possible nested model when inspecting results from SURFACE
library(motmot)

#Extract orthogroup with potential nested model and inspect
og <- colExprsMat[rownames(colExprsMat) == "OG0002894", ]
sortedDat <- sortTraitData(species_tree, y = as.data.frame(og), pass.ultrametric = T, log.trait = F)
plot(sortedDat$phy)
nodelabels() #Check which node the shift might have occured

ou.model <- transformPhylo.ML(sortedDat$trait, phy = sortedDat$phy, 
                              model = "OU") #using OU as comparison as that has been the model with highest support by default
nested.model <- transformPhylo.ML(sortedDat$trait, phy = sortedDat$phy, 
                                 model = "OU", nodeIDs = 28) #Nested model by default considers Brownian model and from given node the given model

#Chi-squared test
1 - pchisq(nested.model$MaximumLikelihood - ou.model$MaximumLikelihood, 
           1)
