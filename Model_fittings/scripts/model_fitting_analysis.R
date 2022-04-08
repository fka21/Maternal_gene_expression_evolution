##### LIBRARYIES ######
library(motmot)
library(castor)
library(MuMIn)
library(tidyverse)
library(geiger)
library(UpSetR)
library(wesanderson)

##### FUNCTION #####
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

###### ANALYZE FITINGS #####
og_cats <- read.table("~/Documents/Gene_expr_evol/Models/OG_categories.tsv", header = T)
og_pres <- read.table("~/Documents/Gene_expr_evol/Models/OG_presence.tbl", header = T, sep = "\t")

#Annotation for orthogroups using a python scripts from:
#https://github.com/davidemms/OrthoFinder/issues/451
annot <- read.table("~/Documents/Gene_expr_evol/Models/N0_blasted.tsv", header = T, sep = "\t", quote = "")
annot <- tibble(annot) %>% dplyr::select(OG, Best_hit, E_val) %>% group_by(OG) %>% dplyr::slice(which.min(E_val))

#Read in model fittings for degraded maternal genes
degr_res <- read.table("~/Documents/Gene_expr_evol/Models/fits_degr_out.tsv", header = F, row.names = 1)
colnames(degr_res) <- degr_res[rownames(degr_res) == "OG", ]
degr_res$Category <- "Degraded"
degr_res$Annot <- annot$Best_hit[match(rownames(degr_res), annot$OG)]

#Read in model fittings for non-degraded genes
pers_res <- read.table("~/Documents/Gene_expr_evol/Models/fits_pers_out.tsv", header = F, row.names = 1)
colnames(pers_res) <- pers_res[rownames(pers_res) == "OG", ]
pers_res$Category <- "Maternal"
pers_res$Annot <- annot$Best_hit[match(rownames(pers_res), annot$OG)]

#Combine results and format
df <- rbind(pers_res, degr_res)
df <-df[!(rownames(df) %in% c("OG", "OG1")), ]
df$SigSQ <- as.numeric(df$SigSQ)
df$Parameter <- as.numeric(df$Parameter)
df$Model_confidence <- as.numeric(df$Model_confidence)
df$Parameter_pval <- as.numeric(df$Parameter_pval)
df$Tree_size <- as.numeric(df$Tree_size)
df$Mean_trait <- as.numeric(df$Mean_trait)
df$z0 <- as.numeric(df$z0)

#Export for later use
write.table(df[df$Model == "OU" & df$Model_confidence > 0.5, ], "~/Documents/Gene_expr_evol/Models/ou_models.tsv", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(df[df$Model == "Brownian" & df$Model_confidence > 0.5, ], "~/Documents/Gene_expr_evol/Models/bm_models.tsv", sep = "\t", quote = F, col.names = T, row.names = T)

##### ORTHOLOGY RELATIONSHIPS #######
#This section was sued to inspect orthology inference results

#Preparing categories table
og_cats <- og_cats %>% tibble() %>% filter(!(is.na(Category)))
og_cats$species <- factor(og_cats$species, levels = species_tree$tip.label)

#Initialize variables
M <- vector(mode = "list", length = 43)
D <- vector(mode = "list", length = 43)

#Looping through each species and extracting orthogroup names belonging to category
for(i in 1:length(unique(og_cats$species))){
  temp_sp <- unique(og_cats$species)[i]
  temp <- filter(og_cats, species == temp_sp)
  
  M[[i]] <- temp$ID[temp$Category == "M"]
  names(M)[i] <- temp_sp
  D[[i]] <- temp$ID[temp$Category == "D"]
  names(D)[i] <- temp_sp
}

#Combining for overlap across species
M <- fromList(M)
D <- fromList(D)

#Prepare data for plotting
plot_df <- data.frame(X = c(rowSums(M), rowSums(D), rowSums(og_pres)),
                      Y = c(rep("Maternal", dim(M)[1]), rep("Degraded", dim(D)[1]), rep("Reference", dim(og_pres)[1])))

#Plot OG presence in a per category basis
svg("~/Desktop/Figures/Overall_OG_presence.svg")
print(ggplot(plot_df, aes(x = X, fill = Y)) +
  geom_bar(alpha = 0.9) +
  facet_grid(rows = vars(Y)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1000)) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$GrandBudapest1) +
  xlab(NULL) + 
  ylab(NULL))
dev.off()

#Plot for OG presence in a per species basis
svg("~/Desktop/Figures/Species_OG_presence.svg")
print(og_cats %>% 
  group_by(Category, species) %>% 
  dplyr::summarise(n = n()) %>% 
  mutate(category = case_when(Category == "D" ~ "Degraded",
                              Category == "M" ~ "Maternal")) %>%
  ggplot(aes(x = species, y = n, fill = category)) + 
  geom_bar(stat = "identity", alpha = 0.9) + 
  facet_grid(rows = vars(category), scales = "free") + 
  theme_bw() + 
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
        legend.position = 'none') +
  scale_fill_manual(values=wes_palettes$GrandBudapest1) +
  xlab(NULL) + 
  ylab(NULL))
dev.off()

##### BM rate heterogeneous #####
#This section was used to fit heterogenus rate Brownian model fits

#Read in data
exprsMat <- read.table("~/Documents/Gene_expr_evol/Models/batch_OG_DF_norm-final.tsv", sep = "\t", header = T, row.names = 1)
geneCat <- read.table("~/Documents/Gene_expr_evol/Models/OG_categories.tsv", header = T, sep = " ", row.names = 1)
species_tree <- read.tree("~/tree_calibration/congr_sp_tree.dated.tre")
species_tree$tip.label[species_tree$tip.label == "Strongylocentrotus_franciscanus"] <- "Mesocentrotus_franciscanus"

#Collapse replicates before fittings
exprsMatSE <- SEcollapseReplicate(exprsMat)
exprsMat <- collapseReplicate(exprsMat)
exprsMat[is.nan(exprsMat)] <- NA

#Extract orthogroup names where Brownian model was winning model
bms_degr <- rownames(df)[df$Model == "Brownian" & df$Category  == "Degraded" & df$Model_confidence > 0.5] 
bms_mate <- rownames(df)[df$Model == "Brownian" & df$Category  == "Maternal" & df$Model_confidence > 0.5] 

#Initialize variables
degr_out <- vector(mode = "list", length = dim(exprsMat)[1])
mate_out <- vector(mode = "list", length =dim(exprsMat)[1])

#Loop through orthogroups for degraded genes
for(k in 1:(dim(exprsMat)[1])){
  #Subsetting OGs
  print(paste0("Processing: ", k))
  temp <- t(exprsMat)[, k, drop = F]
  temp_se <- t(exprsMatSE)[, k, drop = F]

  #Retrieving indexes for gene categories
  index <- geneCat %>% tibble %>% filter(ID == colnames(temp) & Category == "D") %>% unique()
  
  #Subsetting for OGs found in maternal gene category
  temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))
  temp_se_int <- as.vector(subset(temp_se, rownames(temp_se) %in% index$species))
  names(temp_se_int) <- rownames(subset(temp_se, rownames(temp_se) %in% index$species))
    
  #Some species have lost expression values, keeping these will artificially shift model fittings towards OU models
  #No standard errors are attributable to lost expression so using this as a subsetting criteria
  temp <- subset(temp, !(rownames(temp) %in% names(temp_se_int[temp_se_int == 0])))
  temp_se_int <- temp_se_int[temp_se_int != 0]
  
  #Filtering for BM OGs only with tree filtering
  if(!(colnames(temp) %in% bms_degr) | length(rownames(temp)) < 10){next}
  print(paste0("Passed filter: ", k))
  
  #Sorting data
  sortData<- motmot::sortTraitData(species_tree, y = temp, pass.ultrametric = T,
                                log.trait = F)
  
  print(paste0("Fitting: ", k))
  #Fitting heterogeneous rate model
  tm2.ml <- transformPhylo.ML(y = sortData$trait, phy = sortData$phy, 
                              model = "tm2", minCladeSize = 2, nSplits = 4,
                              meserr =  temp_se_int,
                              n.cores =  3,
                              modelCIs = T)
  
  #Summarize model and save to result lise
  sumModel <- motmot::summary.traitMedusa(tm2.ml)
  degr_out[[k]] <- list(colnames(temp), sumModel$ModelFit, sumModel$Rates, sumModel$optimalTree, sumModel$original.phy)
  
  
}

saveRDS(degr_out, "~/Documents/Gene_expr_evol/Models/Results/heterogeneous_model_degr.RDS")

#Loop through orthogroups for non-degraded genes
for(k in 1:(dim(exprsMat)[1])){
  print(paste0("Processing: ", k))
  #Subsetting OGs
  temp <- t(exprsMat)[, k, drop = F]
  temp_se <- t(exprsMatSE)[, k, drop = F]
  
  #Retrieving indexes for gene categories
  index <- geneCat %>% tibble %>% filter(ID == colnames(temp) & Category == "M") %>% unique()
  
  #Subsetting for OGs found in maternal gene category
  temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))
  temp_se_int <- as.vector(subset(temp_se, rownames(temp_se) %in% index$species))
  names(temp_se_int) <- rownames(subset(temp_se, rownames(temp_se) %in% index$species))
  
  #Some species have lost expression values, keeping these will artificially shift model fittings towards OU models
  #No standard errors are attributable to lost expression so using this as a subsetting criteria
  temp <- subset(temp, !(rownames(temp) %in% names(temp_se_int[temp_se_int == 0])))
  temp_se_int <- temp_se_int[temp_se_int != 0]
  
  #Filtering for BM OGs only  
  if(!(colnames(temp) %in% bms_mate) | length(rownames(temp)) < 10){next}
  print(paste0("Passed filter: ", k))
  
  #Sorting data
  sortData<- motmot::sortTraitData(species_tree, y = temp, pass.ultrametric = T,
                                   log.trait = F)
  
  print(paste0("Fitting: ", k))
  #Fitting heterogeneous rate model
  tm2.ml <- transformPhylo.ML(y = sortData$trait, phy = sortData$phy, 
                              model = "tm2", minCladeSize = 2, nSplits = 4,
                              meserr =  temp_se_int,
                              n.cores =  3,
                              modelCIs = T)
  
  #Summarize model and save to result lise
  sumModel <- motmot::summary.traitMedusa(tm2.ml)
  mate_out[[k]] <- list(colnames(temp), sumModel$ModelFit, sumModel$Rates, sumModel$optimalTree, sumModel$original.phy)
  
}

saveRDS(mate_out, "~/Documents/Gene_expr_evol/Models/Results/heterogeneous_model_mate.RDS")

#Filter empty list objects
names(degr_out) <- seq_along(degr_out)
degr_out <- Filter(Negate(is.null), degr_out)
names(mate_out) <- seq_along(mate_out)
mate_out <- Filter(Negate(is.null), mate_out)

#Filter out where no model fit improvement has been achieved by applying a heterogenous model fit
degr_out <- degr_out[sapply(degr_out, function(x) !(x[3] == "Single rate"))]
mate_out <- mate_out[sapply(mate_out, function(x) !(x[3] == "Single rate"))]

#Loop through resuls for plotting
for(i in 1:length(mate_out)){
    #Extract results object
    x <- mate_out[[i]]
    #Which edges belong to which rate
    temp <- x[[5]]$edge[, 1] %in% c(x[[3]][1], getDescendants(x[[5]], as.character(x[[3]][1])))
    temp[temp == T] <- "red"
    temp[temp == F] <- "black"
    
    name <- x[[1]]
    
    #Plot
    tiff(paste0("~/Documents/Gene_expr_evol/Models/Results/Heterogeneous_models/", paste0(name, ".tiff")),
         height = 850, width = 750)
    par(mfrow = c(2,1))
    plot(x[[5]], main = x[[1]], label.offset = 3, cex = 1, edge.width = 1.8)
    axisPhylo()
    plot(x[[4]], label.offset = 3, cex = 1, edge.width = 1.8, edge.color = temp)
    axisPhylo()
    
    dev.off()
    
    svg(paste0("~/Documents/Gene_expr_evol/Models/Results/Heterogeneous_models/", paste0(name, ".svg")),
               height = 14, width = 8)
    par(mfrow = c(2,1))
    plot(x[[5]], main = x[[1]], label.offset = 3, cex = 0.8, edge.width = 1.8)
    axisPhylo()
    plot(x[[4]], label.offset = 3, cex = 0.8, edge.width = 1.8, edge.color = temp)
    axisPhylo()
    
    dev.off()
}

#How many sigma^2 rates can be found across the maternal brownian orthogroups?
svg("~/Desktop/Figures/heterogeneous_rates.svg")
print(c((sapply(mate_out, function(x){ dim(x[[3]])[1]}) + 1), rep(1, 549))
      %>% tibble() %>% ggplot(aes(x = as.character(.))) + 
  geom_bar(width = 0.5, fill = "black", alpha = 0.9) +
  theme_bw())
dev.off()

##### HETEROGENEOUS NESTED MODEL #####
#This chunk of code can be used to test if a nested model improves AICc scores or
#alternatively results in a significant increase in logL

library(motmot)

#Extract orthogroup with potential nested model and inspect
temp <- exprsMat[rownames(exprsMat) == "OG0000715", ]
asd <- t(collapseReplicate(temp))
asd_se <- t(SEcollapseReplicate(temp))
asd_se <- asd_se[!(is.na(asd_se[,1])), , drop = F]
asd_se <- asd_se[!(asd_se[,1] == 0), , drop = F]
asd <- asd[rownames(asd) %in% rownames(asd_se), , drop = F]
sortedDat <- motmot::sortTraitData(species_tree, y = asd[,1, drop = F], log.trait = F,
                                   pass.ultrametric = T)
#Inspecting which node is needed
plot(sortedDat$phy)
nodelabels()

#Fitting models
bm.model <- transformPhylo.ML(sortedDat$trait, phy = sortedDat$phy, 
                              model = "bm", n.cores = 4, meserr = asd_se)  #using Brownina as comparison as that has been the model with highest support by default
tm2.ml <- transformPhylo.ML(y = sortedDat$trait, phy = sortedDat$phy, 
                            model = "tm2", minCladeSize = 2, nSplits = 4, n.cores = 4, meserr = asd_se)
sumModel <- motmot::summary.traitMedusa(tm2.ml) #using heterogenous model also as comparison as it improved model fit
nested.model <- transformPhylo.ML(sortedDat$trait, phy = sortedDat$phy, 
                                  model = "OU", nodeIDs = 46, n.cores = 4, meserr = asd_se)

#Checking the lowest AICc scores
c(bm.model$AICc, sumModel$ModelFit[,4], nested.model$AICc)

1 - pchisq(nested.model$MaximumLikelihood - bm.model$logLikelihood, 
           1)
1 - pchisq(nested.model$MaximumLikelihood - sumModel$ModelFit[,1], 
           1)


##### CONTINUOUS MODELS ANALYSIS #####
library(wesanderson)

#Plots of results model confidence distributions
svg("~/Desktop/Figures/Model_confidence.svg")
png("~/Desktop/Figures/Model_confidence.png", width = 750, height = 550)
print(ggplot(df, aes(x = Model_confidence, fill = Model)) +
  geom_density(alpha = 0.9) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$Royal1) +
  facet_grid(rows = vars(Category)))
dev.off()

#Keeping only support > 0.5 for belowe
df <- df[df$Model_confidence > 0.5, ]
dim(df)

#Plot proportion of winning models
svg("~/Desktop/Figures/Winning_model.svg")
png("~/Desktop/Figures/Winning_model.png", width = 750, height = 550)
print(ggplot(df, aes(x = Category, fill = Model)) +
  geom_bar(position = 'fill', width = 0.5, color = "black", alpha = 0.9) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$Royal1) +
  xlab(NULL) +
  ylab(NULL))
dev.off()

#Plot p-value distributions for each model and category
svg("~/Desktop/Figures/Winning_model_pval.svg")
png("~/Desktop/Figures/Winning_model_pval.png", width = 750, height = 550)
print(ggplot(df, aes(x = Parameter_pval, fill = Model)) +
  geom_density(alpha = 0.8) +
  theme_bw() +
  xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$Royal1) +
  facet_grid(rows = vars(Category)))
dev.off()

#Plot tree sizes for each model and category
svg("~/Desktop/Figures/Winning_model_treesize.svg")
png("~/Desktop/Figures/Winning_model_treesize.png", width = 750, height = 550)
print(ggplot(df, aes(x = Model, y = Tree_size, fill = Model)) +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  xlab(NULL) +
  ylab(NULL) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$Royal1) +
  facet_grid(rows = vars(Category)))
dev.off()

#Subset intermediate tables for parameter comparisons
df_param_brown <- subset(df, df$Model == "Brownian")
df_param_ou <- subset(df, df$Model == "OU")

#Plot Brownian sigsq for each  category
svg("~/Desktop/Figures/Winning_model_sigsq_brown.svg")
png("~/Desktop/Figures/Winning_model_sigsq_brown.png", width = 750, height = 550)
print(ggplot(df_param_brown, aes(x = log(SigSQ), fill = Category)) +
  geom_density(alpha = 0.7) +
  theme_bw() +
  ylab(NULL) +
  theme(text = element_text(size = 15), legend.position = 'none') +
  scale_fill_manual(values=wes_palettes$GrandBudapest1) +
  xlab(expression(log(sigma^2))))
dev.off()

#Plot OU alpha for each category
svg("~/Desktop/Figures/Winning_model_alpha_ou.svg")
png("~/Desktop/Figures/Winning_model_alpha_ou.png", width = 750, height = 550)
print(ggplot(df_param_ou, aes(x = log(Parameter), fill = Category)) +
  geom_density(alpha = 0.7) +
  theme_bw() +
  ylab(NULL) +
  theme(text = element_text(size = 15), legend.position = 'none') +
  scale_fill_manual(values=wes_palettes$GrandBudapest1) +
  xlim(-7, -4) +
  xlab(expression(log(alpha))))
dev.off()

#Plot OU sigsq for each  category
svg("~/Desktop/Figures/Winning_model_sigsq_ou.svg")
png("~/Desktop/Figures/Winning_model_sigsq_ou.png", width = 750, height = 550)
print(ggplot(df_param_ou, aes(x = SigSQ, fill = Category)) +
  geom_density(alpha = 0.7)  +
  theme_bw() +
  ylab(NULL) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$GrandBudapest1)) +
  xlab(expression(log(sigma^2)))
dev.off()

#Plot mean traits for each  category
svg("~/Desktop/Figures/Winning_model_mean_trait.svg")
png("~/Desktop/Figures/Winning_model_mean_trait.png", width = 750, height = 550)
print(ggplot(df, aes(x = Mean_trait, fill = Category)) +
        geom_density(alpha = 0.7)  +
        theme_bw() +
        ylab(NULL) +
        theme(text = element_text(size = 15), legend.position = 'none') +
        scale_fill_manual(values=wes_palettes$GrandBudapest1)) +
  xlab("log(TPM)")
dev.off()

##Testing distributions for mean traits
mean_trait <- wilcox.test(Mean_trait ~ Category, df)
mean_trait_brwn <- wilcox.test(Mean_trait ~ Category, df_param_brown)
mean_trait_ou <- wilcox.test(Mean_trait ~ Category, df_param_ou)

#Inspecting effect sizes
rstatix::wilcox_effsize(df, Mean_trait ~ Category)
rstatix::wilcox_effsize(df_param_brown, Mean_trait ~ Category)
rstatix::wilcox_effsize(df_param_ou, Mean_trait ~ Category)

##Testing distributions for other parameters with effect sizes
brown_sigsq <- wilcox.test(log(SigSQ) ~ Category, df_param_brown)
rstatix::wilcox_effsize(df_param_brown, SigSQ ~ Category)
ou_sigsq <- wilcox.test(log(SigSQ) ~ Category, df_param_ou)
rstatix::wilcox_effsize(df_param_ou, SigSQ ~ Category)
ou_alpha <- wilcox.test(log(Parameter) ~ Category, df_param_ou)

#Testing bimodality and with bimodality
test <- subset(df_param_ou, df_param_ou$Category == "Degraded")
multimode::modetest(log(test$SigSQ))

###### Mk model for losing expression #####
#Import results for binary trait model fittings
mk_mat <- read.table("~/Documents/Gene_expr_evol/Models/binary_fits_pers_out.tsv", header = T, sep = " ")
mk_mat$Category <- "Maternal"
mk_deg <- read.table("~/Documents/Gene_expr_evol/Models/binary_fits_degr_out.tsv", header = T, sep = " ")
mk_deg <- mk_deg[!(is.na(mk_deg$OG)), ]
mk_deg$Category <- "Degraded" 

#Combine data from the 2 categories
mk <- rbind(mk_mat, mk_deg)

#Subset for high confidence and tree size
mk <- mk[(mk$Model_conf > 0.5) & (mk$Tree_size >= 10), ]
mk 
mk <- mk[!(is.na(mk$Model)), ]

#Plot
svg("~/Desktop/Figures/mk_winninng_models.svg")
print(ggplot(mk, aes(x = Category, fill = Model)) +
  geom_bar(stat = "count", alpha = 0.8, position = 'dodge', width = 0.5) +
  #facet_wrap(~Category) +
  theme_bw() +
  xlab(NULL) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$BottleRocket1))
dev.off()

svg("~/Desktop/Figures/mk_transitions.svg")
print(ggplot(mk, aes(x = log(Q12), fill = Category)) +
  geom_density(alpha = 0.6) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$GrandBudapest1))
dev.off()

svg("~/Desktop/Figures/mk_model_conf.svg")
print(ggplot(mk, aes(x = Model_conf, fill = Model)) +
  geom_density(alpha = 0.8) +
  facet_wrap(~Category)+
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$BottleRocket1))
dev.off()

svg("~/Desktop/Figures/mk_tree_sizes.svg")
print(ggplot(mk, aes(y = Tree_size, x = Model ,fill = Model)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~Category) +
  theme_bw() +
  xlab(NULL) +
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values=wes_palettes$BottleRocket1))
dev.off()

