#Script for estimating evolutionary models in gene expression data
#Ferenc.Kagan@uib.no
#23.02.2021

##### LIBRARYIES ######
library(motmot)
library(castor)
library(MuMIn)
library(tidyverse)
library(geiger)
library(phytools)

##### FUNCTION ######
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### READ TREE #####
#For reproducibility
set.seed(7)

tree <- read.tree('/home/s9/tsai/work/Ferenc/Models//congr_sp_tree.dated.tre')
OG_df <- read.table("/home/s9/tsai/work/Ferenc/Models/cont/batch_OG_DF_norm-final.tsv")
OG_cat <- read.table("/home/s9/tsai/work/Ferenc/Models/cont//OG_categories.tsv")

##### MODEL SELECTION #####
variables <- ls(pattern = "*_OG_df")

#Collapsing replicates into single values and calculating standard errors for each species measurement
OG_df_se <- t(SEcollapseReplicate(OG_df))
OG_df <- t(collapseReplicate(OG_df))

#Initialize result tables
res_df_M <- data.frame(OG = character(), 
                          Model = character(),
                          SigSQ  = double(),
                          Parameter = double(), 
                          Model_confidence = double(), 
                          Parameter_pval = double(),
                          Tree_size = double(),
			                    Mean_trait = double(),
			                    z0 = double())
res_df_D <- data.frame(OG = character(), 
                       Model = character(),
                       SigSQ  = double(),
                       Parameter = double(), 
                       Model_confidence = double(), 
                       Parameter_pval = double(),
                       Tree_size = double(),
		                   Mean_trait = double(),
		                   z0 = double())
categories <- c("M", "D")

#Looping through each orthogroup (OG)
for(k in 1:dim(OG_df)[2]){
  print(paste("Processing", k, sep = ": "))
  
  #Extract orthogroup with its standard errors
  temp <- OG_df[, k, drop = F]
  temp_se <- OG_df_se[, k, drop = F]
  
  #Looping through maternal gene categories to prune trees and fit models
  for(i in 1:length(categories)){
    #Subset categories table with current category and OG
    temp_cat <- categories[i]
    index <- OG_cat %>% tibble %>% filter(ID == colnames(temp) & Category == temp_cat) %>% unique()
    
    #Subsetting for OGs found in maternal gene category
    temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))
    temp_se_int <- as.vector(subset(temp_se, rownames(temp_se) %in% index$species))
    names(temp_se_int) <- rownames(subset(temp_se, rownames(temp_se) %in% index$species))
    
    #Some species have lost expression values, keeping these will artificially shift model fittings towards OU models
    #No standard errors are attributable to lost expression so using this as a subsetting criteria
    temp <- subset(temp, !(rownames(temp) %in% names(temp_se_int[temp_se_int == 0])))
    temp_se_int <- temp_se_int[temp_se_int != 0]
    
    #Pointless to use small trees as models will get misspecified
    #Using a tree size of 10 as cutoff
    if(length(temp[, 1]) < 10){next}

    #Preparing data
    sortedData <- sortTraitData(phy = tree, y = temp, log.trait = F, pass.ultrametric = T)
        
    print(paste0("Fitting models ", temp_cat))
    #Model fits
    bm.ml <- fitContinuous(sortedData$phy, sortedData$trait, SE = temp_se_int, model = "BM", ncore = 9)
    ou.ml <- fitContinuous(sortedData$phy, sortedData$trait, SE = temp_se_int, model = "OU", ncore = 9)
    white.ml <- fitContinuous(sortedData$phy, sortedData$trait, SE = temp_se_int, model = "white", ncore = 9)
    
    #Calculate AICc weights
    res_df <- data.frame(AICc = c(bm.ml$opt$aicc, ou.ml$opt$aicc, white.ml$opt$aicc))
    rel_weight <- Weights(res_df$AICc)
    names(rel_weight) <- c("Brownian", "OU", "White")
    
    #Which is the best fitting model
    model_max <- names(which.max(rel_weight))
    
    #Simulate based on winning model parameters
    if(model_max == "Brownian"){
      
      simtrait <- sim.rates(sortedData$phy, bm.ml$opt$sigsq, bm.ml$opt$z0, nsim=100, internal=FALSE, plot=FALSE)
      
    } else if (model_max == "OU") {
      
      simtrait <- transformPhylo.sim(sortedData$phy, sortedData$trait, model = "OU", alpha = ou.ml$opt$alpha, n = 100)
      
    }  else if (model_max == "White") {
      
      simtrait <- matrix(rnorm(length(sortedData$trait) * 100, white.ml$opt$z0), ncol = 100)
      rownames(simtrait) <- sortedData$phy$tip.label
      
    } 
    print(paste0("Fitting simulated data ", temp_cat))
    
    #Fit models to simulated data
    bm.sim <- apply(X = simtrait, 2, FUN = fitContinuous, phy = sortedData$phy, model = "BM", ncore = 9)
    ou.sim <- apply(X = simtrait, 2, FUN = fitContinuous, phy = sortedData$phy, model = "OU", ncore = 9)
    white.sim <- apply(X = simtrait, 2, FUN = fitContinuous, phy = sortedData$phy, model = "white", ncore = 9)
    
    #Initializing intermediate variables
    sim_param <- data.frame(Param = double())
    sim_max <- c()
    
    #Extracting fits and comparing to original model
    for(m in 1:100){
      sim_res <- data.frame(AICc = c(bm.sim[[m]]$opt$aicc, ou.sim[[m]]$opt$aicc, white.sim[[m]]$opt$aicc))
      sim_weight <- Weights(sim_res$AICc)
      names(sim_weight) <- c("Brownian", "OU", "White")
      sim_max <- append(names(which.max(sim_weight)), sim_max)
      
      #Using K-S test to test goodness-of-fit for estimated parameter and simulated parameters
      if(model_max == "Brownian"){
        
        sim_param[m, ] <- bm.sim[[m]]$opt$sigsq
        param_conf <- ks.test( bm.ml$opt$sigsq, sim_param$Param)
        
      } else if (model_max == "OU") {
        
        sim_param[m, ] <- ou.sim[[m]]$opt$alpha
        param_conf <- ks.test( ou.ml$opt$alpha, sim_param$Param)
        
      }  else if (model_max == "White") {
        
        sim_param[m, ] <- white.sim[[m]]$opt$sigsq
        param_conf <- ks.test(white.ml$opt$sigsq, sim_param$Param)
        
      } 
    }
    #What % is the model replicated
    model_conf <- sum(sim_max %in% model_max) / length(sim_max)
    
    #Save everything in output table
    if(model_max == "Brownian"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- bm.ml$opt$sigsq
        res_df_M[k, 4] <- NA
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
	res_df_M[k, 8] <- mean(sortedData$trait)
	res_df_M[k, 9] <- bm.ml$opt$z0
      } else if (i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- bm.ml$opt$sigsq
        res_df_D[k, 4] <- NA
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- param_conf$p.value
        res_df_D[k, 7] <- length(sortedData$phy$tip.label)
	res_df_D[k, 8] <- mean(sortedData$trait)
        res_df_D[k, 9] <- bm.ml$opt$z0}
    } else if(model_max == "OU"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- ou.ml$opt$sigsq
        res_df_M[k, 4] <- ou.ml$opt$alpha
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
	res_df_M[k, 8] <- mean(sortedData$trait)
        res_df_M[k, 9] <- ou.ml$opt$z0
      } else if(i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- ou.ml$opt$sigsq
        res_df_D[k, 4] <- ou.ml$opt$alpha
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- param_conf$p.value
        res_df_D[k, 7] <- length(sortedData$phy$tip.label)
	res_df_D[k, 8] <- mean(sortedData$trait)
        res_df_D[k, 9] <- ou.ml$opt$z0}
    } else if(model_max == "White"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- white.ml$opt$sigsq
        res_df_M[k, 4] <- NA
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
	res_df_M[k, 8] <- mean(sortedData$trait)
        res_df_M[k, 9] <- white.ml$opt$z0
      } else if(i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- white.ml$opt$sigsq
        res_df_D[k, 4] <- NA
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
	res_df_D[k, 8] <- mean(sortedData$trait)
        res_df_D[k, 9] <- white.ml$opt$z0}
    }
    
    
    print(paste("Finished", i, sep = " "))
    
  }

res_df_D <- res_df_D[!(is.na(res_df_D$OG)), ]
res_df_M <- res_df_M[!(is.na(res_df_M$OG)), ]


write.table(res_df_D, "/home/s9/tsai/work/Ferenc/Models/cont/fits_degr_out_1.tsv", quote = F, col.names = T, row.names = F)
write.table(res_df_M, "/home/s9/tsai/work/Ferenc/Models/cont/fits_pers_out_1.tsv", quote = F, col.names = T, row.names = F)

}

       
save.image("/home/s9/tsai/work/Ferenc/Models/cont//continuous_model_fit.RData")

