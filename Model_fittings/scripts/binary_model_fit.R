##### LIBRARIES ######
library(ape)
library(phytools)
library(tidyverse)
library(geiger)
library(MuMIn)
library(castor)
library(motmot)

##### FUNCTION ######
source("~/Documents/Gene_expr_evol/Scripts/functions.R")

##### READ IN DATA #####
tree <- read.tree('/home/s9/tsai/work/Ferenc/Models/congr_sp_tree.dated.tre')
OG_df <- read.table("/home/s9/tsai/work/Ferenc/Models/cont/batch_OG_DF_norm-final.tsv")
OG_cat <- read.table("/home/s9/tsai/work/Ferenc/Models/cont/OG_categories.tsv")

#Collapse replicates to single mean values
OG_df <- t(collapseReplicate(OG_df))

##### MK MODEL FITTING ######
#Initialize results tables
res_df_M <- data.frame(OG = character(),
                       Model = character(),
                       Q12 = numeric(),
                       Q21 = numeric(),
                       Model_conf = numeric(),
                       Parameter_pval = numeric(),
                       Tree_size = numeric())
res_df_D <- data.frame(OG = character(),
                       Model = character(),
                       Q12 = numeric(),
                       Q21 = numeric(),
                       Model_conf = numeric(),
                       Parameter_pval = numeric(),
                       Tree_size = numeric())
categories <- c("M", "D")

#Loop through each orthogroup
for(k in 1:dim(OG_df)[2]){
  print(paste("Processing", k, sep = ": "))
  
  #Extract orthogroup
  temp <- OG_df[, k, drop = F]
  
  #Looping through maternal gene categories to prune trees and fit models
  for(i in 1:length(categories)){
    #Subset categories table with current category and OG
    temp_cat <- categories[i]
    index <- OG_cat %>% tibble %>% filter(ID == colnames(temp) & Category == temp_cat) %>% unique()
    
    #Subsetting for OGs found in maternal gene category
    temp <- as.matrix(subset(temp, rownames(temp) %in% index$species))

    #Discretizing, if expression >= 2, then A (using log(2) as the whole data frame is in log scale), otherwise B
    temp[ temp >= log(2)] <- "A"
    temp[ temp < log(2)] <- "B"
    
    #Skipping homogeneous orthogroups as fitting is pointless
    if(!(("A" %in% temp) & ("B" %in% temp))){next}

    #Adjusting data and tree
    sortedData <- sortTraitData(phy = tree, y = temp, log.trait = F, pass.ultrametric = T)
    
    #Fitting Mk models
    ard.ml <- fitDiscrete(sortedData$phy, sortedData$trait, model = "ARD", ncores = 5)
    er.ml <- fitDiscrete(sortedData$phy, sortedData$trait, model = "ER", ncores = 5)
    white.ml <- fitDiscrete(sortedData$phy, sortedData$trait, model = "ER", transform = "white", ncores = 5)
    
    #Calculate AICc weights
    res_df <- data.frame(AICc = c(ard.ml$opt$aicc, er.ml$opt$aicc, white.ml$opt$aicc))
    rel_weight <- Weights(res_df$AICc)
    names(rel_weight) <- c("ARD", "ER", "ER-white")
    
    #Save winning model for subsampled tree
    model_max <- names(which.max(rel_weight))
    
    #Initialize simulation data.frame
    simtrait <- data.frame(matrix(ncol = 100, nrow = length(sortedData$trait)))
    
    #Simulate data based on winning model
   if(model_max == "ARD"){
        for(sim_rep in 1:100){
          simtrait[,sim_rep] <- rTraitDisc(sortedData$phy, model = "ARD", k = 2, rate = c(ard.ml$opt$q12, ard.ml$opt$q21))
          rownames(simtrait) <- names(rTraitDisc(sortedData$phy, model = "ARD", k = 2, rate = c(ard.ml$opt$q12, ard.ml$opt$q21)))
        }
        
      } else if (model_max == "ER") {
        
        for(sim_rep in 1:100){
          simtrait[,sim_rep] <- rTraitDisc(sortedData$phy, model = "ER", k = 2, rate = er.ml$opt$q12)
          rownames(simtrait) <- names(rTraitDisc(sortedData$phy, model = "ER", k = 2, rate = er.ml$opt$q12))
        }
        
      } else if (model_max == "ER-white"){
         
         #To get a null-hypothesis simulation transforming tree into star phylogeny
         #Following that simulate along this "tree" ER transitions
         temp_tree <- transformPhylo(sortedData$phy, "lambda", lambda = 0)
         
        for(sim_rep in 1:100){
          simtrait[,sim_rep] <- rTraitDisc(temp_tree, model = "ER", k = 2, rate = er.ml$opt$q12)
          rownames(simtrait) <- names(rTraitDisc(temp_tree, model = "ER", k = 2, rate = er.ml$opt$q12))
        }
        
      }
    
    #Fit models to simulated data
    er.sim <- vector(mode = "list", length = 100)
    ard.sim <- vector(mode = "list", length = 100)
    er_white.sim <- vector(mode = "list", length = 100)
    
    for(sims in 1:100){
      if(length(unique(simtrait[,sims])) == 1){next} #Skipping simulated data where only 1 level is represented
      er.sim[[sims]] <- fitDiscrete(tree, dat=simtrait[, sims, drop = F], model = "ER", ncores = 5)
      er_white.sim[[sims]] <- fitDiscrete(tree, dat=simtrait[, sims, drop = F], model = "ER", transform = "white", ncores = 5)
      ard.sim[[sims]] <- fitDiscrete(tree, dat=simtrait[, sims, drop = F], model = "ARD", ncores = 5)
    }
    
    #Initializing intermediate variables
    sim_param <- data.frame(Param = double())
    sim_max <- c()
    
    #Extracting fits and comparing to original model
    for(m in 1:100){
      if(is_empty(er.sim[[m]]) | is_empty(ard.sim[[m]]) | is_empty(er_white.sim[[m]])){next} #Because we skipped traits where only 1 trait was simulated, to avoide errors with AICc comparisons (does not accept NULL)
      
      sim_res <- data.frame(AICc = c(er.sim[[m]]$opt$aicc, ard.sim[[m]]$opt$aicc, er_white.sim[[m]]$opt$aicc))
      sim_weight <- Weights(sim_res$AICc)
      names(sim_weight) <- c("ER", "ARD", "ER-white")
      sim_max <- append(names(which.max(sim_weight)), sim_max)
      
      #Using K-S test to test goodness-of-fit for estimated parameter and simulated parameters
      if(model_max == "ER"){
        
        sim_param[m, ] <- er.sim[[m]]$opt$q12
        param_conf <- ks.test(er.ml$opt$q12, sim_param$Param)
        
      } else if (model_max == "ARD") {
        
        sim_param[m, ] <- ard.sim[[m]]$opt$q12
        param_conf <- ks.test(ard.ml$opt$q12, sim_param$Param)
        
      }
    }
    
    #What % is the model replicated
    model_conf <- sum(sim_max %in% model_max) / length(sim_max)
    
    #Save everything in output table
    if(model_max == "ER"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- er.ml$opt$q12
        res_df_M[k, 4] <- NA
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
      } else if (i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- er.ml$opt$q12
        res_df_D[k, 4] <- NA
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- param_conf$p.value
        res_df_D[k, 7] <- length(sortedData$phy$tip.label)}
    } else if(model_max == "ARD"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- ard.ml$opt$q12
        res_df_M[k, 4] <- ard.ml$opt$q21
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- param_conf$p.value
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
      } else if(i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- ard.ml$opt$q12
        res_df_D[k, 4] <- ard.ml$opt$q21
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- param_conf$p.value
        res_df_D[k, 7] <- length(sortedData$phy$tip.label)}
    } else if(model_max == "ER-white"){ #No simulations done for white noise models
      if(i == 1){
        res_df_M[k, 1] <- colnames(temp)
        res_df_M[k, 2] <- model_max
        res_df_M[k, 3] <- NA
        res_df_M[k, 4] <- NA
        res_df_M[k, 5] <- model_conf
        res_df_M[k, 6] <- NA
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)
      } else if(i == 2){
        res_df_D[k, 1] <- colnames(temp)
        res_df_D[k, 2] <- model_max
        res_df_D[k, 3] <- NA
        res_df_D[k, 4] <- NA
        res_df_D[k, 5] <- model_conf
        res_df_D[k, 6] <- NA
        res_df_M[k, 7] <- length(sortedData$phy$tip.label)}
    }
    
    print(paste0("Done for: ", paste(k, temp_cat, sep = "-")))
      
  }
res_df_M <- res_df_M[!(is.na(res_df_M$OG)), ]
res_df_D <- res_df_D[!(is.na(res_df_D$OG)), ]

write.table(res_df_M, "/home/s9/tsai/work/Ferenc/Models/binary/binary_fits_pers_out.tsv", quote = F, col.names = T, row.names = F)
write.table(res_df_D, "/home/s9/tsai/work/Ferenc/Models/binary/binary_fits_degr_out.tsv", quote = F, col.names = T, row.names = F)

}


save.image("/home/s9/tsai/work/Ferenc/Models/binary/model_select_MK.RData")

