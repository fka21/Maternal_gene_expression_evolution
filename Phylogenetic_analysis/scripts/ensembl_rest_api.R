##### LIBRARIES #####
library(httr)
library(jsonlite)
library(xml2)
library(Biostrings)
library(taxize)
library(phytools)
library(tidyverse)

#Initializing paths for queries
server <- "http://rest.ensembl.org"
ext1 <- "/homology/id/"
ext2 <- c("?compara=metazoa;", "?compara=vertebrates;") #Surveying both databases
ext3 <- "type=orthologues;sequence=cdna;content-type=text/fasta;aligned=0"

ext4 <- "/sequence/id/"
ext5 <- "?type=cds;multiple_sequences=T"

#Change these for queries, generally major lineage representative homologs were chosen
query <- c("FBgn0025936", "WBGene00006868", "CapteG198909" , "LOC763632", "LOC579430",
           "LOC115928778", "LOC105440970", "LOC105439753", "LOC115928770", "LOC763766",
           "Aqu2.1.16689", "Aqu2.1.37706", "Aqu2.1.37708", "Aqu2.1.03862", "Aqu2.1.05648",
          "Aqu2.1.32612", "Ocbimv22021576m.g")

names <- c()

#Looping through and exporting
for(i in 1:length(query)){
  #Initialize empty objects
  string_lst <- list()
  string_names <- c()
  
  #Extract query gene
  temp_id <- query[i]
    
    #Query only vertebrates or metazoa databases as there is some redundacy with IDs with other databases
    for(k in 1:length(ext2)){
      ext <- paste0(ext1, temp_id, ext2[k], ext3)
      r <- content(GET(paste(server, ext, sep = ""), content_type("application/json")))
      
      #If found in database then skip searching the other one as the other one might not have it
      if(!(purrr::is_empty(r$data))){break} 
      }
    
  #Skip if nothing found
  if(purrr::is_empty(r$data[[1]]$homologies)){next}
  
  #Extracting species name
    species_name <- r$data[[1]]$homologies[[1]]$source$species
    id_name <-  r$data[[1]]$homologies[[1]]$source$id
    filename <- paste0(paste(id_name, species_name, sep = "_"), ".fa")

 
  #Getting CDS
    
    for(j in 1:length(r$data[[1]]$homologies)){
      #Retrieving information of each homolog gene (ID + species of origin)
      temp_id <- r$data[[1]]$homologies[[j]]$target$id
      temp_species <- r$data[[1]]$homologies[[j]]$target$species
      
      #Using that information query ensembl rest api sequences for CDS
      ext_seq <- paste0(ext4, temp_id, ext5, ";", "species=", temp_species)
      t <- content(GET(paste0(server, ext_seq), content_type("application/json")))
      
      iso <- c()
      #Extract longest isoform
      for(k in 1:length(t)){
        iso <- c(iso, width(t[[k]]$seq))
      }
      
      iso_longest <- which.max(iso)
      temp_seq <- DNAString(t[[iso_longest]]$seq)
   
      #Skip entries with undetermined characters
      if(alphabetFrequency(temp_seq, baseOnly=TRUE, as.prob=TRUE)[5] > 0){
        next
      }
      
      #Some names are retrieved not following binomial nomenclature, keeping first 2 names only in these cases
      if(length(str_split(temp_species, pattern = "_")[1]) != 2){
        temp_species <- str_extract_all(temp_species, "^[a-z]+_[a-z]+")
      } 
      
      temp_name <- paste(t[[iso_longest]]$id, str_replace_all(temp_species, "_", "-"), sep = "|")
      
      #Building string object
      string_lst <- append(string_lst, list(temp_seq))
      string_names <- c(string_names, temp_name)
    }
    
    LargeDNAStringSet <- DNAStringSet(string_lst) 
    names(LargeDNAStringSet) <- string_names
    
    writeXStringSet(LargeDNAStringSet, paste0("/Users/ferenckagan/Documents/Bioinformatic_analysis/Gene_trees/EphrinR/", filename))
    
}  

