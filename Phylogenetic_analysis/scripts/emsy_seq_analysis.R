##### LIBRARIES #####
library(phangorn)
library(phytools)
library(ape)
library("rjson")
library(tidyverse)
library(viridis)
library(treeio)

##### HYPHY ANALYSIS #####
#Import
test <- fromJSON(file = "~/Documents/Gene_expr_evol/Gene_trees/Emsy_FEL.json")
contrast_fel <- fromJSON(file = "~/Documents/Gene_expr_evol/Gene_trees/Emsy.msa.FEL2.json")

#Initialize variables
estims <- data.frame(alpha = double(),
                     beta = double(),
                     alpha_beta = double(),
                     LRT = double(),
                     p_value = double(),
                     Total_branch_length = double(),
                     q_value_contrast = double(),
                     beta_ref = double(),
                     beta_test = double())

#Extract results from JSON results file and assign to initialized variables
for(i in 1:length(test$MLE$content$`0`)){
  temp <- unlist(test$MLE$content$`0`)
  temp_contrast <- unlist(contrast_fel$MLE$content$`0`)
  
  estims[i, 1] <- temp[seq(1, length(temp), 6)][i]
  estims[i, 2] <- temp[seq(2, length(temp), 6)][i]
  estims[i, 3] <- temp[seq(3, length(temp), 6)][i]
  estims[i, 4] <- temp[seq(4, length(temp), 6)][i]
  estims[i, 5] <- temp[seq(5, length(temp), 6)][i]
  estims[i, 6] <- temp[seq(6, length(temp), 6)][i]
  estims[i, 7] <- temp_contrast[seq(7, length(temp_contrast), 9)][i]
  estims[i, 8] <- temp_contrast[seq(2, length(temp_contrast), 9)][i]
  estims[i, 9] <- temp_contrast[seq(3, length(temp_contrast), 9)][i]
  
}

#Formating estimate tables
estims <- tibble(estims) %>% mutate(site = 1:length(estims$alpha)) %>%
  group_by(site) %>% 
  mutate(category = case_when(alpha > beta ~ "Purifying",
                              beta > alpha ~ "Diversifying")) %>%
  mutate(contrast = case_when(q_value_contrast <= 0.05 ~ "Signif",
                              q_value_contrast > 0.05 ~ "Non-signif")) %>%
  mutate(beta_diff = beta_ref - beta_test)

#Subset and format for plotting
asd <- estims %>% filter(p_value <= 0.001) %>%
  mutate(asd = case_when(category == "Purifying" ~ alpha,
                         category == "Diversifying" ~ -beta))

#Plotting
svg("~/Desktop/Figures/Emsy_FEL_zoom.svg",
    width=12, 
    height=9, 
    pointsize=12)
tiff("~/Desktop/Figures/Emsy_FEL.tiff",
    width=750, height=550)
print(ggplot(asd, aes(x=site, y=asd, color = log(Total_branch_length))) +
        geom_segment( aes(x=site, xend=site, y=0, yend=asd), color = "grey", alpha = 0.8) +
        geom_point(size = 3, aes(color = log(Total_branch_length), shape = contrast), alpha = 0.8) +
        theme_light() +
        theme(
          axis.ticks.x = element_blank(),
          text = element_text(size = 15)
        ) +
        labs(color = "log(TBL)") +
        xlab("Site") +
        ylab("") +
        coord_cartesian(ylim = c(-2, 4), xlim = c(590, 710)) +
        scale_color_viridis_b(option = "magma")+
        #scale_color_viridis() +
  geom_segment(aes(x = 600, xend = 700, y = -1, yend = -1), size = 6) +
  annotate("text", label = "ENT domain", x = 650, y = -1.5))  #These annotation was achieved with hmmscan of the domains against the MSA 
dev.off()

#Extracting and aggregating contrast-FEL results for ENT domain sites
ent <- estims[estims$site <= 700 & estims$site >= 600, ]
ent <- ent %>% dplyr::group_by(contrast) %>% dplyr::summarise(prop = n())
#Extracting and aggregating contrast-FEL results for sites outside ENT domain
non_ent <- estims[estims$site > 700 | estims$site < 600, ]
non_ent <- non_ent %>% dplyr::group_by(contrast) %>% dplyr::summarise(prop = n())

#Plotting contrast-FEL results significance as pie charts 
svg("~/Desktop/Figures/non-ent_prop.svg")
print(ggplot(non_ent, aes(x="", y=prop, fill=contrast)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_brewer(palette="Set1"))
dev.off()

svg("~/Desktop/Figures/ent_prop.svg")
print(ggplot(ent, aes(x="", y=prop, fill=contrast)) +
        geom_bar(stat="identity", width=1, color = "white") +
        coord_polar("y", start=0) +
        theme_void() +
        scale_fill_brewer(palette="Set1"))
dev.off()

##### RELAX #####
#Hyphy vision was used to generate plot for these results

##### COPHYLOGENETIC PLOTS #####
#This section was used to compare chromo domain based gene tree to translated CDS based gene tree
tree_1 <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/ent_domains.msa.contree")
tree_1$tip.label <- str_remove_all(tree_1$tip.label, "\\|priapulus-caudatus$")
tree_2 <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Emsy_hmm-filtered_RMDupl.msa.contree")

asd <- cophylo(tree_1, tree_2, rotate = T)
svg("~/Desktop/Figures/emsy_domain_vs_cds.svg")
print("~/Desktop/Figures/emsy_domain_vs_cds.png", height = 850, width = 750)
print(plot(asd, assoc.type="curved",link.lwd=1,link.lty="solid",
           link.col=make.transparent("blue",0.1),fsize=0.8, ftype = "off"))
dev.off()


###### GENE TREE RECONCILIATION/PLOT #####
#LBA trimming method as suggested by Siddall and Whiting (1999) and Pol and Siddall (2001)
#Same steps were used in the other 2 gene trees

#Import gene tree with outgroups
test.tree <- read.tree("/home/ferenkagan/Emsy_with_outgroup_trim.msa.contree")
#Rooting at outgroup
test.tree <- root(test.tree, outgroup = test.tree$tip.label[grepl("outgroup", test.tree$tip.label)], resolve.root = T)

#Which are the longest branch lengths
head(test.tree$edge[order(test.tree$edge.length, decreasing = T),], 6)

#Labeling nodes to identify which nodes to drop 
test.tree$node.label <- paste0("N", 1:Nnode(test.tree))
tips <- adephylo::listTips(test.tree)
tips_to_drop <- names(c(tips$N255, tips$N248, tips$N36))
tips_to_drop <- names(c(tips$N248, tips$N257))
tips <- test.tree$tip.label[!(test.tree$tip.label %in% tips_to_drop)]

#The tips without longest branches were used to run another gene tree estimate
#and inspected if removing these tips repositions long branches from basal positions
write.table(data.frame(X1 = tips), "~/Emsy_IDs_without_LBA.txt", quote = F, col.names = F, row.names = F)

#Loading long branch trimmed gene tree estimates and inspected the general topology of it
tree <- read.tree("~/Emsy_with_outgroup_trim_LBA_2.msa.contree")
ggtree(tree) + geom_tiplab(size = 0.8)

#Import gene tree without outgroups as outgroups exert LBA
tree_seq_nwk <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Emsy_hmm-filtered_RMDupl.msa.contree")

p <- ggtree(tree_seq_nwk) 

svg("~/Desktop/Figures/emsy_tree.svg")
print(p + theme_tree2() + geom_tiplab(size = 1))
dev.off()

p <- p + geom_tiplab(size = 0.8) +                # labels the tips of all branches with the sample name in the tree file
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = label),
    size = 2,
    color = "darkred",
    hjust = 1,
    vjust = 1)

viewClade(p, node = 357)

p1 <- ggtree(tree_seq_nwk) +
  geom_cladelab(node=357, label="Reference", align=F, angle=270, 
                hjust='center', offset.text=0.03, barsize=1.5, fontsize=5) 

##### RECONCILIATION #####
#Generating species tree for species tree from NCBI Taxonomy database
test <- str_remove_all(tree_seq_nwk$tip.label, ".*\\|") %>%
  str_replace_all("-", "_") %>%
  unique() %>%
  classification(db = "ncbi")

species_tree <- class2tree(test, check = TRUE)
species_tree$phylo$tip.label <- str_replace_all(str_to_lower(species_tree$phylo$tip.label), " ", "-")
#Inspect species tree
ggtree(species_tree$phylo, layout = "circular") +
  geom_tiplab(size = 2)

#Used a rooted gene tree for reconciliation
tree_seq_nwk$tip.label <- str_replace_all(tree_seq_nwk$tip.label, "\\|", "_")
tree_seq_nwk <- drop.tip(tree_seq_nwk, tip = "ENSNVIT00000000553_neovison-vison") #Not in species tree, dropping it
tree_seq_nwk <- root(tree_seq_nwk, outgroup = "ML02756a-RA_mnemiopsis-leidyi", resolve.root = T) #Based on species tree

is.rooted(tree_seq_nwk)
is.binary(species_tree$phylo)
is.rooted(species_tree$phylo)

#Export gene tree and species tree for reconciliations
write.tree(species_tree$phylo, "~/Desktop/sp_tree_notung.nwk")
write.tree(tree_seq_nwk, "~/Desktop/gn_tree_notung.nwk")

#Perform gene tree reconciliation with Notung 2.8
#java -jar Notung/Notung-2.9.1.5.jar -s ~/Desktop/sp_tree_notung.nwk -g ~/Desktop/gn_tree_notung.nwk --treestats --parsable --outputdir ~/Desktop/ --treeoutput newick --speciestag postfix --exact-losses --reconcile

#Import reconciliation table which has been subsetted for summary table
recon <- read_tsv("~/Documents/Gene_expr_evol/Gene_trees/Emsy_recon_summary.tbl")
