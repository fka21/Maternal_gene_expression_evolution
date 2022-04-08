##### LIBRARIES #####
library(phangorn)
library(phytools)
library(ape)
library("rjson")
library(tidyverse)
library(ggmsa)
library(ggtree)
library(viridis)
library(taxize)

##### HYPHY ANALYSIS #####
#Import
test <- fromJSON(file = "~/Documents/Gene_expr_evol/Gene_trees/Nanos_FEL.json")
as.data.frame(test$MLE$content$`0`)

#Initialize variables
estims <- data.frame(alpha = double(),
                     beta = double(),
                     alpha_beta = double(),
                     LRT = double(),
                     p_value = double(),
                     Total_branch_length = double())

#Extract results from JSON results file and assign to initialized variables
for(i in 1:length(test$MLE$content$`0`)){
  temp <- unlist(test$MLE$content$`0`)
  
  estims[i, 1] <- temp[seq(1, length(temp), 6)][i]
  estims[i, 2] <- temp[seq(2, length(temp), 6)][i]
  estims[i, 3] <- temp[seq(3, length(temp), 6)][i]
  estims[i, 4] <- temp[seq(4, length(temp), 6)][i]
  estims[i, 5] <- temp[seq(5, length(temp), 6)][i]
  estims[i, 6] <- temp[seq(6, length(temp), 6)][i]
  
}

#Formating estimate tables
estims <- tibble(estims) %>% mutate(site = 1:length(estims$alpha)) %>%
  group_by(site) %>% 
  mutate(category = case_when(alpha > beta ~ "Purifying",
                              beta > alpha ~ "Diversifying"))

#Subset and format for plotting
asd <- estims %>% filter(p_value <= 0.001) %>%
  mutate(asd = case_when(category == "Purifying" ~ alpha,
                         category == "Diversifying" ~ -beta))

#Plotting
svg("~/Desktop/Figures/Nanos_FEL.svg", width = 750, height = 550)
png("~/Desktop/Figures/Nanos_FEL.png",
    width=750, height=550)
ggplot(asd, aes(x=site, y=asd)) +
  geom_segment( aes(x=site, xend=site, y=0, yend=asd), color = "grey", alpha = 0.8) +
  geom_point(size = 2, aes(color = log(Total_branch_length)), alpha = 0.8) +
  theme_light() +
  theme(
    axis.ticks.x = element_blank(),
    text = element_text(size = 15)
  ) +
  labs(color = "log(TBL)") +
  xlab("Site") +
  ylab("") +
  scale_color_viridis(limits = c(0, 5), oob = scales::squish,
                      labels = c(as.character(0:4), ">= 5"), breaks = 0:5) +
  scale_y_continuous(limits = c(-2,10), oob=scales::squish, 
                     breaks = c(-2, 0, 2.5, 5, 7.5, 10), 
                     labels = c("2", "0", "2.5", "5", "7.5", ">=10")) +
  geom_segment(aes(x = 1480, xend = 1560, y = -1, yend = -1), size = 6) +
  annotate("text", label = "zf-nanos", x = 1520, y = -1.5) #These annotation was achieved with hmmscan of the domains

dev.off()

###### GENE TREE RECONCILIATION/PLOT #####
#Read in tree without outgroups as outgroups exert LBA
tree_seq_nwk <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Nanos/Nanos_RMDupl_hmm-filtered_RMDupl.msa.contree")
p <- ggtree(tree_seq_nwk) 

svg("~/Desktop/Figures/Nanos_tree.svg")
p + theme_tree2()
dev.off()

p <- p + geom_tiplab(size = 0.8) +                # labels the tips of all branches with the sample name in the tree file
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = node),
    size = 2,
    color = "darkred",
    hjust = 1,
    vjust = 1)

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
tree_seq_nwk <- root(tree_seq_nwk, outgroup = "Aqu2.1.39257_001_amphimedon-queenslandica", resolve.root = T) #Based on species tree

is.rooted(tree_seq_nwk)
is.binary(species_tree$phylo)
is.rooted(species_tree$phylo)

#Export gene tree and species tree for reconciliations
write.tree(species_tree$phylo, "~/Desktop/sp_tree_notung.nwk")
write.tree(tree_seq_nwk, "~/Desktop/gn_tree_notung.nwk")

#Perform gene tree reconciliation with Notung 2.8
#java -jar Notung/Notung-2.9.1.5.jar -s ~/Desktop/sp_tree_notung.nwk -g ~/Desktop/gn_tree_notung.nwk --treestats --parsable --outputdir ~/Desktop/ --treeoutput newick --speciestag postfix --exact-losses --reconcile

#Import reconciliation table which has been subsetted for summary table
recon <- read_tsv("~/Documents/Gene_expr_evol/Gene_trees/Nanos_recon_summary.tbl")
