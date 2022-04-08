##### LIBRARIES #####
library(phangorn)
library(phytools)
library(ape)
library("rjson")
library(tidyverse)
library(viridis)

##### HYPHY ANALYSIS #####
#Import
test <- fromJSON(file = "~/Documents/Gene_expr_evol/Gene_trees/Cbx1_FEL.json")
contrast_fel <- fromJSON(file = "~/Documents/Gene_expr_evol/Gene_trees/Cbx1/cbx1_contrast-fel")

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
svg("~/Desktop/Figures/Cbx1_FEL_zoom.svg")
print(ggplot(asd, aes(x=site, y=asd)) +
          geom_segment( aes(x=site, xend=site, y=0, yend=asd), color = "grey", alpha = 0.8) +
          geom_point(size = 2.5, aes(color = log(Total_branch_length), shape = contrast), alpha = 0.8) +
          theme_light() +
          theme(
            axis.ticks.x = element_blank(),
            text = element_text(size = 15)
          ) +
          labs(color = "log(TBL)") +
          xlab("Site") +
          ylab("") +
          coord_cartesian(ylim = c(-2, 4), xlim = c(680, 770)) +
          scale_color_viridis_b(option = "magma") +
          geom_segment(aes(x = 700, xend = 750, y = -1, yend = -1), size = 6) +
          geom_segment(aes(x = 2820, xend = 2880, y = -1, yend = -1), size = 6) +
          annotate("text", label = "Chromo domain 1", x = 725, y = -1.5) + #These annotation was achieved with hmmscan of the domains
          annotate("text", label = "Chromo domain 2", x = 2820, y = -1.5))  #against the MSA 
dev.off()

svg("~/Desktop/Figures/Cbx1_FEL.svg", width = 750, height = 550)
print(ggplot(asd, aes(x=site, y=asd, color = log(Total_branch_length))) +
        geom_bar(stat = 'identity', width = 0.5, alpha = 0.5) +
        theme_light() +
        theme(
          axis.ticks.x = element_blank(),
          text = element_text(size = 15)
        ) +
        labs(color = "log(TBL)") +
        xlab("Site") +
        ylab("") +
        coord_cartesian(ylim = c(-2, 4)) +
        scale_color_viridis() +
        geom_segment(aes(x = 700, xend = 750, y = -1, yend = -1), size = 6) +
        geom_segment(aes(x = 2820, xend = 2880, y = -1, yend = -1), size = 6) +
        annotate("text", label = "Chromo domain 1", x = 725, y = -1.5) + #These annotation was achieved with hmmscan of the domains
        annotate("text", label = "Chromo domain 2", x = 2820, y = -1.5))  #against the MSA 
dev.off()

#Extracting and aggregating contrast-FEL results for chromo domain sites
cbx <- estims[estims$site <= 750 & estims$site >= 700, ]
cbx <- cbx %>% dplyr::group_by(contrast) %>% dplyr::summarise(prop = n())
#Extracting and aggregating contrast-FEL results for sites outside chromo domain
non_cbx <- estims[estims$site > 750 | estims$site < 700, ]
non_cbx <- non_cbx %>% dplyr::group_by(contrast) %>% dplyr::summarise(prop = n())

#Plotting contrast-FEL results significance as pie charts 
svg("~/Desktop/Figures/non-cbx_prop.svg")
print(ggplot(non_cbx, aes(x="", y=prop, fill=contrast)) +
        geom_bar(stat="identity", width=1, color = "white") +
        coord_polar("y", start=0) +
        theme_void() +
        scale_fill_brewer(palette="Set1"))
dev.off()

svg("~/Desktop/Figures/cbx_prop.svg")
print(ggplot(cbx, aes(x="", y=prop, fill=contrast)) +
        geom_bar(stat="identity", width=1, color = "white") +
        coord_polar("y", start=0) +
        theme_void() +
        scale_fill_brewer(palette="Set1"))
dev.off()


##### COPHYLOGENETIC PLOTS #####
#This section was used to compare chromo domain based gene tree to translated CDS based gene tree
tree_1 <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/chromo_hits.msa.contree")
tree_2 <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Cbx1_domain.msa.contree")

asd <- cophylo(tree_1, tree_2, rotate = T)
svg("~/Desktop/Figures/cbx_domain_vs_cds.svg")
print(plot(asd, assoc.type="curved",link.lwd=1,link.lty="solid",
           link.col=make.transparent("blue",0.1),fsize=0.8, ftype = "off"))
dev.off()

###### GENE TREE RECONCILIATION/PLOT #####
#Reading in consensus tree estimates with outgroups
tree_seq_nwk <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Cbx1_with_outgrou_trim.msa.contree")
#Rooting at outgroups
tree_seq_nwk <- root(tree_seq_nwk, outgroup = tree_seq_nwk$tip.label[grepl("outgroup", tree_seq_nwk$tip.label)], resolve.root = T)

#Inspecting tree
ggtree(tree_seq_nwk) + geom_tiplab2(size = 0.7)

#Read in tree without outgroups as outgroups exert LBA
tree_seq_nwk <- read.tree("~/Documents/Gene_expr_evol/Gene_trees/Cbx1_domain.msa.contree")
svg("~/Desktop/Figures/cbx1_tree.svg")
print(ggtree(tree_seq_nwk) + theme_tree2())
dev.off()

p <- ggtree(tree_seq_nwk) 
p <- p + geom_tiplab(size = 0.7) +                # labels the tips of all branches with the sample name in the tree file
  geom_text2(
    mapping = aes(subset = !isTip,
                  label = node),
    size = 1,
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
tree_seq_nwk <- root(tree_seq_nwk, outgroup = "ML019137a-RA_mnemiopsis-leidyi", resolve.root = T) #Based on species tree

is.rooted(tree_seq_nwk)
is.binary(species_tree$phylo)
is.rooted(species_tree$phylo)

#Export gene tree and species tree for reconciliations
write.tree(species_tree$phylo, "~/Desktop/sp_tree_notung.nwk")
write.tree(tree_seq_nwk, "~/Desktop/gn_tree_notung.nwk")

#Perform gene tree reconciliation with Notung 2.8
#java -jar Notung/Notung-2.9.1.5.jar -s ~/Desktop/sp_tree_notung.nwk -g ~/Desktop/gn_tree_notung.nwk --treestats --parsable --outputdir ~/Desktop/ --treeoutput newick --speciestag postfix --exact-losses --reconcile

#Import reconciliation table which has been subsetted for summary table
stats <- read_tsv("~/Documents/Gene_expr_evol/Gene_trees/Cbx1_recon_summary.tbl")
View(stats)

