#Packages needed
library(microViz)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(microbiome)
library(mia)
library(miaViz)
library(scater)
library(ggtree)
library(metacoder)
library(tidyverse)
library(ggtext)
library(forcats)
library(ggh4x)
library(xlsx)
library(ggvenn)

options(scipen=999) 

#create a phyloseq object
#taxonomy file
biom_data <- import_biom(BIOMfilename = "input/taxa_table.biom", treefilename = "input/tree.nwk")

#mapping file
mapping_file <- import_qiime_sample_data(mapfilename = "input/metadata.txt")

# Merge the OTU and mapping data into a phyloseq object
phylo <- merge_phyloseq(biom_data, mapping_file)

#Add names to biom table and check phyloseq objects
colnames(tax_table(phylo))= c("Kingdom","Phylum","Class","Order","Family", "Cluster", "Genus", "Species")
rank_names(phylo)

#cleaning the taxonomy table
tax_table(phylo)[, colnames(tax_table(phylo))]<- gsub(tax_table(phylo)
                                                      [, colnames(tax_table(phylo))], pattern = ".*__", replacement = "")

ps1 <- phylo

#rename ASVs
taxa_names(ps1) <- paste("ASV", 1:ntaxa(ps1), sep="")
taxa_names(ps1)


ps2 <- data.frame(tax_table(ps1))
for (i in 1:8){ ps2[,i] <- as.character(ps2[,i])}
ps2[is.na(ps2)] <- ""

for (i in 1:nrow(ps2)){
  if (ps2[i,2] == ""){
    kingdom <- paste("Unlcassified_", ps2[i,1], sep = "")
    ps2[i, 2:8] <- kingdom
  } else if (ps2[i,3] == ""){
    phylum <- paste("Unclassified_", ps2[i,2], sep = "")
    ps2[i, 3:8] <- phylum
  } else if (ps2[i,4] == ""){
    class <- paste("Unclassified_", ps2[i,3], sep = "")
    ps2[i, 4:8] <- class
  } else if (ps2[i,5] == ""){
    order <- paste("Unclassified_", ps2[i,4], sep = "")
    ps2[i, 5:8] <- order
  } else if (ps2[i,6] == ""){
    family <- paste("Unclassified_", ps2[i,5], sep = "")
    ps2[i, 6:8] <- family
  } else if (ps2[i,7] == ""){
    cluster <- paste("Unclassified_", ps2[i,6], sep = "")
    ps2[i, 7:8] <- cluster
  } else if (ps2[i,8] == ""){
    ps2$Species[i] <- paste("Genus",ps2$Genus[i], sep = "_")
  }
}

tax_table(ps1) <- as.matrix(ps2)

#Examining each ASV with metadata
ntaxa(ps1) #2865
(asv_tab <- data.frame(otu_table(ps1)))

#taxonomic barcharts
#Italian vs Philippine paddy soil
#genus-level
ps1.genus.all <- aggregate_taxa(ps1, "Genus")

ps1.melt.all <- psmelt(ps1.genus.all)

ps1.all.rel_abund <- ps1.melt.all %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

write.csv(ps1.summary.ita,"genus.all.csv", row.names = FALSE)

ps1.summary.all <- ps1.all.rel_abund %>%
  group_by(sample.id, soil, compartment, Genus) %>%
  summarize(rel_abundance = 100*sum(rel_abundance)) %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1"))

taxon_pool.all <- ps1.summary.all %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 3, .groups="drop")

paddy <- inner_join(ps1.summary.all, taxon_pool.all, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(soil, compartment, Genus) %>%
  summarize(rabund = mean(rel_abundance), SD = sd(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus,
                             "Methylococcus", "Methylobacter", "Methylococcaceae", "Methylosarcina",
                             "Methylosarcina,_uncultured", "Methylosoma,_uncultured", "RPC_2a", "RPC_2d",
                             "Unclassified_Methylosarcina,_uncultured",
                             "Methylocystis", "Methylocystaceae_11","Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             "Other"))

paddy %>%
  ggplot(aes(x = Genus, y = rabund, fill = soil, linetype = compartment,
             group = interaction(compartment,soil))) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = rabund-SD, ymax = rabund+SD), linetype = "solid",
                width = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#E3AC37", "#376EAA")) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  ylab("Relative abundance (%)") +
  expand_limits(y = 80) +
  theme_bw() +
  theme(#legend.text = element_text(face = "italic"),
        legend.position = c(0.3, 0.7),
        #strip.placement = "outside",
        #strip.background = element_rect(fill = NA, colour = "white"),
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        #strip.text.x = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 7, angle = 345, hjust = 0.05,
                                   vjust = 0.5, color = "black", face = "bold.italic"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.x = element_blank(),
        strip.background = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.general.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 3000, height = 1500, unit = "px",
       dpi = 300)

#Type I vs type II
#all samples
#genus-level
ps1.class.all <- aggregate_taxa(ps1, "Class")

ps1.melt.all.class <- psmelt(ps1.class.all)

ps1.all.rel_abund.class <- ps1.melt.all.class %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.summary.all.class <- ps1.all.rel_abund.class %>%
  group_by(sample.id, soil, compartment, Class) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

taxon_pool.all.class <- ps1.summary.all.class %>%
  group_by(Class) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.summary.all.class, taxon_pool.all.class, by = "Class") %>%
  mutate(Class = if_else(pool, "Other (<1%)", Class)) %>%
  group_by(soil, Class) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = soil, y= rel_abundance,
             fill = factor(Class, level = c("Other (<1%)",
                                            "Alphaproteobacteria",
                                            "Gammaproteobacteria")))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", "#f67738", "#00898a")) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.spacing = unit(0,"cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(color="black", size=14),
    axis.text.y = element_text(color="black", size = 10),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.general.class.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1800, height = 500, unit = "px",
       dpi = 300)


################################################################################
################################################################################
#Italian paddy soil
#Rhizosphere samples

level_order <- c("IR64", "Kasalath", "Nipponbare", "Rufi", "Bulk soil")

ps1.ita.rz <- subset_samples(ps1,
                             soil == "Italian paddy soil" & compartment %in% c("Rhizosphere"))

ps1.ita.rz.genus <- aggregate_taxa(ps1.ita.rz, "Genus")

ps1.ita.rz.genus.melt <- psmelt(ps1.ita.rz.genus)

ps1.ita.rz.genus.melt.relabund <- ps1.ita.rz.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.ita.rz.genus.melt.relabund.summary <- ps1.ita.rz.genus.melt.relabund %>%
  group_by(sample.id, genotype, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.ita.rz.genus.melt.relabund.summary.pooled <- ps1.ita.rz.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.ita.rz.genus.melt.relabund.summary, ps1.ita.rz.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(genotype, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(genotype, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               #"#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               "#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               #"#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               "#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               #"#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               "#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               #"#00998f", #"Methylomicrobium"
                               "#00898a", #"Methylococcus"
                               #"#007882", #"Methylococcacea_24a"
                               "#176877"#, #"Methylococcaceae"
                               #"#245769", #"Methylocaldum"
                               #"#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.ita.rz.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)

#Root samples
ps1.ita.rt <- subset_samples(ps1,
                             soil == "Italian paddy soil" & compartment %in% c("Root"))

ps1.ita.rt.genus <- aggregate_taxa(ps1.ita.rt, "Genus")

ps1.ita.rt.genus.melt <- psmelt(ps1.ita.rt.genus)

ps1.ita.rt.genus.melt.relabund <- ps1.ita.rt.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.ita.rt.genus.melt.relabund.summary <- ps1.ita.rt.genus.melt.relabund %>%
  group_by(sample.id, genotype, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.ita.rt.genus.melt.relabund.summary.pooled <- ps1.ita.rt.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.ita.rt.genus.melt.relabund.summary, ps1.ita.rt.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(genotype, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(genotype, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               "#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               "#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               #"#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               "#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               #"#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               #"#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               #"#00998f", #"Methylomicrobium"
                               "#00898a", #"Methylococcus"
                               "#007882", #"Methylococcacea_24a"
                               "#176877", #"Methylococcaceae"
                               #"#245769", #"Methylocaldum"
                               "#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.ita.rt.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)

#Bulk soil samples
ps1.ita.blk <- subset_samples(ps1,
                              soil == "Italian paddy soil" & compartment %in% c("Bulk soil"))

ps1.ita.blk.genus <- aggregate_taxa(ps1.ita.blk, "Genus")

ps1.ita.blk.genus.melt <- psmelt(ps1.ita.blk.genus)

ps1.ita.blk.genus.melt.relabund <- ps1.ita.blk.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.ita.blk.genus.melt.relabund.summary <- ps1.ita.blk.genus.melt.relabund %>%
  group_by(sample.id, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.ita.blk.genus.melt.relabund.summary.pooled <- ps1.ita.blk.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.ita.blk.genus.melt.relabund.summary, ps1.ita.blk.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(compartment, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(compartment, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               #"#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               "#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               #"#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               #"#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               #"#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               "#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               "#00998f", #"Methylomicrobium"
                               "#00898a"#, #"Methylococcus"
                               #"#007882", #"Methylococcacea_24a"
                               #"#176877", #"Methylococcaceae"
                               #"#245769", #"Methylocaldum"
                               #"#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.ita.blk.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 600, unit = "px",
       dpi = 300)

################################################################################
################################################################################
#Philippine paddy soil
#rhizosphere
level_order <- c("IR64", "Kasalath", "Nipponbare", "Rufi", "Bulk soil")

ps1.rz.php <- subset_samples(ps1,
                          soil == "Philippine paddy soil" & compartment %in% c("Rhizosphere"))

ps1.rz.php.genus <- aggregate_taxa(ps1.rz.php, "Genus")

ps1.rz.php.genus.melt <- psmelt(ps1.rz.php.genus)

ps1.rz.php.genus.melt.relabund <- ps1.rz.php.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.rz.php.genus.melt.relabund.summary <- ps1.rz.php.genus.melt.relabund %>%
  group_by(sample.id, genotype, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.rz.php.genus.melt.relabund.summary.pooled <- ps1.rz.php.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.rz.php.genus.melt.relabund.summary, ps1.rz.php.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(genotype, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(genotype, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               #"#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               "#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               #"#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               #"#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               #"#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               "#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               "#00998f", #"Methylomicrobium"
                               "#00898a", #"Methylococcus"
                               #"#007882", #"Methylococcacea_24a"
                               "#176877"#, #"Methylococcaceae"
                               #"#245769", #"Methylocaldum"
                               #"#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.php.rz.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)

#Root samples
ps1.rt.php <- subset_samples(ps1,
                             soil == "Philippine paddy soil" & compartment %in% c("Root"))

ps1.rt.php.genus <- aggregate_taxa(ps1.rt.php, "Genus")

ps1.rt.php.genus.melt <- psmelt(ps1.rt.php.genus)

ps1.rt.php.genus.melt.relabund <- ps1.rt.php.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.rt.php.genus.melt.relabund.summary <- ps1.rt.php.genus.melt.relabund %>%
  group_by(sample.id, genotype, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.rt.php.genus.melt.relabund.summary.pooled <- ps1.rt.php.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.rt.php.genus.melt.relabund.summary, ps1.rt.php.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(genotype, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(genotype, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               #"#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               "#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               "#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               #"#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               #"#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               "#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               "#00998f", #"Methylomicrobium"
                               "#00898a", #"Methylococcus"
                               "#007882", #"Methylococcacea_24a"
                               "#176877", #"Methylococcaceae"
                               "#245769", #"Methylocaldum"
                               "#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.php.rt.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)

#Bulk soil samples
ps1.php.blk <- subset_samples(ps1,
                              soil == "Philippine paddy soil" & compartment %in% c("Bulk soil"))

ps1.php.blk.genus <- aggregate_taxa(ps1.php.blk, "Genus")

ps1.php.blk.genus.melt <- psmelt(ps1.php.blk.genus)

ps1.php.blk.genus.melt.relabund <- ps1.php.blk.genus.melt %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

ps1.php.blk.genus.melt.relabund.summary <- ps1.php.blk.genus.melt.relabund %>%
  group_by(sample.id, compartment, Genus) %>%
  summarize(rel_abundance = 100*mean(rel_abundance))

ps1.php.blk.genus.melt.relabund.summary.pooled <- ps1.php.blk.genus.melt.relabund.summary %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 1, .groups="drop")

inner_join(ps1.php.blk.genus.melt.relabund.summary, ps1.php.blk.genus.melt.relabund.summary.pooled, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other (<1%)", Genus)) %>%
  group_by(compartment, Genus) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ggplot(aes(x = factor(compartment, level = level_order), y= rel_abundance,
             fill = factor(Genus, level = c("Other (<1%)",
                                            #type II
                                            "Unclassified_Methylocystis,_uncultured",
                                            "Unclassified_Methylocystis_echinoides",
                                            "Unclassified_Methylocystaceae",
                                            "pmoA2_like_4",
                                            "Methylocystis,_uncultured",
                                            "Methylocystis",
                                            "Methylocystaceae_11",
                                            #type I
                                            "Unclassified_Methylosarcina,_uncultured",
                                            "RPC1_3_like_10,_RPC1",
                                            "RPC_2d",
                                            "RPC_2a",
                                            "Methylosoma,_uncultured",
                                            "Methylosarcina,_uncultured",
                                            "Methylosarcina",
                                            "Methylomicrobium",
                                            "Methylococcus",
                                            "Methylococcaceae_24a",
                                            "Methylococcaceae",
                                            "Methylocaldum",
                                            "Methylobacter"
             )))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(name = "Mean relative abundance (%)", 
                     labels = scales::label_percent()) +
  scale_fill_manual(values = c("black", #"Other (<1%)"
                               #type II
                               "#ff8c00", #"Unclassified_Methylocystis,_uncultured"
                               #"#fc8126", #"Unclassified_Methylocystis_echinoides"
                               "#f67738", #"Unclassified_Methylocystaceae"
                               #"#ee6f46", #"pmoA2_like_4"
                               "#e46851", #"Methylocystis,_uncultured"
                               "#d8635a", #"Methylocystis"
                               #"#ca6061", #"Methylocystaceae_11"
                               #type I
                               "#fafa6e", #"Unclassified_Methylosarcina,_uncultured"
                               "#d1f072", #"RPC1_3_like_10,_RPC1"
                               "#aae479", #"RPC_2d"
                               "#86d780", #"RPC_2a"
                               #"#64c987", #"Methylosoma,_uncultured"
                               "#44b98d", #"Methylosarcina,_uncultured"
                               "#23aa8f", #"Methylosarcina"
                               "#00998f", #"Methylomicrobium"
                               "#00898a", #"Methylococcus"
                               #"#007882", #"Methylococcacea_24a"
                               #"#176877", #"Methylococcaceae"
                               #"#245769", #"Methylocaldum"
                               "#2a4858" #"Methylobacter"
  )) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=14),
        axis.text.y = element_text(color="black", size = 10),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.title.y = element_blank())

ggsave(plot = last_plot(),
       filename = "compositions.php.blk.genus.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 2000, height = 600, unit = "px",
       dpi = 300)
