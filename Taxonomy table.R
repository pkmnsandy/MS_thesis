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

metadata <- psmelt(ps1)

ASV_metadata <- metadata %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1")) %>%
  group_by(OTU, soil, compartment, genotype, nitro) %>%
  summarize(Abundance = sum(Abundance))
  
write.csv(ASV_metadata,"ASV table.csv", row.names = FALSE)


ASV_taxonomy <- metadata %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1")) %>%
  group_by(OTU, Class, Genus) %>%
  summarize(Abundance = sum(Abundance))

write.csv(ASV_metadata,"ASV taxonomy.csv", row.names = FALSE)

summarized_ASV_taxonomy <- ASV_taxonomy %>%
  group_by(Genus) %>%
  summarise_all(funs(paste(na.omit(.), collapse = ", "))) %>%
  select(-c(Class, Abundance))
 
write.csv(summarized_ASV_taxonomy ,"ASV summarized taxonomy.csv", row.names = FALSE)

summarized_ASV_tally <- ASV_taxonomy %>%
  group_by(Genus) %>%
  count(Genus)

summarized_Genus_tally <- ASV_taxonomy %>%
  group_by(Genus) %>%
  summarize(Abundance = sum(Abundance))


write.csv(summarized_ASV_taxonomy ,"ASV tally.csv", row.names = FALSE)

#############################
#############################
#############################
#############################
#all samples
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

inner_join(ps1.summary.all, taxon_pool.all, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other", Genus)) %>%
  group_by(sample.id, soil, compartment, Genus) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Methylocystis", "Methylocystaceae_11","Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             
                             "Methylococcus", "Methylobacter", "Methylococcaceae", "Methylosarcina",
                             "Methylosarcina,_uncultured", "Methylosoma,_uncultured", "RPC_2a", "RPC_2d",
                             "Unclassified_Methylosarcina,_uncultured",
                             "Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = Genus)) +
  geom_col(width = 1.5) +
  facet_grid(~ soil + compartment, space="free_x", scales="free_x") +
  #facet_nested(.~ compartment + nitro, nest_line = element_line(linetype = 2)) +
  scale_fill_manual(values=c("#a48d5b", "#a5a06a", "#a5b27e", "#a5c495", "#a5d6af", #alphaproteobacteria

                             "#5b72a4", "#6c78ab", "#7d7eb1", "#8d84b7", "#9c8abc", "#ab90c0", 
                             "#ba96c5", "#c89dc8", "#d6a4cc", #gammaproteobacteria
                             "#5d5c5d"
  )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "compositions.all.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 1000, unit = "px",
       dpi = 300)


##ASV-level
ps1.melt.all <- psmelt(ps1)

ps1.all.rel_abund <- ps1.melt.all %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

write.csv(ps1.summary.ita,"genus.all.csv", row.names = FALSE)

ps1.summary.all <- ps1.all.rel_abund %>%
  group_by(sample.id, soil, compartment, OTU) %>%
  summarize(rel_abundance = 100*sum(rel_abundance))

taxon_pool.all <- ps1.summary.all %>%
  group_by(OTU) %>%
  summarize(pool = max(rel_abundance) < 10, .groups="drop")

inner_join(ps1.summary.all, taxon_pool.all, by = "OTU") %>%
  mutate(OTU = if_else(pool, "Other", OTU)) %>%
  group_by(sample.id, soil, compartment, OTU) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(OTU, "ASV1", "ASV2459", "ASV2476", "ASV2484", "ASV2772",
                             "ASV618", "ASV138", "ASV515",
                             
                             "ASV2113", "ASV2320", "ASV1642",
                             
                             "Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = OTU)) +
  geom_col(width = 1.5) +
  facet_grid(~ soil + compartment, space="free_x", scales="free_x") +
  #facet_nested(.~ compartment + nitro, nest_line = element_line(linetype = 2)) +
  scale_fill_manual(values=c("#a48d5b", "#a59863", "#a5a26d", "#a5ad78",
                             "#a5b784", "#a5c291", "#a5cca0", "#a5d6af", #alphaproteobacteria
                             
                             "#5b72a4", "#9c8abc", "#d6a4cc", #gammaproteobacteria
  
                             "#5d5c5d"
  )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "compositions.all.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 1000, unit = "px",
       dpi = 300)



#Italian paddy soil subset
ps1.ita <- subset_samples(ps1, soil == "Italian paddy soil")

ps1.genus.ita <- aggregate_taxa(ps1.ita, "Genus")

ps1.melt.ita <- psmelt(ps1.genus.ita)

ps1.ita.rel_abund <- ps1.melt.ita %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

write.csv(ps1.summary.ita,"genus.ita.genotype.csv", row.names = FALSE)

#nitrogen
ps1.summary.ita.nitro <- ps1.ita.rel_abund %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = 100*sum(rel_abundance)) %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1"))

taxon_pool.nitro <- ps1.summary.ita.nitro %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 3, .groups="drop")

inner_join(ps1.summary.ita.nitro, taxon_pool.nitro, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other", Genus)) %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Methylocystis", "Methylocystaceae_11", "Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             
                             "Methylococcus", "Methylococcaceae", "Methylobacter", "Methylosarcina,_uncultured",
                             "Methylosoma,_uncultured", "RPC_2a", "RPC_2d", "Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = Genus)) +
  geom_col(width = 1.1) +
  facet_grid(~ compartment + nitro, space="free_x", scales="free_x") +
  #facet_nested(.~ compartment + nitro, nest_line = element_line(linetype = 2)) +
  scale_fill_manual(values=c("#a48d5b", "#a5a06a", "#a5b27e", "#a5c495", "#a5d6af", #alphaproteobacteria
                             "#5b72a4", "#727aad", "#8782b5", "#9c8abc", "#b092c2", "#c49bc7", "#d6a4cc", #gammaproteobacteria
                             "#5d5c5d" #others
                             )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "bar.italian.nitrogen.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 900, unit = "px",
       dpi = 300)

#gentoype
ps1.summary.ita.genotype <- ps1.ita.rel_abund %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = 100*sum(rel_abundance)) %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1"))

taxon_pool.genotype <- ps1.summary.ita.genotype %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 3, .groups="drop")

inner_join(ps1.summary.ita.genotype, taxon_pool.nitro, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other", Genus)) %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Methylocystis", "Methylocystaceae_11", "Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             
                             "Methylococcus", "Methylococcaceae", "Methylobacter", "Methylosarcina,_uncultured",
                             "Methylosoma,_uncultured", "RPC_2a", "RPC_2d", "Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = Genus)) +
  geom_col(width = 1.1) +
  facet_grid(~ compartment + genotype, space="free_x", scales="free_x") +
  scale_fill_manual(values=c("#a48d5b", "#a5a06a", "#a5b27e", "#a5c495", "#a5d6af", #alphaproteobacteria
                             "#5b72a4", "#727aad", "#8782b5", "#9c8abc", "#b092c2", "#c49bc7", "#d6a4cc", #gammaproteobacteria
                             "#5d5c5d" #others
  )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "bar.italian.genotype.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 900, unit = "px",
       dpi = 300)

taxon.abund.genotype <- ps1.summary.ita.nitro %>%
  group_by(compartment, genotype, Genus) %>%
  summarize(meanR = mean(rel_abundance))

write.csv(taxon.abund.genotype,"genus.ita.genotype.relabund.csv", row.names = FALSE)


taxon.abund.nitrogen <- ps1.summary.ita.nitro %>%
  group_by(compartment, nitro, Genus) %>%
  summarize(meanR = mean(rel_abundance))

write.xlsx(taxon.abund.genotype, file = "taxon_relative_abundance.xlsx",
           sheetName = "ITA_Nitrogen", append = FALSE)


#Philippine paddy soil subset
ps1.ph <- subset_samples(ps1, soil == "Philippine paddy soil")

ps1.genus.ph <- aggregate_taxa(ps1.ph, "Genus")

ps1.melt.ph <- psmelt(ps1.genus.ph)

ps1.ph.rel_abund <- ps1.melt.ph %>%
  group_by(sample.id) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)) %>%
  ungroup()

write.csv(ps1.summary.ph,"genus.ph.genotype.csv", row.names = FALSE)

#nitrogen
ps1.summary.ph.nitro <- ps1.ph.rel_abund %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = 100*sum(rel_abundance)) %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1"))

taxon_pool.nitro <- ps1.summary.ph.nitro %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 3, .groups="drop")

inner_join(ps1.summary.ph.nitro, taxon_pool.nitro, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other", Genus)) %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Methylocystis", "Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             
                             "Methylococcus", "Methylococcaceae", "Methylobacter", "Methylosarcina",
                             "Methylosarcina,_uncultured", "Methylosoma,_uncultured",
                             "RPC_2a", "RPC_2d", "Unclassified_Methylosarcina,_uncultured","Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = Genus)) +
  geom_col(width = 1.1) +
  facet_grid(~ compartment + nitro, space="free_x", scales="free_x") +
  #facet_nested(.~ compartment + nitro, nest_line = element_line(linetype = 2)) +
  scale_fill_manual(values=c("#a48d5b", "#a5a670", "#a5be8d", "#a5d6af", #alphaproteobacteria
                             "#5b72a4", "#6c78ab", "#7d7eb1", "#8d84b7", "#9c8abc",
                             "#ab90c0", "#ba96c5", "#c89dc8", "#d6a4cc", #gammaproteobacteria
                             "#5d5c5d" #others
  )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "bar.ph.nitrogen.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 900, unit = "px",
       dpi = 300)

#gentoype
ps1.summary.ph.genotype <- ps1.ph.rel_abund %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = 100*sum(rel_abundance)) %>%
  mutate(Genus = str_replace(Genus, "Methylococaceae(.*)", "Methylococcaceae\\1"))

taxon_pool.genotype <- ps1.summary.ph.genotype %>%
  group_by(Genus) %>%
  summarize(pool = max(rel_abundance) < 3, .groups="drop")

inner_join(ps1.summary.ph.genotype, taxon_pool.nitro, by = "Genus") %>%
  mutate(Genus = if_else(pool, "Other", Genus)) %>%
  group_by(sample.id, compartment, genotype, nitro, Genus) %>%
  summarize(rel_abundance = sum(rel_abundance)) %>%
  mutate(Genus = fct_relevel(Genus, "Methylocystis", "Methylocystis,_uncultured",
                             "Unclassified_Methylocystaceae", "Unclassified_Methylocystis_echinoides",
                             
                             "Methylococcus", "Methylococcaceae", "Methylobacter", "Methylosarcina",
                             "Methylosarcina,_uncultured", "Methylosoma,_uncultured",
                             "RPC_2a", "RPC_2d", "Unclassified_Methylosarcina,_uncultured","Other")) %>%
  ggplot(aes(x = sample.id, y= rel_abundance, fill = Genus)) +
  geom_col(width = 1.1) +
  facet_grid(~ compartment + genotype, space="free_x", scales="free_x") +
  scale_fill_manual(values=c("#a48d5b", "#a5a670", "#a5be8d", "#a5d6af", #alphaproteobacteria
                             "#5b72a4", "#6c78ab", "#7d7eb1", "#8d84b7", "#9c8abc",
                             "#ab90c0", "#ba96c5", "#c89dc8", "#d6a4cc", #gammaproteobacteria
                             "#5d5c5d" #other
  )) + 
  labs(x = NULL,
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = ggh4x::guide_axis_nested(delim = ".", extend = -1)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "right",
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
        strip.text.x = element_text(color = "black", size = 12),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))

ggsave(plot = last_plot(),
       filename = "bar.ph.genotype.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 3000, height = 900, unit = "px",
       dpi = 300)

taxon.abund.genotype <- ps1.summary.ph.nitro %>%
  group_by(compartment, genotype, Genus) %>%
  summarize(meanR = mean(rel_abundance))

write.csv(taxon.abund.genotype,"genus.ph.genotype.relabund.csv", row.names = FALSE)


taxon.abund.nitrogen <- ps1.summary.ph.nitro %>%
  group_by(compartment, nitro, Genus) %>%
  summarize(meanR = mean(rel_abundance))

write.xlsx(taxon.abund.genotype, file = "taxon_relative_abundance.xlsx",
           sheetName = "ITA_Nitrogen", append = FALSE)

#phylogenetic tree
amplicon_seqs_genus <- tax_glom(ps1, taxrank="Species")

ggtree(amplicon_seqs_genus) +
  geom_tippoint(aes(color=Class), size=1.5) +
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3)
