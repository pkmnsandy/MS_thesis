library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(ape)
library(picante)
library(car)
library(mia)
library(vegan)
library(factoextra)
library(pairwiseAdonis)
library(tidyverse)
library(metagMisc)
library(microViz)
library(ggplotify)

set.seed(999)

#making a phyloseq object
f <- function(z)sd(z)/sqrt(length(z)) #standard error
#taxonomy file
biom_data <- import_biom(BIOMfilename = "input/taxa_table.biom", treefilename = "input/tree.nwk")

#mapping file
mapping_file <- import_qiime_sample_data(mapfilename = "input/metadata.txt")

#reading tree file
tree_file <- read.tree("input/tree.nwk")

# Merge the OTU and mapping data into a phyloseq object
phylo <- merge_phyloseq(biom_data, mapping_file)

#Add names to biom table and check phyloseq objects
colnames(tax_table(phylo))= c("Kingdom","phylum","Class","Order","Family", "Cluster", "Genus", "Species")
rank_names(phylo)

print(phylo)

#cleaning the taxonomy table
tax_table(phylo)[, colnames(tax_table(phylo))]<- gsub(tax_table(phylo)
                                                      [, colnames(tax_table(phylo))], pattern = ".*__", replacement = "")

datatable(tax_table(phylo))

ps1 <- phylo

#rename ASVs
taxa_names(ps1) <- paste("ASV", 1:ntaxa(ps1), sep="")
datatable(tax_table(ps1))


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

#general examination
compartment %in% c("Rhizosphere", "Bulk soil"))

phylo.general.rclr <- microbiome::transform(phylo, "rclr")

#extract metadata
phylo.general.meta <- meta(phylo)

#Generate distance matrix/calculate robust Aitchison distance
phylo.general.dist <- phyloseq::distance(phylo.general.rclr, method = "euclidean")

#ADONIS test
set.seed(999)
phylo.general.permanova <- vegan::adonis2(phylo.general.dist ~ 
                                            phyloseq::sample_data(phylo.general.rclr)
                                          $soil*phyloseq::sample_data(phylo.general.rclr)$compartment,
                                          permutations = 9999)
phylo.general.permanova


#pairwise adonis
pairwise.adonis(x = phylo.general.dist,
                factors = sample_data(phylo.general.rclr)$genotype,
                p.adjust.m = "bonferroni",
                perm = 9999)

#calculate beta dispersion
#genotype
set.seed(999)
dispr.general.gen <- vegan::betadisper(phylo.general.dist,
                                       phyloseq::sample_data(phylo.general.rclr)$compartment,
                                       type = "centroid")

#permutest
permutest(dispr.general.gen, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.general.gen.dist <- data.frame(dispr.general.gen$distances)

#combine data with metadata
dispr.general.gen.dist <- tibble::rownames_to_column(dispr.general.gen.dist, "sample.id")
dispr.general.gen.plot <- dispr.general.gen.dist %>% left_join(phylo.general.meta)


#Italian soil subset
#rhizosphere
phylo.ita.rz <- subset_samples(ps1, soil == "Italian paddy soil" & 
                              compartment %in% c("Rhizosphere", "Bulk soil"))

phylo.ita.rz.rclr <- microbiome::transform(phylo.ita.rz, "rclr")

#extract metadata
phylo.ita.rz.meta <- meta(phylo.ita.rz)

#Generate distance matrix/calculate robust Aitchison distance
phylo.ita.rz.dist <- phyloseq::distance(phylo.ita.rz.rclr, method = "euclidean")

#ADONIS test
set.seed(999)
phylo.ita.rz.permanova <- vegan::adonis2(phylo.ita.rz.dist ~ 
                                     phyloseq::sample_data(phylo.ita.rz.rclr)
                                     $genotype*phyloseq::sample_data(phylo.ita.rz.rclr)$nitro,
                                     permutations = 9999)
phylo.ita.rz.permanova


#pairwise adonis
pairwise.adonis(x = phylo.ita.rz.dist,
                factors = sample_data(phylo.ita.rz.rclr)$genotype,
                p.adjust.m = "bonferroni",
                perm = 9999)

#calculate beta dispersion
#genotype
set.seed(999)
dispr.ita.rz.gen <- vegan::betadisper(phylo.ita.rz.dist,
                                      phyloseq::sample_data(phylo.ita.rz.rclr)$genotype,
                                      type = "centroid")

#permutest
permutest(dispr.ita.rz.gen, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.ita.rz.gen.dist <- data.frame(dispr.ita.rz.gen$distances)

#combine data with metadata
dispr.ita.rz.gen.dist <- tibble::rownames_to_column(dispr.ita.rz.gen.dist, "sample.id")
dispr.ita.rz.gen.plot <- dispr.ita.rz.gen.dist %>% left_join(phylo.ita.rz.meta)

#plot a boxplot
ggplot(dispr.ita.rz.gen.plot,
       aes(x = factor(genotype,
                      level=c("Bulk soil",
                              "Rufi",
                              "Nipponbare",
                              "Kasalath",
                              "IR64")),
           y = dispr.ita.rz.gen.distances,
           fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 22) +
  scale_y_continuous(limits = c(8, 27),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.ita.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#nitrogen
set.seed(999)
dispr.ita.rz.nit <- vegan::betadisper(phylo.ita.rz.dist,
                                      phyloseq::sample_data(phylo.ita.rz.rclr)$nitro,
                                      type = "centroid")

#anova test
permutest(dispr.ita.rz.nit, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.ita.rz.nit.dist <- data.frame(dispr.ita.rz.nit$distances)

#combine data with metadata
dispr.ita.rz.nit.dist <- tibble::rownames_to_column(dispr.ita.rz.nit.dist, "sample.id")
dispr.ita.rz.nit.plot <- dispr.ita.rz.nit.dist %>% left_join(phylo.ita.rz.meta)

#plot a boxplot
ggplot(dispr.ita.rz.nit.plot,
       aes(x = factor(nitro,
                      level=c("N0", "N50")),
           y = dispr.ita.rz.nit.distances,
           fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 23) +
  scale_y_continuous(limits = c(8, 27),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.ita.rz.nit.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#Root
phylo.ita.rt <- subset_samples(ps1, soil == "Italian paddy soil" & 
                                 compartment %in% c("Root", "Bulk soil"))

phylo.ita.rt.rclr <- microbiome::transform(phylo.ita.rt, "rclr")

#extract metadata
phylo.ita.rt.meta <- meta(phylo.ita.rt)

#Generate distance matrix/calculate robust Aitchison distance
phylo.ita.rt.dist <- phyloseq::distance(phylo.ita.rt.rclr, method = "euclidean")

#ADONIS test
set.seed(999)
phylo.ita.rt.permanova <- vegan::adonis2(phylo.ita.rt.dist ~ 
                                           phyloseq::sample_data(phylo.ita.rt.rclr)$
                                           genotype*phyloseq::sample_data(phylo.ita.rt.rclr)$nitro,
                                         permutations = 9999)
phylo.ita.rt.permanova

#pairwise adonis
pairwise.adonis(x = phylo.ita.rt.dist,
                factors = sample_data(phylo.ita.rt.rclr)$genotype,
                p.adjust.m = "bonferroni",
                perm = 9999)

#calculate beta dispersion
#genotype
set.seed(999)
dispr.ita.rt.gen <- vegan::betadisper(phylo.ita.rt.dist,
                                      phyloseq::sample_data(phylo.ita.rt.rclr)$genotype,
                                      type = "centroid")

#anova test
permutest(dispr.ita.rt.gen, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.ita.rt.gen.dist <- data.frame(dispr.ita.rt.gen$distances)

#combine data with metadata
dispr.ita.rt.gen.dist <- tibble::rownames_to_column(dispr.ita.rt.gen.dist, "sample.id")
dispr.ita.rt.gen.plot <- dispr.ita.rt.gen.dist %>% left_join(phylo.ita.rt.meta)

#plot a boxplot
ggplot(dispr.ita.rt.gen.plot,
       aes(x = factor(genotype,
                      level=c("Bulk soil",
                              "Rufi",
                              "Nipponbare",
                              "Kasalath",
                              "IR64")),
           y = dispr.ita.rt.gen.distances,
           fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 20) +
  scale_y_continuous(limits = c(8, 30),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.ita.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#nitrogen
set.seed(999)
dispr.ita.rt.nit <- vegan::betadisper(phylo.ita.rt.dist,
                                      phyloseq::sample_data(phylo.ita.rt.rclr)$nitro,
                                      type = "centroid")

#anova test
permutest(dispr.ita.rt.nit, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.ita.rt.nit.dist <- data.frame(dispr.ita.rt.nit$distances)

#combine data with metadata
dispr.ita.rt.nit.dist <- tibble::rownames_to_column(dispr.ita.rt.nit.dist, "sample.id")
dispr.ita.rt.nit.plot <- dispr.ita.rt.nit.dist %>% left_join(phylo.ita.rt.meta)

#plot a boxplot
ggplot(dispr.ita.rt.nit.plot,
       aes(x = factor(nitro,
                      level=c("N0", "N50")),
           y = dispr.ita.rt.nit.distances,
           fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 21) +
  scale_y_continuous(limits = c(8, 30),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.ita.rt.nit.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)

#Philippine soil subset
#rhizosphere
phylo.php.rz <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                 compartment %in% c("Rhizosphere", "Bulk soil"))

phylo.php.rz.rclr <- microbiome::transform(phylo.php.rz, "rclr")

#extract metadata
phylo.php.rz.meta <- meta(phylo.php.rz)

#Generate distance matrix/calculate robust Aitchison distance
phylo.php.rz.dist <- phyloseq::distance(phylo.php.rz.rclr, method = "euclidean")

#ADONIS test
set.seed(999)
phylo.php.rz.permanova <- vegan::adonis2(phylo.php.rz.dist ~ 
                                           phyloseq::sample_data(phylo.php.rz.rclr)$
                                           genotype*phyloseq::sample_data(phylo.php.rz.rclr)$nitro,
                                         permutations = 9999)
phylo.php.rz.permanova

#pairwise adonis
pairwise.adonis(x = phylo.php.rz.dist,
                factors = sample_data(phylo.php.rz.rclr)$genotype,
                p.adjust.m = "bonferroni",
                perm = 9999)

#calculate beta dispersion
#genotype
set.seed(999)
dispr.php.rz.gen <- vegan::betadisper(phylo.php.rz.dist,
                                      phyloseq::sample_data(phylo.php.rz.rclr)$genotype,
                                      type = "centroid")

#anova test
permutest(dispr.php.rz.gen, permutations = 9999)
TukeyHSD(dispr.php.rz.gen)

#extract data for plotting
#distance to centroid
dispr.php.rz.gen.dist <- data.frame(dispr.php.rz.gen$distances)

#combine data with metadata
dispr.php.rz.gen.dist <- tibble::rownames_to_column(dispr.php.rz.gen.dist, "sample.id")
dispr.php.rz.gen.plot <- dispr.php.rz.gen.dist %>% left_join(phylo.php.rz.meta)

#plot a boxplot
ggplot(dispr.php.rz.gen.plot,
       aes(x = factor(genotype,
                      level=c("Bulk soil",
                              "Rufi",
                              "Nipponbare",
                              "Kasalath",
                              "IR64")),
           y = dispr.php.rz.gen.distances,
           fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "b", "a\nb", "a\nb"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 25) +
  scale_y_continuous(limits = c(8, 35),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.php.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#nitrogen
set.seed(999)
dispr.php.rz.nit <- vegan::betadisper(phylo.php.rz.dist,
                                      phyloseq::sample_data(phylo.php.rz.rclr)$nitro,
                                      type = "centroid")

#anova test
permutest(dispr.php.rz.nit, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.php.rz.nit.dist <- data.frame(dispr.php.rz.nit$distances)

#combine data with metadata
dispr.php.rz.nit.dist <- tibble::rownames_to_column(dispr.php.rz.nit.dist, "sample.id")
dispr.php.rz.nit.plot <- dispr.php.rz.nit.dist %>% left_join(phylo.php.rz.meta)

#plot a boxplot
ggplot(dispr.php.rz.nit.plot,
       aes(x = factor(nitro,
                      level=c("N0", "N50")),
           y = dispr.php.rz.nit.distances,
           fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 26) +
  scale_y_continuous(limits = c(8, 35),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.php.rz.nit.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#Root
phylo.php.rt <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                 compartment %in% c("Root", "Bulk soil"))

phylo.php.rt.rclr <- microbiome::transform(phylo.php.rt, "rclr")

#extract metadata
phylo.php.rt.meta <- meta(phylo.php.rt)

#Generate distance matrix/calculate robust Aitchison distance
phylo.php.rt.dist <- phyloseq::distance(phylo.php.rt.rclr, method = "euclidean")

#ADONIS test
set.seed(999)
phylo.php.rt.permanova <- vegan::adonis2(phylo.php.rt.dist ~ 
                                           phyloseq::sample_data(phylo.php.rt.rclr)$
                                           genotype*phyloseq::sample_data(phylo.php.rt.rclr)$nitro,
                                         permutations = 9999)
phylo.php.rt.permanova

#pairwise adonis
pairwise.adonis(x = phylo.php.rt.dist,
                factors = sample_data(phylo.php.rt.rclr)$genotype,
                p.adjust.m = "bonferroni",
                perm = 9999)

#calculate beta dispersion
#genotype
dispr.php.rt.gen <- vegan::betadisper(phylo.php.rt.dist,
                                      phyloseq::sample_data(phylo.php.rt.rclr)$genotype,
                                      type = "centroid",
                                      bias.adjust = TRUE)

#anova test
permutest(dispr.php.rt.gen, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.php.rt.gen.dist <- data.frame(dispr.php.rt.gen$distances)

#combine data with metadata
dispr.php.rt.gen.dist <- tibble::rownames_to_column(dispr.php.rt.gen.dist, "sample.id")
dispr.php.rt.gen.plot <- dispr.php.rt.gen.dist %>% left_join(phylo.php.rt.meta)

#plot a boxplot
ggplot(dispr.php.rt.gen.plot,
       aes(x = factor(genotype,
                      level=c("Bulk soil",
                              "Rufi",
                              "Nipponbare",
                              "Kasalath",
                              "IR64")),
           y = dispr.php.rt.gen.distances,
           fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 20) +
  scale_y_continuous(limits = c(8, 25),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.php.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)


#nitrogen
set.seed(999)
dispr.php.rt.nit <- vegan::betadisper(phylo.php.rt.dist,
                                      phyloseq::sample_data(phylo.php.rt.rclr)$nitro,
                                      type = "centroid")

#anova test
permutest(dispr.php.rt.nit, permutations = 9999)

#extract data for plotting
#distance to centroid
dispr.php.rt.nit.dist <- data.frame(dispr.php.rt.nit$distances)

#combine data with metadata
dispr.php.rt.nit.dist <- tibble::rownames_to_column(dispr.php.rt.nit.dist, "sample.id")
dispr.php.rt.nit.plot <- dispr.php.rt.nit.dist %>% left_join(phylo.php.rt.meta)

#plot a boxplot
ggplot(dispr.php.rt.nit.plot,
       aes(x = factor(nitro,
                      level=c("N0", "N50")),
           y = dispr.php.rt.nit.distances,
           fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 21) +
  scale_y_continuous(limits = c(8, 25),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Distance to centroid") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")

ggsave(plot = last_plot(),
       filename = "betadispr.php.rt.nit.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1282.67, unit = "px",
       dpi = 300)

