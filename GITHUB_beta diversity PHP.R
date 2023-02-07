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

#Philippine soil subset
phylo.php <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                              compartment %in% c("Rhizosphere", "Bulk soil"))

#Aitchison distance
rhizosphere.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.php)

#transform counts  by clr transformation
rhizosphere.tse.clr <- transformCounts(rhizosphere.tse, method = "rclr", abund_values = "counts")

#view clr table
rhizosphere.clr_assay <- assays(rhizosphere.tse.clr)$rclr

#transpose clr table
rhizosphere.clr_assay <- t(rhizosphere.clr_assay)

#calcualte Euclidean distance between samples
rhizosphere.euclidean_dist <- vegan::vegdist(rhizosphere.clr_assay, method="euclidean")

#PERMANOVA
#getting relative abundance
rhizosphere.tax_abundance <- assays(rhizosphere.tse.clr)$rclr
rhizosphere.tax_abundance <- t(rhizosphere.tax_abundance)

rhizosphere.pca <- prcomp(rhizosphere.clr_assay, center=T, scale = F)

#PCA result for variables
var.res.rhizosphere <- get_pca(rhizosphere.pca, "var")

#screeplot
rhizosphere.sd <- rhizosphere.pca$sdev
rhizosphere.loadings <- rhizosphere.pca$rotation
rownames(rhizosphere.loadings) <- colnames(rhizosphere.euclidean_dist)
rhizosphere.scores <- rhizosphere.pca$x
rhizosphere.var <- rhizosphere.sd^2
rhizosphere.varPercent <- (rhizosphere.var/sum(rhizosphere.var))*100

rhizosphere.pca.scree <- data.frame(rhizosphere.sd, rhizosphere.var, rhizosphere.varPercent)

#explained variance
sum(rhizosphere.varPercent[1:14]) #60.79858% contribution of the first five PCs

rhizosphere.loadings <- rhizosphere.pca$rotation # Loadings Matrix
rhizosphere.scores <- predict(rhizosphere.pca) # Scores Matrix
rhizosphere.importance <- summary(rhizosphere.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
rhizosphere.plotDataScores <- cbind(rhizosphere.scores[, c("PC1", "PC2")],
                                    colData(rhizosphere.tse.clr)) %>%
  as.data.frame()

#genotype
rhizosphere.plotDataScores.gen <- rhizosphere.plotDataScores
centroid.php.rz.gen <- aggregate(cbind(PC1, PC2)~genotype, data = rhizosphere.plotDataScores.gen, mean)
rhizosphere.plotDataScores.gen <- merge(rhizosphere.plotDataScores.gen,
                                        centroid.php.rz.gen,
                                        by = "genotype",
                                        suffixes = c("", ".centroid"))

ggplot(rhizosphere.plotDataScores.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = centroid.php.rz.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = rhizosphere.plotDataScores.gen,
             aes(fill = genotype),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (13.7% explained variance)") +
  ylab("PC2 (6.5% explained variance)") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = c(0.8, 0.78))

ggsave(plot = last_plot(),
       filename = "beta diversity.php.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#nitrogen
rhizosphere.plotDataScores.nit <- rhizosphere.plotDataScores
centroid.php.rz.nit <- aggregate(cbind(PC1, PC2)~nitro, data = rhizosphere.plotDataScores.nit, mean)
rhizosphere.plotDataScores.nit <- merge(rhizosphere.plotDataScores.nit,
                                        centroid.php.rz.nit,
                                        by = "nitro",
                                        suffixes = c("", ".centroid"))

ggplot(rhizosphere.plotDataScores.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = centroid.php.rz.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = rhizosphere.plotDataScores.nit,
             aes(fill = nitro),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (13.7% explained variance)") +
  ylab("PC2 (6.5% explained variance)") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = c(0.8, 0.78))

ggsave(plot = last_plot(),
       filename = "beta diversity.php.rz.nitrogen.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#Root compartment
phylo.php <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                              compartment %in% c("Root", "Bulk soil"))

#Aitchison distance
Root.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.php)

#transform counts  by clr transformation
Root.tse.clr <- transformCounts(Root.tse, method = "rclr", abund_values = "counts")

#view clr table
Root.clr_assay <- assays(Root.tse.clr)$rclr

#transpose clr table
Root.clr_assay <- t(Root.clr_assay)

#calcualte Euclidean distance between samples
Root.euclidean_dist <- vegan::vegdist(Root.clr_assay, method="euclidean")

#PERMANOVA
#getting relative abundance
Root.tax_abundance <- assays(Root.tse.clr)$rclr
Root.tax_abundance <- t(Root.tax_abundance)

Root.pca <- prcomp(Root.clr_assay, center=T, scale = F)

#PCA result for variables
var.res.Root <- get_pca(Root.pca, "var")

#screeplot
Root.sd <- Root.pca$sdev
Root.loadings <- Root.pca$rotation
rownames(Root.loadings) <- colnames(Root.euclidean_dist)
Root.scores <- Root.pca$x
Root.var <- Root.sd^2
Root.varPercent <- (Root.var/sum(Root.var))*100

Root.pca.scree <- data.frame(Root.sd, Root.var, Root.varPercent)

#explained variance
sum(Root.varPercent[1:14]) #56.08445% contribution of the first five PCs

Root.loadings <- Root.pca$rotation # Loadings Matrix
Root.scores <- predict(Root.pca) # Scores Matrix
Root.importance <- summary(Root.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
Root.plotDataScores <- cbind(Root.scores[, c("PC1", "PC2")],
                             colData(Root.tse.clr)) %>%
  as.data.frame()

#genotype
Root.plotDataScores.gen <- Root.plotDataScores
centroid.php.rt.gen <- aggregate(cbind(PC1, PC2)~genotype, data = Root.plotDataScores.gen, mean)
Root.plotDataScores.gen <- merge(Root.plotDataScores.gen,
                                 centroid.php.rt.gen,
                                 by = "genotype",
                                 suffixes = c("", ".centroid"))

ggplot(Root.plotDataScores.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = centroid.php.rt.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = Root.plotDataScores.gen,
             aes(fill = genotype),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (9.3% explained variance)") +
  ylab("PC2 (6.4% explained variance)") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = c(0.8, 0.78))

ggsave(plot = last_plot(),
       filename = "beta diversity.php.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#nitrogen
Root.plotDataScores.nit <- Root.plotDataScores
centroid.php.rt.nit <- aggregate(cbind(PC1, PC2)~nitro, data = Root.plotDataScores.nit, mean)
Root.plotDataScores.nit <- merge(Root.plotDataScores.nit,
                                 centroid.php.rt.nit,
                                 by = "nitro",
                                 suffixes = c("", ".centroid"))

ggplot(Root.plotDataScores.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = centroid.php.rt.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = Root.plotDataScores.nit,
             aes(fill = nitro),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (9.3% explained variance)") +
  ylab("PC2 (6.4% explained variance)") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = c(0.8, 0.78))

ggsave(plot = last_plot(),
       filename = "beta diversity.php.rt.nitrogen.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)
