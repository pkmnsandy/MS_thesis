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

#Italian paddy soil
#rhizosphere
phylo.ita.rz <- subset_samples(ps1, soil == "Italian paddy soil" & 
                              compartment %in% c("Rhizosphere", "Bulk soil"))

#Aitchison distance
phylo.ita.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.ita.rz)

#transform counts  by rclr transformation
phylo.ita.rz.rclr <- transformCounts(phylo.ita.rz.tse, method = "rclr", abund_values = "counts")

#view clr table
phylo.ita.rz.rclr.assay <- assays(phylo.ita.rz.rclr)$rclr

#transpose clr table
phylo.ita.rz.rclr.assay <- t(phylo.ita.rz.rclr.assay)

#PCA calculation
phylo.ita.rz.pca <- prcomp(phylo.ita.rz.rclr.assay, center=T, scale = F)

#PC contributions
phylo.ita.rz.pca.sd <- phylo.ita.rz.pca$sdev
phylo.ita.rz.pca.loadings <- phylo.ita.rz.pca$rotation
phylo.ita.rz.pca.scores <- phylo.ita.rz.pca$x
phylo.ita.rz.pca.var <- phylo.ita.rz.pca.sd^2
phylo.ita.rz.pca.varPercent <- (phylo.ita.rz.pca.var/sum(phylo.ita.rz.pca.var))*100

phylo.ita.rz.pca.scree <- data.frame(phylo.ita.rz.pca.sd, phylo.ita.rz.pca.var, phylo.ita.rz.pca.varPercent)

#isolation of PCA data for plotting
phylo.ita.rz.pca.loadings <- phylo.ita.rz.pca$rotation # Loadings Matrix
phylo.ita.rz.pca.scores <- predict(phylo.ita.rz.pca) # Scores Matrix
phylo.ita.rz.pca.importance <- summary(phylo.ita.rz.pca)$importance # Explained Variance

#first two PC components.
phylo.ita.rz.pca.plotDataScores <- cbind(phylo.ita.rz.pca.scores[, c("PC1", "PC2")],
                                    colData(phylo.ita.rz.rclr)) %>%
  as.data.frame()

#genotype
phylo.ita.rz.pca.centroid.gen <- aggregate(cbind(PC1, PC2)~genotype, data = phylo.ita.rz.pca.plotDataScores, mean)
phylo.ita.rz.pca.gen <- merge(phylo.ita.rz.pca.plotDataScores,
                                        phylo.ita.rz.pca.centroid.gen,
                                           by = "genotype",
                                           suffixes = c("", ".centroid"))

ggplot(phylo.ita.rz.pca.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.ita.rz.pca.centroid.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.ita.rz.pca.gen,
             aes(fill = genotype),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (14.0% explained variance)") +
  ylab("PC2 (8.5% explained variance)") +
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
       filename = "beta diversity.ita.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1300, height = 1300, unit = "px",
       dpi = 300)

#nitrogen
#genotype
phylo.ita.rz.pca.centroid.nit <- aggregate(cbind(PC1, PC2)~nitro, data = phylo.ita.rz.pca.plotDataScores, mean)
phylo.ita.rz.pca.nit <- merge(phylo.ita.rz.pca.plotDataScores,
                              phylo.ita.rz.pca.centroid.nit,
                              by = "nitro",
                              suffixes = c("", ".centroid"))

ggplot(phylo.ita.rz.pca.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.ita.rz.pca.centroid.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.ita.rz.pca.nit,
             aes(fill = nitro),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (14.0% explained variance)") +
  ylab("PC2 (8.5% explained variance)") +
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
       filename = "beta diversity.ita.rz.nitrogen.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1300, height = 1300, unit = "px",
       dpi = 300)


#Root
phylo.ita.rt <- subset_samples(ps1, soil == "Italian paddy soil" & 
                                 compartment %in% c("Root", "Bulk soil"))

#Aitchison distance
phylo.ita.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.ita.rt)

#transform counts  by rclr transformation
phylo.ita.rt.rclr <- transformCounts(phylo.ita.rt.tse, method = "rclr", abund_values = "counts")

#view clr table
phylo.ita.rt.rclr.assay <- assays(phylo.ita.rt.rclr)$rclr

#transpose clr table
phylo.ita.rt.rclr.assay <- t(phylo.ita.rt.rclr.assay)

#PCA calculation
phylo.ita.rt.pca <- prcomp(phylo.ita.rt.rclr.assay, center=T, scale = F)

#PC contributions
phylo.ita.rt.pca.sd <- phylo.ita.rt.pca$sdev
phylo.ita.rt.pca.loadings <- phylo.ita.rt.pca$rotation
phylo.ita.rt.pca.scores <- phylo.ita.rt.pca$x
phylo.ita.rt.pca.var <- phylo.ita.rt.pca.sd^2
phylo.ita.rt.pca.varPercent <- (phylo.ita.rt.pca.var/sum(phylo.ita.rt.pca.var))*100

phylo.ita.rt.pca.scree <- data.frame(phylo.ita.rt.pca.sd, phylo.ita.rt.pca.var, phylo.ita.rt.pca.varPercent)

#isolation of PCA data for plotting
phylo.ita.rt.pca.loadings <- phylo.ita.rt.pca$rotation # Loadings Matrix
phylo.ita.rt.pca.scores <- predict(phylo.ita.rt.pca) # Scores Matrix
phylo.ita.rt.pca.importance <- summary(phylo.ita.rt.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
phylo.ita.rt.pca.plotDataScores <- cbind(phylo.ita.rt.pca.scores[, c("PC1", "PC2")],
                                         colData(phylo.ita.rt.rclr)) %>%
  as.data.frame()

#genotype
phylo.ita.rt.pca.centroid.gen <- aggregate(cbind(PC1, PC2)~genotype, data = phylo.ita.rt.pca.plotDataScores, mean)
phylo.ita.rt.pca.gen <- merge(phylo.ita.rt.pca.plotDataScores,
                              phylo.ita.rt.pca.centroid.gen,
                              by = "genotype",
                              suffixes = c("", ".centroid"))

ggplot(phylo.ita.rt.pca.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.ita.rt.pca.centroid.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.ita.rt.pca.gen,
             aes(fill = genotype),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (12.4% explained variance)") +
  ylab("PC2 (8.8% explained variance)") +
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
       filename = "beta diversity.ita.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1300, height = 1300, unit = "px",
       dpi = 300)

#nitrogen
phylo.ita.rt.pca.centroid.nit <- aggregate(cbind(PC1, PC2)~nitro, data = phylo.ita.rt.pca.plotDataScores, mean)
phylo.ita.rt.pca.nit <- merge(phylo.ita.rt.pca.plotDataScores,
                              phylo.ita.rt.pca.centroid.nit,
                              by = "nitro",
                              suffixes = c("", ".centroid"))

ggplot(phylo.ita.rt.pca.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.ita.rt.pca.centroid.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.ita.rt.pca.nit,
             aes(fill = nitro),size = 2, alpha = 1, color = "black", shape = 21) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  xlab("PC1 (12.4% explained variance)") +
  ylab("PC2 (8.8% explained variance)") +
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
       filename = "beta diversity.ita.rt.nitrogen.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1300, height = 1300, unit = "px",
       dpi = 300)


#Philippine paddy soil
#rhizosphere
phylo.php.rz <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                 compartment %in% c("Rhizosphere", "Bulk soil"))

#Aitchison distance
phylo.php.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.php.rz)

#transform counts  by rclr transformation
phylo.php.rz.rclr <- transformCounts(phylo.php.rz.tse, method = "rclr", abund_values = "counts")

#view clr table
phylo.php.rz.rclr.assay <- assays(phylo.php.rz.rclr)$rclr

#transpose clr table
phylo.php.rz.rclr.assay <- t(phylo.php.rz.rclr.assay)

#PCA calculation
phylo.php.rz.pca <- prcomp(phylo.php.rz.rclr.assay, center=T, scale = F)

#PC contributions
phylo.php.rz.pca.sd <- phylo.php.rz.pca$sdev
phylo.php.rz.pca.loadings <- phylo.php.rz.pca$rotation
phylo.php.rz.pca.scores <- phylo.php.rz.pca$x
phylo.php.rz.pca.var <- phylo.php.rz.pca.sd^2
phylo.php.rz.pca.varPercent <- (phylo.php.rz.pca.var/sum(phylo.php.rz.pca.var))*100

phylo.php.rz.pca.scree <- data.frame(phylo.php.rz.pca.sd, phylo.php.rz.pca.var, phylo.php.rz.pca.varPercent)

#isolation of PCA data for plotting
phylo.php.rz.pca.loadings <- phylo.php.rz.pca$rotation # Loadings Matrix
phylo.php.rz.pca.scores <- predict(phylo.php.rz.pca) # Scores Matrix
phylo.php.rz.pca.importance <- summary(phylo.php.rz.pca)$importance # Explained Variance

#first two PC components.
phylo.php.rz.pca.plotDataScores <- cbind(phylo.php.rz.pca.scores[, c("PC1", "PC2")],
                                         colData(phylo.php.rz.rclr)) %>%
  as.data.frame()

#genotype
phylo.php.rz.pca.centroid.gen <- aggregate(cbind(PC1, PC2)~genotype, data = phylo.php.rz.pca.plotDataScores, mean)
phylo.php.rz.pca.gen <- merge(phylo.php.rz.pca.plotDataScores,
                              phylo.php.rz.pca.centroid.gen,
                              by = "genotype",
                              suffixes = c("", ".centroid"))

ggplot(phylo.php.rz.pca.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.php.rz.pca.centroid.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.php.rz.pca.gen,
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
       width = 1300, height = 1300, unit = "px",
       dpi = 300)

#nitrogen
phylo.php.rz.pca.centroid.nit <- aggregate(cbind(PC1, PC2)~nitro, data = phylo.php.rz.pca.plotDataScores, mean)
phylo.php.rz.pca.nit <- merge(phylo.php.rz.pca.plotDataScores,
                              phylo.php.rz.pca.centroid.nit,
                              by = "nitro",
                              suffixes = c("", ".centroid"))

ggplot(phylo.php.rz.pca.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.php.rz.pca.centroid.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.php.rz.pca.nit,
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
       width = 1300, height = 1300, unit = "px",
       dpi = 300)


#Root
phylo.php.rt <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                 compartment %in% c("Root", "Bulk soil"))

#Aitchison distance
phylo.php.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.php.rt)

#transform counts  by rclr transformation
phylo.php.rt.rclr <- transformCounts(phylo.php.rt.tse, method = "rclr", abund_values = "counts")

#view clr table
phylo.php.rt.rclr.assay <- assays(phylo.php.rt.rclr)$rclr

#transpose clr table
phylo.php.rt.rclr.assay <- t(phylo.php.rt.rclr.assay)

#PCA calculation
phylo.php.rt.pca <- prcomp(phylo.php.rt.rclr.assay, center=T, scale = F)

#PC contributions
phylo.php.rt.pca.sd <- phylo.php.rt.pca$sdev
phylo.php.rt.pca.loadings <- phylo.php.rt.pca$rotation
phylo.php.rt.pca.scores <- phylo.php.rt.pca$x
phylo.php.rt.pca.var <- phylo.php.rt.pca.sd^2
phylo.php.rt.pca.varPercent <- (phylo.php.rt.pca.var/sum(phylo.php.rt.pca.var))*100

phylo.php.rt.pca.scree <- data.frame(phylo.php.rt.pca.sd, phylo.php.rt.pca.var, phylo.php.rt.pca.varPercent)

#isolation of PCA data for plotting
phylo.php.rt.pca.loadings <- phylo.php.rt.pca$rotation # Loadings Matrix
phylo.php.rt.pca.scores <- predict(phylo.php.rt.pca) # Scores Matrix
phylo.php.rt.pca.importance <- summary(phylo.php.rt.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
phylo.php.rt.pca.plotDataScores <- cbind(phylo.php.rt.pca.scores[, c("PC1", "PC2")],
                                         colData(phylo.php.rt.rclr)) %>%
  as.data.frame()

#genotype
phylo.php.rt.pca.centroid.gen <- aggregate(cbind(PC1, PC2)~genotype, data = phylo.php.rt.pca.plotDataScores, mean)
phylo.php.rt.pca.gen <- merge(phylo.php.rt.pca.plotDataScores,
                              phylo.php.rt.pca.centroid.gen,
                              by = "genotype",
                              suffixes = c("", ".centroid"))

ggplot(phylo.php.rt.pca.gen, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.php.rt.pca.centroid.gen, aes(x = PC1, y = PC2, fill = genotype),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.php.rt.pca.gen,
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
       width = 1300, height = 1300, unit = "px",
       dpi = 300)

#nitrogen
phylo.php.rt.pca.centroid.nit <- aggregate(cbind(PC1, PC2)~nitro, data = phylo.php.rt.pca.plotDataScores, mean)
phylo.php.rt.pca.nit <- merge(phylo.php.rt.pca.plotDataScores,
                              phylo.php.rt.pca.centroid.nit,
                              by = "nitro",
                              suffixes = c("", ".centroid"))

ggplot(phylo.php.rt.pca.nit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = PC1.centroid, y = PC2.centroid, xend = PC1, yend = PC2),
               color = "black", alpha = 0.3) +
  geom_point(data = phylo.php.rt.pca.centroid.nit, aes(x = PC1, y = PC2, fill = nitro),
             size = 4, color = "black", shape = 21) +
  geom_point(data = phylo.php.rt.pca.nit,
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
       width = 1300, height = 1300, unit = "px",
       dpi = 300)