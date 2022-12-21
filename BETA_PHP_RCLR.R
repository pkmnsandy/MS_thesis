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

ggplot(rhizosphere.pca.scree,
       aes(x = 1:nrow(rhizosphere.pca.scree), y=rhizosphere.varPercent)) +
  geom_col() +
  geom_hline(yintercept = (1/nrow(rhizosphere.pca.scree))*100,
             linetype = "solid", color = "red", size = 0.5) + 
  xlab("Principal components") +
  ylab("Variance (%)") +
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
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "PCA.php.rhizosphere.scree.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)


#explained variance
sum(rhizosphere.varPercent[1:13]) #56.4159% contribution of the first five PCs

rhizosphere.loadings <- rhizosphere.pca$rotation # Loadings Matrix
rhizosphere.scores <- predict(rhizosphere.pca) # Scores Matris
rhizosphere.importance <- summary(rhizosphere.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
rhizosphere.plotDataScores <- cbind(rhizosphere.scores[, c("PC1", "PC2")],
                                    colData(rhizosphere.tse.clr)) %>%
  as.data.frame()

# We want to use min-max feature scaling on the scores so they are between zero and one, the same as the loadings.
normalize <- function(x) return ((x - min(x)) / (max(x) - min(x)))

rhizosphere.plotDataScores[, "PC1"] <- scale(normalize(rhizosphere.plotDataScores[, "PC1"]), center = TRUE, scale = FALSE)
rhizosphere.plotDataScores[, "PC2"] <- scale(normalize(rhizosphere.plotDataScores[, "PC2"]), center = TRUE, scale = FALSE)

rhizosphere.plotDataLoadings <- as.data.frame(rhizosphere.loadings) %>%
  filter(row.names(rhizosphere.loadings) %in% c("ASV1642", "ASV1621", "ASV1910", "ASV1364", "ASV1331",
                                                "ASV1632", "ASV515", "ASV2322", "ASV2320", "ASV2127"
                                                )) %>%
  select(c(PC1, PC2))

ggplot(rhizosphere.plotDataScores, aes(x = PC1[,1], y = PC2[,1])) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(data = rhizosphere.plotDataScores, aes(fill = genotype, shape = nitro), size = 3, alpha = 1, color = "black") +
  geom_segment(data = rhizosphere.plotDataLoadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.03, "npc")), alpha = 0.5) + 
  geom_text(data = rhizosphere.plotDataLoadings,
            mapping = aes(x = PC1, y = PC2,label = rownames(rhizosphere.plotDataLoadings)),
            hjust = 1, vjust = -0.2, colour = "darkred", size = 4, check_overlap = FALSE) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  scale_shape_manual(values = c(21, 24)) +
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
        legend.position = c(0.2, 0.3))

ggsave(plot = last_plot(),
       filename = "PCA.php.rhizosphere.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 1500, height = 2000, unit = "px",
       dpi = 300)

var.dat <- fviz_contrib(rhizosphere.pca, choice = "var", axes = 1:2, top = 25) +
  xlab("Amplicon sequence variants") +
  scale_x_discrete(labels = c(
    "ASV1642: RPC 2a",
    "ASV1621: RPC 2d",
    "ASV1910: Methylococcus",
    "ASV1364: Methylosarcina",
    "ASV1331: Methylomicrobium",
    "ASV1632: RPC 2a",
    "ASV515: Unclassified Methylocystaceae",
    "ASV2322: Methylococcus",
    "ASV2320: Methylococcus",
    "ASV2127: Methylococcus",
    "ASV2772: Methylocystis",
    "ASV1325: Methylosarcina, uncultured",
    "ASV1992: Methylococcus",
    "ASV53: Methylocystis",
    "ASV2139: Methylococcus",
    "ASV2459: Methylocystis",
    "ASV2135: Methylococcus",
    "ASV1633: RPC 2a",
    "ASV877: Methylocystis",
    "ASV1: Methylocystis",
    "ASV2129: Methylococcus",
    "ASV1678: RPC 2a",
    "ASV2159: Methylococcus",
    "ASV2134: Methylococcus",
    "ASV301: Methylocystis"
  )) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        plot.title = element_blank(),
        axis.title.x = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10, face = "italic",
                                   angle = 60, vjust = 1, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = "right")

var.dat

ggsave(plot = last_plot(),
       filename = "PCA.php.rhizosphere.contrib.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1500, unit = "px",
       dpi = 300)

pg <- ggplot_build(var.dat)
pg.df <- pg$plot$data


######################################################################################################################
######################################################################################################################
######################################################################################################################

#Root subset
phylo.php.root <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                   compartment %in% c("Root", "Bulk soil"))

#Aitchison distance
Root.tse <- makeTreeSummarizedExperimentFromPhyloseq(phylo.php.root)

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


root.pca <- prcomp(Root.clr_assay, center=T, scale = F)

#PCA result for variables
var.res.root <- get_pca(root.pca, "var")

#screeplot
root.sd <- root.pca$sdev
root.loadings <- root.pca$rotation
rownames(root.loadings) <- colnames(Root.euclidean_dist)
root.scores <- root.pca$x
root.var <- root.sd^2
root.varPercent <- (root.var/sum(root.var))*100

root.pca.scree <- data.frame(root.sd, root.var, root.varPercent)

ggplot(root.pca.scree,
       aes(x = 1:nrow(root.pca.scree), y=root.varPercent)) +
  geom_col() +
  geom_hline(yintercept = (1/nrow(root.pca.scree))*100,
             linetype = "solid", color = "red", size = 0.5) + 
  xlab("Principal components") +
  ylab("Variance (%)") +
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
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "PCA.php.root.scree.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1000, unit = "px",
       dpi = 300)


#explained variance
sum(root.varPercent[1:17]) #56.71781% contribution of the first five PCs

#exract the principal components
Root.loadings <- Root.pca$rotation # Loadings Matrix
Root.scores <- predict(Root.pca) # Scores Matris
Root.importance <- summary(Root.pca)$importance # Explained Variance

# The plot data includes occupation and the first two PC components.
Root.plotDataScores <- cbind(Root.scores[, c("PC1", "PC2")],
                             colData(Root.tse.clr)) %>%
  as.data.frame()

# We want to use min-max feature scaling on the scores so they are between zero and one, the same as the loadings.
normalize <- function(x) return ((x - min(x)) / (max(x) - min(x)))

Root.plotDataScores[, "PC1"] <- scale(normalize(Root.plotDataScores[, "PC1"]), center = TRUE, scale = FALSE)
Root.plotDataScores[, "PC2"] <- scale(normalize(Root.plotDataScores[, "PC2"]), center = TRUE, scale = FALSE)

fviz_contrib(root.pca, choice = "var", axes = 1:2, top = 25) +
  xlab("Amplicon sequence variants")


Root.plotDataLoadings <- as.data.frame(Root.loadings) %>%
  filter(row.names(Root.loadings) %in% c("ASV1642", "ASV1621", "ASV2127", "ASV2320","ASV2139",
                                         "ASV2322", "ASV1364", "ASV2134", "ASV2113", "ASV618")) %>%
  select(c(PC1, PC2))

ggplot(Root.plotDataScores, aes(x = PC1[,1], y = PC2[,1])) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(data = Root.plotDataScores, aes(fill = genotype, shape = nitro), size = 3, alpha = 1, color = "black") +
  geom_segment(data = Root.plotDataLoadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.03, "npc")), alpha = 0.5) + 
  geom_text(data = Root.plotDataLoadings,
            mapping = aes(x = PC1, y = PC2,label = rownames(Root.plotDataLoadings)),
            hjust = 1, vjust = -0.2, colour = "darkred", size = 4, check_overlap = FALSE) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  scale_shape_manual(values = c(21, 24)) +
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
       filename = "PCA.php.root.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 1500, height = 2000, unit = "px",
       dpi = 300)

fviz_contrib(root.pca, choice = "var", axes = 1:2, top = 25) +
  xlab("Amplicon sequence variants") +
  scale_x_discrete(labels = c(
    "ASV1642: RPC 2a",
    "ASV1621: RPC 2d",
    "ASV2127: Methylococcus",
    "ASV2320: Methylococcus",
    "ASV2139: Methylococcus",
    "ASV2322: Methylococcus",
    "ASV1364: Methylosarcina, uncultured",
    "ASV2134: Methylococcus",
    "ASV2113: Methylococcus",
    "ASV618: Methylocystis",
    "ASV2129: Methylococcus",
    "ASV2120: Methylococcus",
    "ASV2528: Methylocystis",
    "ASV2132: Methylococcus",
    "ASV647: Methylocystis",
    "ASV53: Methylocystis",
    "ASV2135: Methylococcus",
    "ASV760: Unclassified Methylocystis echinoides",
    "ASV1442: Methylobacter",
    "ASV1632: RPC 2a",
    "ASV2524: Methylocystis",
    "ASV682: Methylocystis",
    "ASV2159: Methylococcus",
    "ASV302: Methylocystis",
    "ASV1505: Methylosoma, uncultured")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        plot.title = element_blank(),
        axis.title.x = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 10, face = "italic",
                                   angle = 68, vjust = 1, hjust=1),
        axis.text.y = element_text(color="black", size = 10),
        legend.box.background = element_rect(linetype="solid", size = 0.5, color="white"),
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "PCA.php.root.contrib.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1500, unit = "px",
       dpi = 300)

pg <- ggplot_build(var.dat)
pg.df <- pg$plot$data
