library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling
library(ape)
library(picante)
library(car)
library(scales)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)
library(lmerTest)
library(tidyverse)
library(cowplot)
library(ggstance)
library(ggpubr)
library(sjPlot)
library(sjmisc)
library(xtable)
library(htmltools)
library(htmlTable)
library(lmerTest)
library(rstatix)
library(broom)
library(rcompanion)
library(FSA)

#making a phyloseq object

#taxonomy file
biom_data <- import_biom(BIOMfilename = "input/taxa_table.biom", treefilename = "input/tree.nwk")

#mapping file
mapping_file <- import_qiime_sample_data(mapfilename = "input/metadata.txt")

#reading tree file
tree_file <- read.tree("input/tree.nwk")

# Merge the OTU and mapping data into a phyloseq object
phylo <- merge_phyloseq(biom_data, mapping_file)

#Add names to biom table and check phyloseq objects
colnames(tax_table(phylo))= c("Kingdom","Phylum","Class","Order","Family", "Cluster", "Genus", "Species")
rank_names(phylo)

print(phylo)

#cleaning the taxonomy table
tax_table(phylo)[, colnames(tax_table(phylo))] <- gsub(tax_table(phylo)
                                                       [, colnames(tax_table(phylo))], pattern = ".*__", replacement = "")
tax_table(phylo)[tax_table(phylo)[, "Species"] == "", "Species"] <- "Unidentified"

datatable(tax_table(phylo))

#general examination
phylo.general.asvtab <- as.data.frame(phylo@otu_table)

phylo.general.tree <- phylo@phy_tree

# phylo.general.meta from previous code chunks

# We first need to check if the tree is rooted or not 

phylo@phy_tree
# it is an unrooted tree

phylo.general.pd <- pd(t(phylo.general.asvtab), phylo.general.tree, include.root=F)
datatable(phylo.general.pd)

#isolate the metadata file from the phyloseq object
phylo.general.meta <- meta(phylo)

#custom sort levels
phylo.general.meta$soil <- factor(phylo.general.meta$soil, levels=c('Italian paddy soil', 'Philippine paddy soil'))
phylo.general.meta$compartment <- factor(phylo.general.meta$compartment, levels=c('Bulk soil', 'Rhizosphere', 'Root'))

#add PD to the metadata file
phylo.general.meta$Phylogenetic_Diversity <- phylo.general.pd$PD

#Hypothesis testing
set.seed(999)

#ANOVA procedure on original data
phylo.general.aov.org <- aov(
  (Phylogenetic_Diversity) ~ soil*compartment, data = phylo.general.meta,
  contrasts = list(
    soil = 'contr.sum',
    compartment = 'contr.sum'
  )
)

#rank transformed
phylo.general.aov <- aov(
  rank(Phylogenetic_Diversity) ~ soil*compartment, data = phylo.general.meta,
  contrasts = list(
    soil = 'contr.sum',
    compartment = 'contr.sum'
  )
)

Anova(phylo.general.aov.org, type = 'III')
Anova(phylo.general.aov, type = 'III')


#means
model.tables(phylo.general.aov.org, type = "means", se = TRUE)

#Residuals
phylo.general.aov.res = phylo.general.aov$resid
#phylo.general.aov.res = phylo.general.aov.org$resid

#QQplot
qqnorm(
  phylo.general.aov.res, pch = 20, main = "Rank-transformed data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(phylo.general.aov.res)

#test for normality
shapiro.test(phylo.general.aov.res)

#homogeneity of variance
plot(phylo.general.aov.org, 1, main = "Rank-transformed data")
#plot(phylo.general.aov.org, 1, main = "Original data")

#test for homogeneity of variance
leveneTest(phylo.general.aov)
#leveneTest(phylo.general.aov.org)

#Non-parametric test
#Scheirer-Ray-Hare test
set.seed(999)
scheirerRayHare(Phylogenetic_Diversity ~ soil + compartment,
                data = phylo.general.meta)

ph <- phylo.general.meta %>%
  subset(soil == "Philippine paddy soil",
         select = c("Phylogenetic_Diversity", "compartment"))

#posthoc
phylo.general.ph.post <- dunnTest(Phylogenetic_Diversity ~ compartment,
                                  data = ph,
                                  method = "bonferroni")
phylo.general.ph.post

phylo.general.ph <- phylo.general.ph.post$res
cldList(P.adj ~ Comparison,
        data = phylo.general.ph,
        threshold = 0.05)


ita <- phylo.general.meta %>%
  subset(soil == "Italian paddy soil",
         select = c("Phylogenetic_Diversity", "compartment"))

#posthoc
phylo.general.ita.post <- dunnTest(Phylogenetic_Diversity ~ compartment,
                                   data = ita,
                                   method = "bonferroni")
phylo.general.ita.post

phylo.general.ita <- phylo.general.ita.post$res
cldList(P.adj ~ Comparison,
        data = phylo.general.ita,
        threshold = 0.05)

########################
########################
########################
#Italian paddy soil
#rhizosphere compartment
phylo.ita.rz <- subset_samples(phylo, soil=="Italian paddy soil" &
                                 compartment %in% c("Bulk soil", "Rhizosphere"))

phylo.ita.rz.asvtab <- as.data.frame(phylo.ita.rz@otu_table)

phylo.ita.rz.tree <- phylo.ita.rz@phy_tree

# phylo.ita.rz.meta from previous code chunks

# We first need to check if the tree is rooted or not 

phylo.ita.rz@phy_tree
# it is an unrooted tree

phylo.ita.rz.pd <- pd(t(phylo.ita.rz.asvtab), phylo.ita.rz.tree, include.root=F)
datatable(phylo.ita.rz.pd)

#isolate the metadata file from the phyloseq object
phylo.ita.rz.meta <- meta(phylo.ita.rz)

#custom sort levels
phylo.ita.rz.meta$genotype <- factor(phylo.ita.rz.meta$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'Kasalath', 'IR64'))
phylo.ita.rz.meta$nitro <- factor(phylo.ita.rz.meta$nitro, levels=c('N0', 'N50'))

#add PD to the metadata file
phylo.ita.rz.meta$Phylogenetic_Diversity <- phylo.ita.rz.pd$PD

#Hypothesis testing
set.seed(999)

#ANOVA procedure on original data
phylo.ita.rz.aov.org <- aov(
  (Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.ita.rz.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

#rank transformed
phylo.ita.rz.aov <- aov(
  rank(Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.ita.rz.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

Anova(phylo.ita.rz.aov, type = 'II')
Anova(phylo.ita.rz.aov.org, type = 'II')

#means
model.tables(phylo.ita.rz.aov, type = "means", se = TRUE)

#Residuals
phylo.ita.rz.aov.res = phylo.ita.rz.aov$resid
#phylo.ita.rz.aov.res = phylo.ita.rz.aov.org$resid

#QQplot
qqnorm(
  phylo.ita.rz.aov.res, pch = 20, main = "Rank-transformed data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(phylo.ita.rz.aov.res)

#test for normality
shapiro.test(phylo.ita.rz.aov.res)

#homogeneity of variance
plot(phylo.ita.rz.aov.org, 1, main = "Rank-transformed data")
#plot(phylo.ita.rz.aov.org, 1, main = "Original data")

#test for homogeneity of variance
leveneTest(phylo.ita.rz.aov)
#leveneTest(phylo.ita.rz.aov.org)

#Compute estimated marginal means for factor combinations
emmeans(phylo.ita.rz.aov, pairwise ~ genotype | nitro)
#emmeans(phylo.ita.rz.aov.org, pairwise ~ genotype | nitro)

phylo.ita.rz.aov.emmeans <- emmeans(phylo.ita.rz.aov,  ~ genotype | nitro)
#phylo.ita.rz.aov.emmeans <- emmeans(phylo.ita.rz.aov.org,  ~ genotype | nitro)
print(phylo.ita.rz.aov.emmeans)

phylo.ita.rz.aov.emmeans %>% 
  pairs() %>% 
  test(joint = TRUE)

pairs(phylo.ita.rz.aov.emmeans)

#multiple comparisons
summary(glht(phylo.ita.rz.aov, linfct = mcp(genotype = "Tukey")))
#summary(glht(phylo.ita.rz.aov.org, linfct = mcp(genotype = "Tukey")))

#visualization
phylo.ita.rz.lm <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = phylo.ita.rz.meta)

#genotype
phylo.ita.rz.aov.gen <- cld(emmeans(phylo.ita.rz.lm, c("genotype")), alpha = 0.05, Letters = "letters")

phylo.ita.rz.gen.mean<- phylo.ita.rz.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.ita.rz.aov.gen) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

phylo.ita.rz.gen.jitter <- phylo.ita.rz.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.ita.rz.meta) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

ggplot(phylo.ita.rz.meta, aes(x = factor(genotype,
                                         level=c("Bulk soil",
                                                 "Rufi",
                                                 "Nipponbare",
                                                 "Kasalath",
                                                 "IR64")),
                              y = Phylogenetic_Diversity,
                              fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 7.2) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.ita.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#nitrogen
ggplot(phylo.ita.rz.meta, aes(x = factor(nitro,
                                         level=c("N0", "N50")),
                              y = Phylogenetic_Diversity,
                              fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 7.2) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.ita.rz.nitro.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#root compartment
#Italian paddy soil
#Root compartment
phylo.ita.rt <- subset_samples(phylo, soil=="Italian paddy soil" &
                                 compartment %in% c("Bulk soil", "Root"))

phylo.ita.rt.asvtab <- as.data.frame(phylo.ita.rt@otu_table)

phylo.ita.rt.tree <- phylo.ita.rt@phy_tree

# phylo.ita.rt.meta from previous code chunks

# We first need to check if the tree is rooted or not 

phylo.ita.rt@phy_tree
# it is an unrooted tree

phylo.ita.rt.pd <- pd(t(phylo.ita.rt.asvtab), phylo.ita.rt.tree, include.root=F)
datatable(phylo.ita.rt.pd)

#isolate the metadata file from the phyloseq object
phylo.ita.rt.meta <- meta(phylo.ita.rt)

#custom sort levels
phylo.ita.rt.meta$genotype <- factor(phylo.ita.rt.meta$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'Kasalath', 'IR64'))
phylo.ita.rt.meta$nitro <- factor(phylo.ita.rt.meta$nitro, levels=c('N0', 'N50'))

#add PD to the metadata file
phylo.ita.rt.meta$Phylogenetic_Diversity <- phylo.ita.rt.pd$PD

#Hypothesis testing
set.seed(999)

#ANOVA procedure on original data
phylo.ita.rt.aov.org <- aov(
  (Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.ita.rt.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

#rank transformed
phylo.ita.rt.aov <- aov(
  rank(Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.ita.rt.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

set.seed(999)
Anova(phylo.ita.rt.aov, type = 'II')
Anova(phylo.ita.rt.aov.org, type = 'II')

#means
#model.tables(phylo.ita.rt.aov, type = "means", se = TRUE)
model.tables(phylo.ita.rt.aov.org, type = "means", se = TRUE)

#Residuals
phylo.ita.rt.aov.res = phylo.ita.rt.aov$resid
#phylo.ita.rt.aov.res = phylo.ita.rt.aov.org$resid

#QQplot
qqnorm(
  phylo.ita.rt.aov.res, pch = 20, main = "Rank-transformed data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(phylo.ita.rt.aov.res)

#test for normality
shapiro.test(phylo.ita.rt.aov.res)

#homogeneity of variance
plot(phylo.ita.rt.aov, 1, main = "Rank-transformed data")
#plot(phylo.ita.rt.aov.org, 1, main = "Original data")

#test for homogeneity of variance
leveneTest(phylo.ita.rt.aov)
#leveneTest(phylo.ita.rt.aov.org)

#Compute estimated marginal means for factor combinations
emmeans(phylo.ita.rt.aov, pairwise ~ genotype | nitro)
##emmeans(phylo.ita.rt.aov.org, pairwise ~ genotype | nitro)

phylo.ita.rt.aov.emmeans <- emmeans(phylo.ita.rt.aov,  ~ genotype | nitro)
#phylo.ita.rt.aov.emmeans <- emmeans(phylo.ita.rt.aov.org,  ~ genotype | nitro)
print(phylo.ita.rt.aov.emmeans)

phylo.ita.rt.aov.emmeans %>% 
  pairs() %>% 
  test(joint = TRUE)

pairs(phylo.ita.rt.aov.emmeans)

#multiple comparisons
summary(glht(phylo.ita.rt.aov, linfct = mcp(genotype = "Tukey")))
#summary(glht(phylo.ita.rt.aov.org, linfct = mcp(genotype = "Tukey")))

#visualization
phylo.ita.rt.lm <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = phylo.ita.rt.meta)

#genotype
phylo.ita.rt.aov.gen <- cld(emmeans(phylo.ita.rt.lm, c("genotype")), alpha = 0.05, Letters = "letters")

phylo.ita.rt.gen.mean<- phylo.ita.rt.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.ita.rt.aov.gen) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

phylo.ita.rt.gen.jitter <- phylo.ita.rt.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.ita.rt.meta) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

ggplot(phylo.ita.rt.meta, aes(x = factor(genotype,
                                         level=c("Bulk soil",
                                                 "Rufi",
                                                 "Nipponbare",
                                                 "Kasalath",
                                                 "IR64")),
                              y = Phylogenetic_Diversity,
                              fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 6.5) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.ita.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#nitrogen
#genotype
phylo.ita.rt.aov.nitro <- cld(emmeans(phylo.ita.rt.lm, c("nitro")), alpha = 0.05, Letters = "letters")

ggplot(phylo.ita.rt.meta, aes(x = factor(nitro,
                                         level=c("N0", "N50")),
                              y = Phylogenetic_Diversity,
                              fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "b"),
           size = 4,
           x = c(1,2),
           y = 6.5) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.ita.rt.nitro.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#Philippine paddy soil
#Rhizosphere compartment
phylo.php.rz <- subset_samples(phylo, soil=="Philippine paddy soil" &
                                 compartment %in% c("Bulk soil", "Rhizosphere"))

phylo.php.rz.asvtab <- as.data.frame(phylo.php.rz@otu_table)

phylo.php.rz.tree <- phylo.php.rz@phy_tree

# phylo.php.rz.meta from previous code chunks

# We first need to check if the tree is Rhizosphereed or not 

phylo.php.rz@phy_tree
# it is an unRhizosphereed tree

phylo.php.rz.pd <- pd(t(phylo.php.rz.asvtab), phylo.php.rz.tree, include.root=F)
datatable(phylo.php.rz.pd)

#isolate the metadata file from the phyloseq object
phylo.php.rz.meta <- meta(phylo.php.rz)

#custom sort levels
phylo.php.rz.meta$genotype <- factor(phylo.php.rz.meta$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'Kasalath', 'IR64'))
phylo.php.rz.meta$nitro <- factor(phylo.php.rz.meta$nitro, levels=c('N0', 'N50'))

#add PD to the metadata file
phylo.php.rz.meta$Phylogenetic_Diversity <- phylo.php.rz.pd$PD

#Hypothesis testing
set.seed(999)

#ANOVA procedure on original data
phylo.php.rz.aov.org <- aov(
  (Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.php.rz.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

#rank transformed
phylo.php.rz.aov <- aov(
  rank(Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.php.rz.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

Anova(phylo.php.rz.aov, type = 'II')
Anova(phylo.php.rz.aov.org, type = 'II')

#means
#model.tables(phylo.php.rz.aov, type = "means", se = TRUE)
model.tables(phylo.php.rz.aov.org, type = "means", se = TRUE)

#Residuals
#phylo.php.rz.aov.res = phylo.php.rz.aov$resid
phylo.php.rz.aov.res = phylo.php.rz.aov.org$resid

#QQplot
qqnorm(
  phylo.php.rz.aov.res, pch = 20, main = "Rank-transformed data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(phylo.php.rz.aov.res)

#test for normality
shapiro.test(phylo.php.rz.aov.res)

#homogeneity of variance
#plot(phylo.php.rz.aov, 1, main = "Rank-transformed data")
plot(phylo.php.rz.aov.org, 1, main = "Original data")

#test for homogeneity of variance
#leveneTest(phylo.php.rz.aov)
leveneTest(phylo.php.rz.aov.org)

#Compute estimated marginal means for factor combinations
#emmeans(phylo.php.rz.aov, pairwise ~ genotype | nitro)
emmeans(phylo.php.rz.aov.org, pairwise ~ genotype | nitro)

#phylo.php.rz.aov.emmeans <- emmeans(phylo.php.rz.aov,  ~ genotype | nitro)
phylo.php.rz.aov.emmeans <- emmeans(phylo.php.rz.aov.org,  ~ genotype)
print(phylo.php.rz.aov.emmeans)

phylo.php.rz.aov.emmeans %>% 
  pairs() %>% 
  test(joint = TRUE)

pairs(phylo.php.rz.aov.emmeans)

#multiple comparisons
#summary(glht(phylo.php.rz.aov, linfct = mcp(genotype = "Tukey")))
summary(glht(phylo.php.rz.aov.org, linfct = mcp(genotype = "Tukey")))

#visualization
phylo.php.rz.lm <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = phylo.php.rz.meta)

#genotype
phylo.php.rz.aov.gen <- cld(emmeans(phylo.php.rz.lm, c("genotype")), alpha = 0.05, Letters = "letters")

phylo.php.rz.gen.mean<- phylo.php.rz.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.php.rz.aov.gen) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

phylo.php.rz.gen.jitter <- phylo.php.rz.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.php.rz.meta) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

ggplot(phylo.php.rz.meta, aes(x = factor(genotype,
                                         level=c("Bulk soil",
                                                 "Rufi",
                                                 "Nipponbare",
                                                 "Kasalath",
                                                 "IR64")),
                              y = Phylogenetic_Diversity,
                              fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "ab", "b", "ab", "ab"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 12) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.php.rz.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#nitrogen
ggplot(phylo.php.rz.meta, aes(x = factor(nitro,
                                         level=c("N0", "N50")),
                              y = Phylogenetic_Diversity,
                              fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 10) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.php.rz.nitro.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#root
#Philippine paddy soil
#Root compartment
phylo.php.rt <- subset_samples(phylo, soil=="Philippine paddy soil" &
                                 compartment %in% c("Bulk soil", "Root"))

phylo.php.rt.asvtab <- as.data.frame(phylo.php.rt@otu_table)

phylo.php.rt.tree <- phylo.php.rt@phy_tree

# phylo.php.rt.meta from previous code chunks

# We first need to check if the tree is Rooted or not 

phylo.php.rt@phy_tree
# it is an unRooted tree

phylo.php.rt.pd <- pd(t(phylo.php.rt.asvtab), phylo.php.rt.tree, include.root=F)
datatable(phylo.php.rt.pd)

#isolate the metadata file from the phyloseq object
phylo.php.rt.meta <- meta(phylo.php.rt)

#custom sort levels
phylo.php.rt.meta$genotype <- factor(phylo.php.rt.meta$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'Kasalath', 'IR64'))
phylo.php.rt.meta$nitro <- factor(phylo.php.rt.meta$nitro, levels=c('N0', 'N50'))

#add PD to the metadata file
phylo.php.rt.meta$Phylogenetic_Diversity <- phylo.php.rt.pd$PD

#Hypothesis testing
set.seed(999)

#ANOVA procedure on original data
phylo.php.rt.aov.org <- aov(
  (Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.php.rt.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

#rank transformed
phylo.php.rt.aov <- aov(
  rank(Phylogenetic_Diversity) ~ genotype*nitro, data = phylo.php.rt.meta,
  contrasts = list(
    genotype = 'contr.sum',
    nitro = 'contr.sum'
  )
)

Anova(phylo.php.rt.aov, type = 'II')
Anova(phylo.php.rt.aov.org, type = 'II')

#means
model.tables(phylo.php.rt.aov, type = "means", se = TRUE)
#model.tables(phylo.php.rt.aov.org, type = "means", se = TRUE)

#Residuals
phylo.php.rt.aov.res = phylo.php.rt.aov$resid
#phylo.php.rt.aov.res = phylo.php.rt.aov.org$resid

#QQplot
qqnorm(
  phylo.php.rt.aov.res, pch = 20, main = "Rank-transformed data",
  cex.lab = 1, cex.axis = 0.7, cex.main = 1
)
qqline(phylo.php.rt.aov.res)

#test for normality
shapiro.test(phylo.php.rt.aov.res)

#homogeneity of variance
plot(phylo.php.rt.aov, 1, main = "Rank-transformed data")
#plot(phylo.php.rt.aov.org, 1, main = "Original data")

#test for homogeneity of variance
leveneTest(phylo.php.rt.aov)
#leveneTest(phylo.php.rt.aov.org)

#Compute estimated marginal means for factor combinations
emmeans(phylo.php.rt.aov, pairwise ~ genotype | nitro)
#emmeans(phylo.php.rt.aov.org, pairwise ~ genotype | nitro)

phylo.php.rt.aov.emmeans <- emmeans(phylo.php.rt.aov,  ~ genotype | nitro)
#phylo.php.rt.aov.emmeans <- emmeans(phylo.php.rt.aov.org,  ~ genotype | nitro)
print(phylo.php.rt.aov.emmeans)

phylo.php.rt.aov.emmeans %>% 
  pairs() %>% 
  test(joint = TRUE)

pairs(phylo.php.rt.aov.emmeans)

#multiple comparisons
summary(glht(phylo.php.rt.aov, linfct = mcp(genotype = "Tukey")))
#summary(glht(phylo.php.rt.aov.org, linfct = mcp(genotype = "Tukey")))

#visualization
phylo.php.rt.lm <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = phylo.php.rt.meta)

#genotype
phylo.php.rt.aov.gen <- cld(emmeans(phylo.php.rt.lm, c("genotype")), alpha = 0.05, Letters = "letters")

phylo.php.rt.gen.mean<- phylo.php.rt.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.php.rt.aov.gen) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

phylo.php.rt.gen.jitter <- phylo.php.rt.meta %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(phylo.php.rt.meta) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

ggplot(phylo.php.rt.meta, aes(x = factor(genotype,
                                         level=c("Bulk soil",
                                                 "Rufi",
                                                 "Nipponbare",
                                                 "Kasalath",
                                                 "IR64")),
                              y = Phylogenetic_Diversity,
                              fill = genotype, group = genotype)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a", "a", "a", "a"),
           size = 4,
           x = c(1,2,3,4,5),
           y = 9.5) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#503D2E", "#048789", "#E2A72E", "#D44D27", "#EFEBC8")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.php.rt.genotype.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)

#nitrogen
ggplot(phylo.php.rt.meta, aes(x = factor(nitro,
                                         level=c("N0", "N50")),
                              y = Phylogenetic_Diversity,
                              fill = nitro, group = nitro)) +
  geom_boxplot() +
  annotate(geom="text", label = c("a", "a"),
           size = 4,
           x = c(1,2),
           y = 10) +
  scale_y_continuous(limits = c(3, 12),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#57A4B1", "#F95355")) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)") +
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
       filename = "alpha diversity.php.rt.nitro.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 500, height = 1480, unit = "px",
       dpi = 300)
