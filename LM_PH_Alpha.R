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

phylo2 = subset_samples(phylo, soil=="Philippine paddy soil" & 
                          compartment %in% c("Bulk soil", "Rhizosphere"))

set.seed(999)

ps0.rar.asvtab <- as.data.frame(phylo2@otu_table)

ps0.rar.tree <- phylo2@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 

phylo2@phy_tree

# it is an unrooted tree
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=F)
datatable(df.pd)

#isolate the metadata file from the phyloseq object
hmp.meta.orig <- meta(phylo2)

#add PD to the metadata file
hmp.meta.orig$Phylogenetic_Diversity <- df.pd$PD

rhizosphere.df <- hmp.meta.orig

#convert chr to factor and custom sort levels
rhizosphere.df$genotype <- factor(rhizosphere.df$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'IR64', 'Kasalath'))
rhizosphere.df$nitro <- factor(rhizosphere.df$nitro, levels=c('N0', 'N50'))

#assigning contrasts
#genotype
levels<-5
categories<-seq(levels, 2)
basematrix=matrix(-1, nrow=levels, ncol=levels)
diag(basematrix[1:levels, 2:levels])<-seq(levels-1, 1)
sub.basematrix<-basematrix[,2:levels]
sub.basematrix[upper.tri(sub.basematrix-1)]<-0
contrasts_5<-sub.basematrix %*% diag(1/categories)
contrasts(rhizosphere.df$genotype) <- contrasts_5

#nitrogen
nitrogen_con <- matrix(c(0.5, -0.5), ncol=1)
contrasts(rhizosphere.df$nitro) <- nitrogen_con


#Hypothesis testing
set.seed(999)

lm.rhizosphere <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = rhizosphere.df)
aov.rhizosphere <- anova_test(Phylogenetic_Diversity ~ genotype*nitro, data = rhizosphere.df)
aov.rhizosphere

#posthoc
pwc.rhizosphere <- rhizosphere.df %>%
  emmeans_test(
    Phylogenetic_Diversity ~ genotype, p.adjust.method = "bonferroni",
    model = lm.rhizosphere
  )

#assumption for normality
shapiro.test(rhizosphere.residuals)
#W = 0.85781, p-value = 5.854e-05

#assumption for homoscedasticity
leveneTest(residuals(rhizosphere.lmer) ~ genotype*nitrogen, data = rhizosphere.df)
#group: df=9, F=0.7827, P=0.6335

#post-hoc test
rhizosphere.lmer.posthoc <- cld(emmeans(lm.rhizosphere, c("genotype")), alpha = 0.05, Letters = "letters")

#visualize soil
mean.GN.rhizosphere <- rhizosphere.df %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(rhizosphere.lmer.posthoc) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

jitter.GN.rhizosphre <- rhizosphere.df %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(rhizosphere.df) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

position.rhizosphere <- pwc.rhizosphere %>%
  mutate(group1 = str_replace(group1, "Bulk soil", "Bulk soil\n(n=10)")) %>%
  mutate(group1 = str_replace(group1, "Rufi", "Rufi\n(n=9)")) %>%
  mutate(group1 = str_replace(group1, "Nipponbare", "Nipponbare\n(n=9)")) %>%
  mutate(group1 = str_replace(group1, "Kasalath", "Kasalath\n(n=8)")) %>%
  mutate(group1 = str_replace(group1, "IR64", "IR64\n(n=10)")) %>%
  mutate(group2 = str_replace(group2, "Bulk soil", "Bulk soil\n(n=10)")) %>%
  mutate(group2 = str_replace(group2, "Rufi", "Rufi\n(n=9)")) %>%
  mutate(group2 = str_replace(group2, "Nipponbare", "Nipponbare\n(n=9)")) %>%
  mutate(group2 = str_replace(group2, "Kasalath", "Kasalath\n(n=8)")) %>%
  mutate(group2 = str_replace(group2, "IR64", "IR64\n(n=10)"))

ggplot(mean.GN.rhizosphere, aes(x = factor(genotype,
                                                level=c("Bulk soil\n(n=10)",
                                                        "Rufi\n(n=9)",
                                                        "Nipponbare\n(n=9)",
                                                        "Kasalath\n(n=8)",
                                                        "IR64\n(n=10)")),
                                y = emmean, fill = nitro)) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                position=position_dodge(0.3), width = 0.1, size = 0.8) +
  geom_point(size = 5, shape = 22, stroke = 1, position=position_dodge(0.3)) +
  geom_jitter(data = jitter.GN.rhizosphre, aes(x = factor(genotype,
                                                             level=c("Bulk soil\n(n=10)",
                                                                     "Rufi\n(n=9)",
                                                                     "Nipponbare\n(n=9)",
                                                                     "Kasalath\n(n=8)",
                                                                     "IR64\n(n=10)")), 
                                           y = Phylogenetic_Diversity,
                                           fill = nitro),
              shape = 21, color = "black", alpha=0.8,
              position=position_jitter(0.2)) +
  scale_y_continuous(limits = c(3, 11.8),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#A89957", "#5766A8")) +
  scale_shape_manual(values=c("24", "21")) +
  annotate(geom="text", label = c("a", "ab", "bc", "ab", "ab"), size = 6,
           x = c(1,2,3,4,5),
           y = c(8, 9, 10, 9.5, 8.8)) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)",
       subtitle = get_test_label(aov.rhizosphere, detailed = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 12),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 12),
        legend.position = "top")

ggsave(plot = last_plot(),
       filename = "Alpha.ph.rhizosphere.pdf",
       device = "pdf",
       path = "final_figures/",
       width = 1500, height = 1800, unit = "px",
       dpi = 300)

#Root subset
phylo3 = subset_samples(phylo, soil=="Philippine paddy soil" & 
                          compartment %in% c("Bulk soil", "Root"))

set.seed(999)

ps0.rar.asvtab <- as.data.frame(phylo3@otu_table)

ps0.rar.tree <- phylo3@phy_tree

# hmp.meta from previous code chunks

# We first need to check if the tree is rooted or not 

phylo3@phy_tree

# it is an unrooted tree
df.pd <- pd(t(ps0.rar.asvtab), ps0.rar.tree,include.root=F)
datatable(df.pd)

#isolate the metadata file from the phyloseq object
hmp.meta.orig <- meta(phylo3)

#add PD to the metadata file
hmp.meta.orig$Phylogenetic_Diversity <- df.pd$PD

Root.df <- hmp.meta.orig

#convert chr to factor and custom sort levels
Root.df$genotype <- factor(Root.df$genotype, levels=c('Bulk soil', 'Rufi', 'Nipponbare', 'IR64', 'Kasalath'))
Root.df$nitro <- factor(Root.df$nitro, levels=c('N0', 'N50'))

#assigning contrasts
#genotype
levels<-5
categories<-seq(levels, 2)
basematrix=matrix(-1, nrow=levels, ncol=levels)
diag(basematrix[1:levels, 2:levels])<-seq(levels-1, 1)
sub.basematrix<-basematrix[,2:levels]
sub.basematrix[upper.tri(sub.basematrix-1)]<-0
contrasts_5<-sub.basematrix %*% diag(1/categories)
contrasts(Root.df$genotype) <- contrasts_5

#nitrogen
nitrogen_con <- matrix(c(0.5, -0.5), ncol=1)
contrasts(Root.df$nitro) <- nitrogen_con

#Hypothesis testing
set.seed(999)

lm.root <- lm(Phylogenetic_Diversity ~ genotype*nitro, data = Root.df)
aov.root <- anova_test(Phylogenetic_Diversity ~ genotype*nitro, data = Root.df)
aov.root

#posthoc
pwc.root <- Root.df %>%
  emmeans_test(
    Phylogenetic_Diversity ~ genotype, p.adjust.method = "bonferroni",
    model = lm.root
  )

#assumption for normality
shapiro.test(rhizosphere.residuals)
#W = 0.85781, p-value = 5.854e-05

#assumption for homoscedasticity
leveneTest(residuals(rhizosphere.lmer) ~ genotype*nitrogen, data = rhizosphere.df)
#group: df=9, F=0.7827, P=0.6335

#post-hoc test
root.lmer.posthoc <- cld(emmeans(lm.root, c("genotype")), alpha = 0.05, Letters = "letters")

#visualize soil
mean.GN.root <- Root.df %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(root.lmer.posthoc) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

jitter.GN.root <- Root.df %>% group_by(genotype) %>% summarize(num=n()) %>%
  left_join(Root.df) %>%
  mutate(genotype = paste0(genotype, "\n", "(n=", num, ")"))

position.root <- pwc.root %>%
  mutate(group1 = str_replace(group1, "Bulk soil", "Bulk soil\n(n=10)")) %>%
  mutate(group1 = str_replace(group1, "Rufi", "Rufi\n(n=9)")) %>%
  mutate(group1 = str_replace(group1, "Nipponbare", "Nipponbare\n(n=10)")) %>%
  mutate(group1 = str_replace(group1, "Kasalath", "Kasalath\n(n=11)")) %>%
  mutate(group1 = str_replace(group1, "IR64", "IR64\n(n=10)")) %>%
  mutate(group2 = str_replace(group2, "Bulk soil", "Bulk soil\n(n=10)")) %>%
  mutate(group2 = str_replace(group2, "Rufi", "Rufi\n(n=9)")) %>%
  mutate(group2 = str_replace(group2, "Nipponbare", "Nipponbare\n(n=10)")) %>%
  mutate(group2 = str_replace(group2, "Kasalath", "Kasalath\n(n=11)")) %>%
  mutate(group2 = str_replace(group2, "IR64", "IR64\n(n=10)"))

ggplot(mean.GN.root, aes(x = factor(genotype,
                                           level=c("Bulk soil\n(n=10)",
                                                   "Rufi\n(n=9)",
                                                   "Nipponbare\n(n=10)",
                                                   "Kasalath\n(n=11)",
                                                   "IR64\n(n=10)")),
                                y = emmean, fill = nitro)) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                position=position_dodge(0.3), width = 0.1, size = 0.8) +
  geom_point(size = 5, shape = 22, stroke = 1, position=position_dodge(0.3)) +
  geom_jitter(data = jitter.GN.root, aes(x = factor(genotype,
                                                          level=c("Bulk soil\n(n=10)",
                                                                  "Rufi\n(n=9)",
                                                                  "Nipponbare\n(n=10)",
                                                                  "Kasalath\n(n=11)",
                                                                  "IR64\n(n=10)")), 
                                               y = Phylogenetic_Diversity,
                                               fill = nitro),
              shape = 21, color = "black", alpha=0.8,
              position=position_jitter(0.2)) +
  scale_y_continuous(limits = c(3, 11.8),
                     breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values=c("#A89957", "#5766A8")) +
  scale_shape_manual(values=c("24", "21")) +
  annotate(geom="text", label = c("a", "ab", "ab", "ab", "bc"), size = 6,
           x = c(1,2,3,4,5),
           y = c(7.7, 9, 8.8, 8.5, 8.8)) +
  labs(x = "",
       y = "Phylogenetic Diversity (PD)",
       subtitle = get_test_label(aov.root, detailed = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size = 12),
        axis.text.y = element_text(color="black", size = 10),
        axis.text = element_text(size = 12),
        legend.position = "top")

ggsave(plot = last_plot(),
       filename = "Alpha.ph.root.pdf",
       device = "pdf",
       path = "final_figures/",
       width = 1500, height = 1800, unit = "px",
       dpi = 300)
