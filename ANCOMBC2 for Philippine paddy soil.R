#libraries needed
library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(DT)
library(tidyverse)
library(ggpubr)
library(mia)

#creating a phyloseq object
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
tax_table(phylo)[, colnames(tax_table(phylo))] <- gsub(tax_table(phylo)[, colnames(tax_table(phylo))], pattern = ".*__", replacement = "")

#renaming unclassified taxas with the last highest classified rank
ps1 <- phylo

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
datatable(tax_table(ps1))

#ANCOMBC2 for Nitrogen
#subset
ancombc2.Philippine.bulk <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                              compartment %in% c("Bulk soil"))

ancombc2.Philippine.bulk.tse <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.bulk)

#run ancombc2
set.seed(999)
ancombc2.Philippine.bulk.tse.output <- ancombc2(data = ancombc2.Philippine.bulk.tse,
                                            assay_name = "counts", tax_level = "Genus",
                  fix_formula = "nitro", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))

#structural zeros
tab_zero.Philippine.bulk = ancombc2.Philippine.bulk.tse.output$zero_ind
tab_zero.Philippine.bulk %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.bulk = ancombc2.Philippine.bulk.tse.output$pseudo_sens_tab
tab_sens.Philippine.bulk %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.bulk), digits = 2)

res_prim.Philippine.bulk <- ancombc2.Philippine.bulk.tse.output$res

df_Philippine.bulk.nitro = res_prim.Philippine.bulk %>%
  rownames_to_column("tax_id")

df_fig_Philippine.bulk.nitro = df_Philippine.bulk.nitro %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


df_fig_Philippine.bulk.nitro$tax_id = factor(df_fig_Philippine.bulk.nitro$tax_id,
                                          levels = df_fig_Philippine.bulk.nitro$tax_id)
df_fig_Philippine.bulk.nitro$direct = factor(df_fig_Philippine.bulk.nitro$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

df_fig_Philippine.bulk.nitro %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
                                   hjust = 1#, face = "italic"
                                   ),
          axis.text.y = element_text(#angle = 60,
            hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Bulk.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

#coord_flip
ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Bulk.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)


#subset: Rhizosphere
ancombc2.Philippine.Rhizosphere <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                                 compartment %in% c("Rhizosphere"))

ancombc2.Philippine.Rhizosphere.tse <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.Rhizosphere)

#run ancombc2
set.seed(999)
ancombc2.Philippine.Rhizosphere.tse.output <- ancombc2(data = ancombc2.Philippine.Rhizosphere.tse,
                                                    assay_name = "counts", tax_level = "Genus",
                                                    fix_formula = "nitro", rand_formula = NULL,
                                                    p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                                    prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                                    group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                                                    alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                    global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                                    iter_control = list(tol = 1e-2, max_iter = 20, 
                                                                        verbose = TRUE),
                                                    em_control = list(tol = 1e-5, max_iter = 100),
                                                    lme_control = lme4::lmerControl(),
                                                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                                    trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                                                nrow = 2, 
                                                                                                byrow = TRUE),
                                                                                         matrix(c(-1, 0, 1, -1),
                                                                                                nrow = 2, 
                                                                                                byrow = TRUE)),
                                                                         node = list(2, 2),
                                                                         solver = "ECOS",
                                                                         B = 100))

#structural zeros
tab_zero.Philippine.Rhizosphere = ancombc2.Philippine.Rhizosphere.tse.output$zero_ind
tab_zero.Philippine.Rhizosphere %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.Rhizosphere = ancombc2.Philippine.Rhizosphere.tse.output$pseudo_sens_tab
tab_sens.Philippine.Rhizosphere %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.Rhizosphere), digits = 2)

res_prim.Philippine.Rhizosphere <- ancombc2.Philippine.Rhizosphere.tse.output$res

df_Philippine.Rhizosphere.nitro = res_prim.Philippine.Rhizosphere %>%
  rownames_to_column("tax_id")

df_fig_Philippine.Rhizosphere.nitro = df_Philippine.Rhizosphere.nitro %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


df_fig_Philippine.Rhizosphere.nitro$tax_id = factor(df_fig_Philippine.Rhizosphere.nitro$tax_id,
                                                 levels = df_fig_Philippine.Rhizosphere.nitro$tax_id)
df_fig_Philippine.Rhizosphere.nitro$direct = factor(df_fig_Philippine.Rhizosphere.nitro$direct, 
                                                 levels = c("Positive LFC", "Negative LFC"))

df_fig_Philippine.Rhizosphere.nitro %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Rhizosphere.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

#coord_flip
ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Rhizosphere.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)

#subset: Root
ancombc2.Philippine.Root <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                          compartment %in% c("Root"))

ancombc2.Philippine.Root.tse <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.Root)

#run ancombc2
set.seed(999)
ancombc2.Philippine.Root.tse.output <- ancombc2(data = ancombc2.Philippine.Root.tse,
                                             assay_name = "counts", tax_level = "Genus",
                                             fix_formula = "nitro", rand_formula = NULL,
                                             p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                             prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                             group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                                             global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                             iter_control = list(tol = 1e-2, max_iter = 20, 
                                                                 verbose = TRUE),
                                             em_control = list(tol = 1e-5, max_iter = 100),
                                             lme_control = lme4::lmerControl(),
                                             mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                             trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                                         nrow = 2, 
                                                                                         byrow = TRUE),
                                                                                  matrix(c(-1, 0, 1, -1),
                                                                                         nrow = 2, 
                                                                                         byrow = TRUE)),
                                                                  node = list(2, 2),
                                                                  solver = "ECOS",
                                                                  B = 100))

#structural zeros
tab_zero.Philippine.Root = ancombc2.Philippine.Root.tse.output$zero_ind
tab_zero.Philippine.Root %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.Root = ancombc2.Philippine.Root.tse.output$pseudo_sens_tab
tab_sens.Philippine.Root %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.Root), digits = 2)

res_prim.Philippine.Root <- ancombc2.Philippine.Root.tse.output$res

df_Philippine.Root.nitro = res_prim.Philippine.Root %>%
  rownames_to_column("tax_id")

df_fig_Philippine.Root.nitro = df_Philippine.Root.nitro %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


df_fig_Philippine.Root.nitro$tax_id = factor(df_fig_Philippine.Root.nitro$tax_id,
                                          levels = df_fig_Philippine.Root.nitro$tax_id)
df_fig_Philippine.Root.nitro$direct = factor(df_fig_Philippine.Root.nitro$direct, 
                                          levels = c("Positive LFC", "Negative LFC"))

df_fig_Philippine.Root.nitro %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Root.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

#coord_flip
ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.Nitrogen.Root.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)

#genotype
#subset: Rufi
ancombc2.Philippine.B.Rufi <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                            genotype %in% c("Bulk soil", "Rufi"))

tse.Philippine.B.Rufi <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.B.Rufi)

#run ancombc2
set.seed(999)
output.Philippine.B.Rufi <- ancombc2(data = tse.Philippine.B.Rufi,
                                  assay_name = "counts", tax_level = "Genus",
                                  fix_formula = "genotype", rand_formula = NULL,
                                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                  group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                              nrow = 2, 
                                                                              byrow = TRUE),
                                                                       matrix(c(-1, 0, 1, -1),
                                                                              nrow = 2, 
                                                                              byrow = TRUE)),
                                                       node = list(2, 2),
                                                       solver = "ECOS",
                                                       B = 100))

#structural zeros
tab_zero.Philippine.B.Rufi = output.Philippine.B.Rufi$zero_ind
tab_zero.Philippine.B.Rufi %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.B.Rufi = output.Philippine.B.Rufi$pseudo_sens_tab
tab_sens.Philippine.B.Rufi %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.B.Rufi), digits = 2)

res_prim.Philippine.B.Rufi <- output.Philippine.B.Rufi$res

df.Philippine.B.Rufi = res_prim.Philippine.B.Rufi %>%
  rownames_to_column("tax_id")

df_fig.Philippine.B.Rufi = df.Philippine.B.Rufi %>%
  arrange(desc(lfc_genotypeRufi)) %>%
  mutate(direct = ifelse(lfc_genotypeRufi > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_genotypeRufi < 0.05, "Significant", "Not significant"))


df_fig.Philippine.B.Rufi$tax_id = factor(df_fig.Philippine.B.Rufi$tax_id,
                                      levels = df_fig.Philippine.B.Rufi$tax_id)
df_fig.Philippine.B.Rufi$direct = factor(df_fig.Philippine.B.Rufi$direct, 
                                      levels = c("Positive LFC", "Negative LFC"))

df_fig.Philippine.B.Rufi %>%
  ggplot(aes(x = tax_id, y = lfc_genotypeRufi, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_genotypeRufi - se_genotypeRufi, ymax = lfc_genotypeRufi + se_genotypeRufi), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Rufi.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

#coord_flip
ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Rufi.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)

#subset: Nipponbare
ancombc2.Philippine.B.Nipponbare <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                                  genotype %in% c("Bulk soil", "Nipponbare"))

tse.Philippine.B.Nipponbare <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.B.Nipponbare)

#run ancombc2
set.seed(999)
output.Philippine.B.Nipponbare <- ancombc2(data = tse.Philippine.B.Nipponbare,
                                        assay_name = "counts", tax_level = "Genus",
                                        fix_formula = "genotype", rand_formula = NULL,
                                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                        group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                                        alpha = 0.05, n_cl = 2, verbose = TRUE,
                                        global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                        iter_control = list(tol = 1e-2, max_iter = 20, 
                                                            verbose = TRUE),
                                        em_control = list(tol = 1e-5, max_iter = 100),
                                        lme_control = lme4::lmerControl(),
                                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                                    nrow = 2, 
                                                                                    byrow = TRUE),
                                                                             matrix(c(-1, 0, 1, -1),
                                                                                    nrow = 2, 
                                                                                    byrow = TRUE)),
                                                             node = list(2, 2),
                                                             solver = "ECOS",
                                                             B = 100))

#structural zeros
tab_zero.Philippine.B.Nipponbare = output.Philippine.B.Nipponbare$zero_ind
tab_zero.Philippine.B.Nipponbare %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.B.Nipponbare = output.Philippine.B.Nipponbare$pseudo_sens_tab
tab_sens.Philippine.B.Nipponbare %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.B.Nipponbare), digits = 2)

res_prim.Philippine.B.Nipponbare <- output.Philippine.B.Nipponbare$res

df.Philippine.B.Nipponbare = res_prim.Philippine.B.Nipponbare %>%
  rownames_to_column("tax_id")

df_fig.Philippine.B.Nipponbare = df.Philippine.B.Nipponbare %>%
  arrange(desc(lfc_genotypeNipponbare)) %>%
  mutate(direct = ifelse(lfc_genotypeNipponbare > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_genotypeNipponbare < 0.05, "Significant", "Not significant"))


df_fig.Philippine.B.Nipponbare$tax_id = factor(df_fig.Philippine.B.Nipponbare$tax_id,
                                            levels = df_fig.Philippine.B.Nipponbare$tax_id)
df_fig.Philippine.B.Nipponbare$direct = factor(df_fig.Philippine.B.Nipponbare$direct, 
                                            levels = c("Positive LFC", "Negative LFC"))

df_fig.Philippine.B.Nipponbare %>%
  ggplot(aes(x = tax_id, y = lfc_genotypeNipponbare, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_genotypeNipponbare - se_genotypeNipponbare, ymax = lfc_genotypeNipponbare + se_genotypeNipponbare), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Nipponbare.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Nipponbare.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)

#subset: Kasalath
ancombc2.Philippine.B.Kasalath <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                                genotype %in% c("Bulk soil", "Kasalath"))

tse.Philippine.B.Kasalath <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.B.Kasalath)

#run ancombc2
set.seed(999)
output.Philippine.B.Kasalath <- ancombc2(data = tse.Philippine.B.Kasalath,
                                      assay_name = "counts", tax_level = "Genus",
                                      fix_formula = "genotype", rand_formula = NULL,
                                      p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                      group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                                      alpha = 0.05, n_cl = 2, verbose = TRUE,
                                      global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                      iter_control = list(tol = 1e-2, max_iter = 20, 
                                                          verbose = TRUE),
                                      em_control = list(tol = 1e-5, max_iter = 100),
                                      lme_control = lme4::lmerControl(),
                                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                      trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                                  nrow = 2, 
                                                                                  byrow = TRUE),
                                                                           matrix(c(-1, 0, 1, -1),
                                                                                  nrow = 2, 
                                                                                  byrow = TRUE)),
                                                           node = list(2, 2),
                                                           solver = "ECOS",
                                                           B = 100))

#structural zeros
tab_zero.Philippine.B.Kasalath = output.Philippine.B.Kasalath$zero_ind
tab_zero.Philippine.B.Kasalath %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.B.Kasalath = output.Philippine.B.Kasalath$pseudo_sens_tab
tab_sens.Philippine.B.Kasalath %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.B.Kasalath), digits = 2)

res_prim.Philippine.B.Kasalath <- output.Philippine.B.Kasalath$res

df.Philippine.B.Kasalath = res_prim.Philippine.B.Kasalath %>%
  rownames_to_column("tax_id")

df_fig.Philippine.B.Kasalath = df.Philippine.B.Kasalath %>%
  arrange(desc(lfc_genotypeKasalath)) %>%
  mutate(direct = ifelse(lfc_genotypeKasalath > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_genotypeKasalath < 0.05, "Significant", "Not significant"))


df_fig.Philippine.B.Kasalath$tax_id = factor(df_fig.Philippine.B.Kasalath$tax_id,
                                          levels = df_fig.Philippine.B.Kasalath$tax_id)
df_fig.Philippine.B.Kasalath$direct = factor(df_fig.Philippine.B.Kasalath$direct, 
                                          levels = c("Positive LFC", "Negative LFC"))

df_fig.Philippine.B.Kasalath %>%
  ggplot(aes(x = tax_id, y = lfc_genotypeKasalath, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_genotypeKasalath - se_genotypeKasalath, ymax = lfc_genotypeKasalath + se_genotypeKasalath), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Kasalath.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_Kasalath.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)

#subset: IR64
ancombc2.Philippine.B.IR64 <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                            genotype %in% c("Bulk soil", "IR64"))

tse.Philippine.B.IR64 <- makeTreeSummarizedExperimentFromPhyloseq(ancombc2.Philippine.B.IR64)

#run ancombc2
set.seed(999)
output.Philippine.B.IR64 <- ancombc2(data = tse.Philippine.B.IR64,
                                  assay_name = "counts", tax_level = "Genus",
                                  fix_formula = "genotype", rand_formula = NULL,
                                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                  group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                              nrow = 2, 
                                                                              byrow = TRUE),
                                                                       matrix(c(-1, 0, 1, -1),
                                                                              nrow = 2, 
                                                                              byrow = TRUE)),
                                                       node = list(2, 2),
                                                       solver = "ECOS",
                                                       B = 100))

#structural zeros
tab_zero.Philippine.B.IR64 = output.Philippine.B.IR64$zero_ind
tab_zero.Philippine.B.IR64 %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
tab_sens.Philippine.B.IR64 = output.Philippine.B.IR64$pseudo_sens_tab
tab_sens.Philippine.B.IR64 %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(tab_sens.Philippine.B.IR64), digits = 2)

res_prim.Philippine.B.IR64 <- output.Philippine.B.IR64$res

df.Philippine.B.IR64 = res_prim.Philippine.B.IR64 %>%
  rownames_to_column("tax_id")

df_fig.Philippine.B.IR64 = df.Philippine.B.IR64 %>%
  arrange(desc(lfc_genotypeIR64)) %>%
  mutate(direct = ifelse(lfc_genotypeIR64 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_genotypeIR64 < 0.05, "Significant", "Not significant"))


df_fig.Philippine.B.IR64$tax_id = factor(df_fig.Philippine.B.IR64$tax_id,
                                      levels = df_fig.Philippine.B.IR64$tax_id)
df_fig.Philippine.B.IR64$direct = factor(df_fig.Philippine.B.IR64$direct, 
                                      levels = c("Positive LFC", "Negative LFC"))

df_fig.Philippine.B.IR64 %>%
  ggplot(aes(x = tax_id, y = lfc_genotypeIR64, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_genotypeIR64 - se_genotypeIR64, ymax = lfc_genotypeIR64 + se_genotypeIR64), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Bulk soil: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(#angle = 60,
          hjust = 1#, face = "italic"
        ),
        axis.text.y = element_text(#angle = 60,
          hjust = 1, face = "italic"))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_IR64.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2500, height = 1800, unit = "px",
       dpi = 300)

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.Philippine.genotype.Bulk_IR64.pdf",
       device = "pdf",
       path = "Final_figures/",
       width = 2000, height = 1200, unit = "px",
       dpi = 300)
