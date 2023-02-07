#libraries needed
library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(DT)
library(tidyverse)
library(ggpubr)
library(mia)
library(gghighlight)

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
datatable(tax_table(ps1))


#Italian paddy soil
#rhizosphere compartment
#bulk vs rufi

ita.rz <- subset_samples(ps1, soil == "Italian paddy soil" &
                           compartment %in% c("Bulk soil", "Rhizosphere"))

ita.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(ita.rz)

#run ancombc2
set.seed(999)
ita.rz.output <- ancombc2(data = ita.rz.tse,
                          assay_name = "counts", tax_level = "Genus",
                          fix_formula = "genotype", rand_formula = NULL,
                          p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = FALSE, dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
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
ita.rz.zero = ita.rz.output$zero_ind
ita.rz.zero %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
ita.rz.scores = ita.rz.output$pseudo_sens_tab
ita.rz.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(ita.rz.scores), digits = 2)

#primary results
ita.rz.prime_res <- ita.rz.output$res

#convert prime_res to dataframe
ita.rz.df <- ita.rz.prime_res %>%
  rownames_to_column("tax_id")

#dunnet type comparisons
ita.rz.dunn <- ita.rz.output$res_dunn

ita.rz.fig <- ita.rz.dunn %>%
  rownames_to_column("taxon") %>%
  transmute(taxon,
            `lfc_IR64` = round(lfc_genotypeIR64, 2),
            `lfc_Kasalath` = round(lfc_genotypeKasalath, 2),
            `lfc_Nipponbare` = round(lfc_genotypeNipponbare, 2),
            `lfc_Rufi` = round(lfc_genotypeRufi, 2),
            `sig_IR64` = round(q_genotypeIR64, 2),
            `sig_Kasalath` = round(q_genotypeKasalath, 2),
            `sig_Nipponbare` = round(q_genotypeNipponbare, 2),
            `sig_Rufi` = round(q_genotypeRufi, 2)) %>%
  pivot_longer(cols = -taxon,
               names_to = c(".value", "group"),
               names_sep = "_") %>%
  mutate(sig2 = ifelse(sig < 0.05, 1, 0)) %>%
  arrange(taxon)

lo = -3
up = 3
mid = 0

ita.rz.fig %>%
  ggplot(aes(x = factor(group, level = c("Rufi", "Nipponbare", "Kasalath", "IR64")),
             y = taxon, fill = lfc)) + 
  geom_tile(color = "black", alpha = 0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(group, taxon, label = lfc), color = "black", size = 3) +
  labs(x = NULL, y = NULL, subtitle = "Rhizosphere") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic", color = "black", size = 9),
        axis.text = element_text(size = 10),
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.ITA.rz.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1418, height = 1500, unit = "px",
       dpi = 300)

#Italian paddy soil
#Root compartment

ita.rt <- subset_samples(ps1, soil == "Italian paddy soil" &
                           compartment %in% c("Bulk soil", "Root"))

ita.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(ita.rt)

#run ancombc2
set.seed(999)
ita.rt.output <- ancombc2(data = ita.rt.tse,
                          assay_name = "counts", tax_level = "Genus",
                          fix_formula = "genotype", rand_formula = NULL,
                          p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = FALSE, dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
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
ita.rt.zero = ita.rt.output$zero_ind
ita.rt.zero %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
ita.rt.scores = ita.rt.output$pseudo_sens_tab
ita.rt.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(ita.rt.scores), digits = 2)

#primary results
ita.rt.prime_res <- ita.rt.output$res

#convert prime_res to dataframe
ita.rt.df <- ita.rt.prime_res %>%
  rownames_to_column("tax_id")

#dunnet type comparisons
ita.rt.dunn <- ita.rt.output$res_dunn

ita.rt.fig <- ita.rt.dunn %>%
  rownames_to_column("taxon") %>%
  transmute(taxon,
            `lfc_IR64` = round(lfc_genotypeIR64, 2),
            `lfc_Kasalath` = round(lfc_genotypeKasalath, 2),
            `lfc_Nipponbare` = round(lfc_genotypeNipponbare, 2),
            `lfc_Rufi` = round(lfc_genotypeRufi, 2),
            `sig_IR64` = round(q_genotypeIR64, 2),
            `sig_Kasalath` = round(q_genotypeKasalath, 2),
            `sig_Nipponbare` = round(q_genotypeNipponbare, 2),
            `sig_Rufi` = round(q_genotypeRufi, 2)) %>%
  pivot_longer(cols = -taxon,
               names_to = c(".value", "group"),
               names_sep = "_") %>%
  arrange(taxon)

lo = -3
up = 3
mid = 0

ita.rt.fig %>%
  mutate(sig2 = ifelse(sig < 0.05, 1, 0)) %>%
  ggplot(aes(x = factor(group, level = c("Rufi", "Nipponbare", "Kasalath", "IR64")),
             y = taxon, fill = lfc)) + 
  geom_tile(color = "black", alpha = 0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(group, taxon, label = lfc), color = "black", size = 3) +
  labs(x = NULL, y = NULL, subtitle = "Root") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic", color = "black", size = 9),
        axis.text = element_text(size = 10),
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.ITA.rt.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1480, height = 1500, unit = "px",
       dpi = 300)

#Philippine paddy soil
#Rhizosphere compartment

php.rz <- subset_samples(ps1, soil == "Philippine paddy soil" &
                           compartment %in% c("Bulk soil", "Rhizosphere"))

php.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(php.rz)

#run ancombc2
set.seed(999)
php.rz.output <- ancombc2(data = php.rz.tse,
                          assay_name = "counts", tax_level = "Genus",
                          fix_formula = "genotype", rand_formula = NULL,
                          p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = FALSE, dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
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
php.rz.zero = php.rz.output$zero_ind
php.rz.zero %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
php.rz.scores = php.rz.output$pseudo_sens_tab
php.rz.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(php.rz.scores), digits = 2)

#primary results
php.rz.prime_res <- php.rz.output$res

#convert prime_res to dataframe
php.rz.df <- php.rz.prime_res %>%
  rownames_to_column("tax_id")

#dunnet type comparisons
php.rz.dunn <- php.rz.output$res_dunn

php.rz.fig <- php.rz.dunn %>%
  rownames_to_column("taxon") %>%
  transmute(taxon,
            `lfc_IR64` = round(lfc_genotypeIR64, 2),
            `lfc_Kasalath` = round(lfc_genotypeKasalath, 2),
            `lfc_Nipponbare` = round(lfc_genotypeNipponbare, 2),
            `lfc_Rufi` = round(lfc_genotypeRufi, 2),
            `sig_IR64` = round(q_genotypeIR64, 2),
            `sig_Kasalath` = round(q_genotypeKasalath, 2),
            `sig_Nipponbare` = round(q_genotypeNipponbare, 2),
            `sig_Rufi` = round(q_genotypeRufi, 2)) %>%
  pivot_longer(cols = -taxon,
               names_to = c(".value", "group"),
               names_sep = "_") %>%
  arrange(taxon)

lo = -3
up = 3
mid = 0

php.rz.fig %>%
  mutate(sig2 = ifelse(sig < 0.05, 1, 0)) %>%
  ggplot(aes(x = factor(group, level = c("Rufi", "Nipponbare", "Kasalath", "IR64")),
             y = taxon, fill = lfc)) + 
  geom_tile(aes(alpha = sig2), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(group, taxon, label = lfc), color = "black", size = 3) +
  labs(x = NULL, y = NULL, subtitle = "Rhizosphere") +
  scale_alpha_continuous(guide="none") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic", color = "black", size = 9),
        axis.text = element_text(size = 10),
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.php.rz.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#Philippine paddy soil
#Root compartment

php.rt <- subset_samples(ps1, soil == "Philippine paddy soil" &
                           compartment %in% c("Bulk soil", "Root"))

php.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(php.rt)

#run ancombc2
set.seed(999)
php.rt.output <- ancombc2(data = php.rt.tse,
                          assay_name = "counts", tax_level = "Genus",
                          fix_formula = "genotype", rand_formula = NULL,
                          p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                          prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "genotype", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = FALSE, dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
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
php.rt.zero = php.rt.output$zero_ind
php.rt.zero %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
php.rt.scores = php.rt.output$pseudo_sens_tab
php.rt.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(php.rt.scores), digits = 2)

#primary results
php.rt.prime_res <- php.rt.output$res

#convert prime_res to dataframe
php.rt.df <- php.rt.prime_res %>%
  rownames_to_column("tax_id")

#dunnet type comparisons
php.rt.dunn <- php.rt.output$res_dunn

php.rt.fig <- php.rt.dunn %>%
  rownames_to_column("taxon") %>%
  transmute(taxon,
            `lfc_IR64` = round(lfc_genotypeIR64, 2),
            `lfc_Kasalath` = round(lfc_genotypeKasalath, 2),
            `lfc_Nipponbare` = round(lfc_genotypeNipponbare, 2),
            `lfc_Rufi` = round(lfc_genotypeRufi, 2),
            `sig_IR64` = round(q_genotypeIR64, 2),
            `sig_Kasalath` = round(q_genotypeKasalath, 2),
            `sig_Nipponbare` = round(q_genotypeNipponbare, 2),
            `sig_Rufi` = round(q_genotypeRufi, 2)) %>%
  pivot_longer(cols = -taxon,
               names_to = c(".value", "group"),
               names_sep = "_") %>%
  arrange(taxon)

lo = -3
up = 3
mid = 0

php.rt.fig %>%
  mutate(sig2 = ifelse(sig <= 0.05, 1, 0)) %>%
  ggplot(aes(x = factor(group, level = c("Rufi", "Nipponbare", "Kasalath", "IR64")),
             y = taxon, fill = lfc)) + 
  geom_tile(aes(alpha = sig2), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(group, taxon, label = lfc), color = "black", size = 3) +
  labs(x = NULL, y = NULL, subtitle = "Root") +
  scale_alpha_continuous(guide="none") +
  theme_bw() +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(color="black", size=10),
        axis.text.x = element_text(color="black", size = 10,
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic", color = "black", size = 9),
        axis.text = element_text(size = 10),
        legend.position = "right")

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.php.rt.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)


#######
#######
#######
#Italian paddy soil
#nitro
#bulk soil
ita.nit.blk <- subset_samples(ps1, soil == "Italian paddy soil" & 
                                compartment %in% c("Bulk soil"))

ita.nit.blk.tse <- makeTreeSummarizedExperimentFromPhyloseq(ita.nit.blk)

#run ancombc2
set.seed(999)
ita.nit.blk.output <- ancombc2(data = ita.nit.blk.tse,
                               assay_name = "counts", tax_level = "Genus",
                               fix_formula = "nitro", rand_formula = NULL,
                               p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                           nrow = 2, byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1),
                                                                           nrow = 2, byrow = TRUE)),
                                                    node = list(2, 2),
                                                    solver = "ECOS", B = 100))

#structural zeros
ita.nit.blk.zeros = ita.nit.blk.output$zero_ind
ita.nit.blk.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
ita.nit.blk.scores = ita.nit.blk.output$pseudo_sens_tab
ita.nit.blk.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(ita.nit.blk.scores), digits = 2)

ita.nit.blk.prime_res <- ita.nit.blk.output$res

ita.nit.blk.prime_res.df = ita.nit.blk.prime_res %>%
  rownames_to_column("tax_id")

ita.nit.blk.plot = ita.nit.blk.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


ita.nit.blk.plot$tax_id = factor(ita.nit.blk.plot$tax_id,
                                 levels = ita.nit.blk.plot$tax_id)
ita.nit.blk.plot$direct = factor(ita.nit.blk.plot$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

ita.nit.blk.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2, 2.6),
                     breaks = scales::pretty_breaks(n = 5)) +
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
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.ITA.nit.blk.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#Rhizosphere
ita.nit.rz <- subset_samples(ps1, soil == "Italian paddy soil" & 
                               compartment %in% c("Rhizosphere"))

ita.nit.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(ita.nit.rz)

#run ancombc2
set.seed(999)
ita.nit.rz.output <- ancombc2(data = ita.nit.rz.tse,
                              assay_name = "counts", tax_level = "Genus",
                              fix_formula = "nitro", rand_formula = NULL,
                              p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1),
                                                                          nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2),
                                                   solver = "ECOS", B = 100))

#structural zeros
ita.nit.rz.zeros = ita.nit.rz.output$zero_ind
ita.nit.rz.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
ita.nit.rz.scores = ita.nit.rz.output$pseudo_sens_tab
ita.nit.rz.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(ita.nit.rz.scores), digits = 2)

ita.nit.rz.prime_res <- ita.nit.rz.output$res

ita.nit.rz.prime_res.df = ita.nit.rz.prime_res %>%
  rownames_to_column("tax_id")

ita.nit.rz.plot = ita.nit.rz.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


ita.nit.rz.plot$tax_id = factor(ita.nit.rz.plot$tax_id,
                                levels = ita.nit.rz.plot$tax_id)
ita.nit.rz.plot$direct = factor(ita.nit.rz.plot$direct, 
                                levels = c("Positive LFC", "Negative LFC"))

ita.nit.rz.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Rhizosphere: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.ITA.nit.rz.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)


#Root
ita.nit.rt <- subset_samples(ps1, soil == "Italian paddy soil" & 
                               compartment %in% c("Root"))

ita.nit.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(ita.nit.rt)

#run ancombc2
set.seed(999)
ita.nit.rt.output <- ancombc2(data = ita.nit.rt.tse,
                              assay_name = "counts", tax_level = "Genus",
                              fix_formula = "nitro", rand_formula = NULL,
                              p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1),
                                                                          nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2),
                                                   solver = "ECOS", B = 100))

#structural zeros
ita.nit.rt.zeros = ita.nit.rt.output$zero_ind
ita.nit.rt.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
ita.nit.rt.scores = ita.nit.rt.output$pseudo_sens_tab
ita.nit.rt.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(ita.nit.rt.scores), digits = 2)

ita.nit.rt.prime_res <- ita.nit.rt.output$res

ita.nit.rt.prime_res.df = ita.nit.rt.prime_res %>%
  rownames_to_column("tax_id")

ita.nit.rt.plot = ita.nit.rt.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


ita.nit.rt.plot$tax_id = factor(ita.nit.rt.plot$tax_id,
                                levels = ita.nit.rt.plot$tax_id)
ita.nit.rt.plot$direct = factor(ita.nit.rt.plot$direct, 
                                levels = c("Positive LFC", "Negative LFC"))

ita.nit.rt.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Root: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.ITA.nit.rt.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#####
####
#Philippine paddy soil
#nitro
#bulk soil
php.nit.blk <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                                compartment %in% c("Bulk soil"))

php.nit.blk.tse <- makeTreeSummarizedExperimentFromPhyloseq(php.nit.blk)

#run ancombc2
set.seed(999)
php.nit.blk.output <- ancombc2(data = php.nit.blk.tse,
                               assay_name = "counts", tax_level = "Genus",
                               fix_formula = "nitro", rand_formula = NULL,
                               p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                               prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                           nrow = 2, byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1),
                                                                           nrow = 2, byrow = TRUE)),
                                                    node = list(2, 2),
                                                    solver = "ECOS", B = 100))

#structural zeros
php.nit.blk.zeros = php.nit.blk.output$zero_ind
php.nit.blk.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
php.nit.blk.scores = php.nit.blk.output$pseudo_sens_tab
php.nit.blk.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(php.nit.blk.scores), digits = 2)

php.nit.blk.prime_res <- php.nit.blk.output$res

php.nit.blk.prime_res.df = php.nit.blk.prime_res %>%
  rownames_to_column("tax_id")

php.nit.blk.plot = php.nit.blk.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


php.nit.blk.plot$tax_id = factor(php.nit.blk.plot$tax_id,
                                 levels = php.nit.blk.plot$tax_id)
php.nit.blk.plot$direct = factor(php.nit.blk.plot$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

php.nit.blk.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2.5, 2.5),
                     breaks = scales::pretty_breaks(n = 5)) +
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
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.PHP.nit.blk.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)

#Rhizosphere
php.nit.rz <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                               compartment %in% c("Rhizosphere"))

php.nit.rz.tse <- makeTreeSummarizedExperimentFromPhyloseq(php.nit.rz)

#run ancombc2
set.seed(999)
php.nit.rz.output <- ancombc2(data = php.nit.rz.tse,
                              assay_name = "counts", tax_level = "Genus",
                              fix_formula = "nitro", rand_formula = NULL,
                              p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1),
                                                                          nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2),
                                                   solver = "ECOS", B = 100))

#structural zeros
php.nit.rz.zeros = php.nit.rz.output$zero_ind
php.nit.rz.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
php.nit.rz.scores = php.nit.rz.output$pseudo_sens_tab
php.nit.rz.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(php.nit.rz.scores), digits = 2)

php.nit.rz.prime_res <- php.nit.rz.output$res

php.nit.rz.prime_res.df = php.nit.rz.prime_res %>%
  rownames_to_column("tax_id")

php.nit.rz.plot = php.nit.rz.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


php.nit.rz.plot$tax_id = factor(php.nit.rz.plot$tax_id,
                                levels = php.nit.rz.plot$tax_id)
php.nit.rz.plot$direct = factor(php.nit.rz.plot$direct, 
                                levels = c("Positive LFC", "Negative LFC"))

php.nit.rz.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Rhizosphere: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.PHP.nit.rz.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)


#Root
php.nit.rt <- subset_samples(ps1, soil == "Philippine paddy soil" & 
                               compartment %in% c("Root"))

php.nit.rt.tse <- makeTreeSummarizedExperimentFromPhyloseq(php.nit.rt)

#run ancombc2
set.seed(999)
php.nit.rt.output <- ancombc2(data = php.nit.rt.tse,
                              assay_name = "counts", tax_level = "Genus",
                              fix_formula = "nitro", rand_formula = NULL,
                              p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "nitro", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                              iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2, byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1),
                                                                          nrow = 2, byrow = TRUE)),
                                                   node = list(2, 2),
                                                   solver = "ECOS", B = 100))

#structural zeros
php.nit.rt.zeros = php.nit.rt.output$zero_ind
php.nit.rt.zeros %>%
  datatable(caption = "The detection of structural zeros")

#sensitivity scores
php.nit.rt.scores = php.nit.rt.output$pseudo_sens_tab
php.nit.rt.scores %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(php.nit.rt.scores), digits = 2)

php.nit.rt.prime_res <- php.nit.rt.output$res

php.nit.rt.prime_res.df = php.nit.rt.prime_res %>%
  rownames_to_column("tax_id")

php.nit.rt.plot = php.nit.rt.prime_res.df %>%
  arrange(desc(lfc_nitroN50)) %>%
  mutate(direct = ifelse(lfc_nitroN50 > 0, "Positive LFC", "Negative LFC")) %>%
  mutate(signif = ifelse(q_nitroN50 < 0.05, "Significant", "Not significant"))


php.nit.rt.plot$tax_id = factor(php.nit.rt.plot$tax_id,
                                levels = php.nit.rt.plot$tax_id)
php.nit.rt.plot$direct = factor(php.nit.rt.plot$direct, 
                                levels = c("Positive LFC", "Negative LFC"))

php.nit.rt.plot %>%
  ggplot(aes(x = tax_id, y = lfc_nitroN50, fill = signif)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_nitroN50 - se_nitroN50, ymax = lfc_nitroN50 + se_nitroN50), 
                width = 0.2, position = position_dodge(0.05), color = "black") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values=c("#73777B", "#FE7E6D"),
                    name = NULL) +
  labs(x = NULL, y = "Log fold change (LFC)", 
       subtitle = "Root: N0 vs N50") +
  coord_flip() +
  scale_color_discrete(name = NULL) +
  geom_hline(yintercept = 0) +
  theme_bw() + 
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 12),
        axis.text.x = element_text( color = "black"),
        legend.position = c(0.8,0.8))

ggsave(plot = last_plot(),
       filename = "ANCOMBC2.PHP.nit.rt.pdf",
       device = "pdf",
       path = "GITHUB_figures/",
       width = 1500, height = 1500, unit = "px",
       dpi = 300)