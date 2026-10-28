# AIM ---------------------------------------------------------------------
# Pseudobulk DE on nerve EC + FIBRO aggregated counts, with the donor (Sample_Name) as
# blocking factor:
#   gene ~ Sample_Name + CellType * Compartment
# Contrasts (Epineurial vs Endoneurial):
#   1. Compartment_EC: effect in EC
#   2. Compartment_FIBRO: effect in FIBRO
#   3. Compartment_avg: effect averaged over EC and FIBRO (= adjusted for CellType)
#   4. Interaction_FIBRO_vs_EC: does the compartment effect differ between FIBRO and EC?
# The genes DE in both EC and FIBRO (same direction) are compared with Compartment_avg.
#
# Notes
# - Perineurial exists only in FIBRO (not in EC), so it cannot be part of a joint CellType * Compartment model (its main effect would be collinear with the FIBRO:Perineurial interaction). The joint model uses Endoneurial + Epineurial only.
# - With treatment coding (reference = EC / Endoneurial) the "Compartment" coefficient alone is the effect in EC only; FIBRO is Compartment + CellType:Compartment; the average over EC and FIBRO is Compartment + 0.5 * CellType:Compartment.
# - A sample that has only one of the two compartments within a cell type adds no information to the compartment contrasts. See the pairing table printed below.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(ashr)
library(RNAseqQC)
library(ggrepel)
library(patchwork)
library(ComplexUpset)

# parameters --------------------------------------------------------------
data_dir <- "../../data"
out_dir <- "../../out"

MIN_CELLS <- 20             # drop pseudobulks with fewer cells than this (0 = keep all)
padj_thr <- 0.05
lfc_thr <- 1               # used for volcano/heatmap highlighting
max_hm_genes <- 200        # cap on genes drawn in heatmaps

# colours for CellType and Compartment (plots and heatmap annotations)
col_hm <- list(CellType = c(EC = "firebrick", FIBRO = "steelblue"),
               Compartment = c(Endoneurial = "gold", Epineurial = "purple"))

# read in the data --------------------------------------------------------
# load the count matrices (first column = gene symbol, other columns = PB_ID)
ec_df <- read_csv(file.path(data_dir, "pseudobulk_EC_counts.csv"), show_col_types = FALSE) %>%
  dplyr::rename(gene = 1)
fib_df <- read_csv(file.path(data_dir, "pseudobulk_FIBRO_counts.csv"), show_col_types = FALSE) %>%
  dplyr::rename(gene = 1)
stopifnot(!anyDuplicated(ec_df$gene), !anyDuplicated(fib_df$gene))

ec_mat <- as.matrix(ec_df[, -1])
rownames(ec_mat) <- ec_df$gene
fib_mat <- as.matrix(fib_df[, -1])
rownames(fib_mat) <- fib_df$gene

# load the meta
ec_meta <- read_csv(file.path(data_dir, "pseudobulk_EC_metadata.csv"), show_col_types = FALSE)
fib_meta <- read_csv(file.path(data_dir, "pseudobulk_FIBRO_metadata.csv"), show_col_types = FALSE)

# every count column must have a metadata row (and vice versa)
stopifnot(setequal(colnames(ec_mat), ec_meta$PB_ID),
          setequal(colnames(fib_mat), fib_meta$PB_ID))

# make sure the genes are the same and in the same order
genes <- intersect(rownames(ec_mat), rownames(fib_mat))
message("genes EC: ", nrow(ec_mat), " | FIBRO: ", nrow(fib_mat), " | shared: ", length(genes))

counts <- cbind(ec_mat[genes, ], fib_mat[genes, ])
stopifnot(all(counts == round(counts)))
storage.mode(counts) <- "integer"

colData <- bind_rows(ec_meta, fib_meta) %>%
  # joint model: only the compartments shared by EC and FIBRO
  filter(Compartment %in% c("Endoneurial", "Epineurial")) %>%
  mutate(CellType = factor(CellType, levels = c("EC", "FIBRO")),
         Compartment = factor(Compartment, levels = c("Endoneurial", "Epineurial")),
         Sample_Name = factor(Sample_Name)) %>%
  as.data.frame() %>%
  # add the covariate for the celltype and compartment
  mutate(celltype_compartment = paste(CellType,Compartment,sep = "|"))
rownames(colData) <- colData$PB_ID

if (MIN_CELLS > 0) {
  message("dropping ", sum(colData$N_cells < MIN_CELLS), " pseudobulks with < ", MIN_CELLS, " cells")
  colData <- colData[colData$N_cells >= MIN_CELLS, ]
}
counts <- counts[, rownames(colData)]

# wrangling ---------------------------------------------------------------
# design balance
colData %>%
  dplyr::count(CellType, Compartment) %>%
  pivot_wider(names_from = Compartment, values_from = n, values_fill = 0)

# which samples contribute to which cell type / compartment
colData %>%
  dplyr::count(Sample_Name, CellType, Compartment) %>%
  pivot_wider(names_from = c(CellType, Compartment), values_from = n, values_fill = 0) %>%
  print(n = Inf)

# "informative pairs" = samples with BOTH compartments within a cell type.
# This is what the blocked model actually uses to estimate the compartment effect.
colData %>%
  group_by(CellType, Sample_Name) %>%
  summarise(paired = n_distinct(Compartment) == 2, .groups = "drop") %>%
  group_by(CellType) %>%
  summarise(n_samples = n(), n_paired = sum(paired)) %>%
  print()

# design matrix ----------------------------------------------------------
# visualise the model matrix (Sample_Name as blocking factor)
design_block <- ~ Sample_Name + CellType * Compartment
mat_design <- model.matrix(design_block, colData)

# observation ordering built from the factors in colData (not from the dummy columns),
# so that reference-level rows (EC / Endoneurial) are grouped with their peers
LUT_sample <- colData %>%
  arrange(CellType, Compartment, Sample_Name) %>%
  mutate(Observation = PB_ID,
         Label = fct_inorder(paste(CellType, Compartment, Sample_Name, sep = "-"))) %>%
  dplyr::select(Observation, Label)

term_order <- colnames(mat_design) %>%
  str_replace("\\(Intercept\\)", "Int") %>%
  str_remove_all("CellType|Compartment|Sample_Name")

as.data.frame(mat_design) %>%
  mutate(Observation = rownames(mat_design)) %>%
  pivot_longer(cols = -Observation, names_to = "Term", values_to = "Value") %>%
  mutate(Term = str_replace(Term, "\\(Intercept\\)", "Int"),
         Term = str_remove_all(Term, "CellType|Compartment|Sample_Name"),
         Term = factor(Term, levels = term_order)) %>%
  left_join(LUT_sample, by = "Observation") %>%
  ggplot(aes(x = Term, y = Label, fill = as.factor(Value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  theme_minimal() +
  labs(title = "Design Matrix Visualization",
       subtitle = deparse(design_block),
       x = "Model Terms", y = "Observation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave(file.path(out_dir, "plot", "design_matrix.pdf"),
       width = 4 + 0.25 * ncol(mat_design), height = 8)

saveRDS(counts, file.path(out_dir, "object", "counts_EC_FIBRO_pBulk.rds"))
saveRDS(colData, file.path(out_dir, "object", "colData_EC_FIBRO_pBulk.rds"))

# filter low abundant genes ----
dds_base <- DESeqDataSetFromMatrix(countData = counts, colData = colData,
                                   design = ~ CellType * Compartment)
# group = the 4 CellType x Compartment combinations.
# filterByExpr keeps a gene only if it passes the CPM cutoff in at least as many samples as the smallest group (n, scaled down when > 10).
# Without group all samples count as one group and a gene would have to be expressed in ~70% of them, dropping genes specific to one cell type / compartment.
feat_keep <- edgeR::filterByExpr(counts(dds_base),
                                 group = colData$celltype_compartment)
table(feat_keep)
dds_base <- dds_base[feat_keep, ]

# QC ----------------------------------------------------------------------
colData %>%
  mutate(tot_reads = colSums(counts(dds_base))[rownames(.)]) %>%
  ggplot(aes(x = PB_ID, y = tot_reads, fill = CellType)) + geom_col() +
  theme_bw() +
  facet_wrap(~Compartment,scales = "free") +
  theme(axis.text.x = element_text(hjust = 1, angle = 90))
ggsave(file.path(out_dir, "plot", "libsize.pdf"), width = 10, height = 5)

vds <- vst(dds_base, blind = FALSE)

pdf(file.path(out_dir, "plot", "heatmap_SampleCluster.pdf"), width = 12, height = 7)
set.seed(1)
hm <- plot_sample_clustering(vds, anno_vars = c("CellType", "Compartment"),
                             distance = "euclidean")
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

pca <- plotPCA(vds, intgroup = c("CellType", "Compartment", "Sample_Name", "N_cells"),
               returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(x = PC1, y = PC2, col = Compartment, shape = CellType)) +
  geom_point(size = 3, alpha = 0.7) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw()
ggsave(file.path(out_dir, "plot", "PCA_byCompartment.pdf"), width = 6, height = 4)

# same PCA coloured by sample: shows how much the sample effect drives the variance
ggplot(pca, aes(x = PC1, y = PC2, col = Sample_Name, shape = CellType)) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = Sample_Name)) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw()
ggsave(file.path(out_dir, "plot", "PCA_bySample.pdf"), width = 12, height = 8)

# batch removal -----------------------------------------------------------
# remove the sample (batch) effect from the vst values, for visualisation only:
# it mimics on the plots what Sample_Name in the DE design does (DE uses the raw counts, never these values).
# The CellType * Compartment design is passed so that the donor effect is estimated together with the biology, as in the DE model. Without it removeBatchEffect subtracts each donor's mean profile and nothing else: donors with a single cell type (12a, 14 = EC only; 12b, 2, 5 = FIBRO only) lose the EC vs FIBRO difference too and collapse to the centre of PC1.
vds_batch <- vds

mat_vst_batch <- limma::removeBatchEffect(assay(vds),
                                          batch = vds$Sample_Name,
                                          design = model.matrix(~ CellType * Compartment, colData))
assay(vds_batch) <- mat_vst_batch

saveRDS(vds_batch, file.path(out_dir, "object", "vst_filterBatch_EC_FIBRO_pBulk.rds"))

# clustering samples
pdf(file.path(out_dir, "plot", "heatmap_SampleCluster_batch.pdf"), width = 12, height = 7)
set.seed(1)
hm_batch <- plot_sample_clustering(vds_batch, anno_vars = c("CellType", "Compartment", "Sample_Name"),
                                   distance = "euclidean")
draw(hm_batch, heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA
pca_batch <- plotPCA(vds_batch, intgroup = c("CellType", "Compartment", "Sample_Name", "N_cells"),
                     returnData = TRUE)
pv_batch <- round(100 * attr(pca_batch, "percentVar"))

p1_batch <- ggplot(pca_batch, aes(x = PC1, y = PC2, col = Compartment, shape = CellType)) +
  geom_point(size = 3, alpha = 0.7) +
  xlab(paste0("PC1: ", pv_batch[1], "% variance")) + ylab(paste0("PC2: ", pv_batch[2], "% variance")) +
  theme_bw()

# same PCA coloured by sample: the sample effect should be gone
p2_batch <- ggplot(pca_batch, aes(x = PC1, y = PC2, col = Sample_Name, shape = CellType)) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = Sample_Name)) +
  xlab(paste0("PC1: ", pv_batch[1], "% variance")) + ylab(paste0("PC2: ", pv_batch[2], "% variance")) +
  theme_bw()

p1_batch + p2_batch + plot_annotation("batch (Sample_Name) removal")
ggsave(file.path(out_dir, "plot", "PCA_batch.pdf"), width = 16, height = 8)

# compartment distance ----------------------------------------------------
# quantify the Endoneurial - Epineurial distance per cell type (visible in the PCA as a larger separation in FIBRO than in EC).
# For each donor having both compartments within a cell type, the euclidean distance between the two pseudobulks in the vst space (all filtered genes).
# Within-donor pairs cancel the donor effect, and the donor is the unit of the test.
dist_vst <- as.matrix(dist(t(assay(vds))))

df_dist <- colData %>%
  dplyr::select(PB_ID, Sample_Name, CellType, Compartment) %>%
  pivot_wider(names_from = Compartment, values_from = PB_ID) %>%
  filter(!is.na(Endoneurial), !is.na(Epineurial)) %>%
  mutate(distance = dist_vst[cbind(Endoneurial, Epineurial)])
write_tsv(df_dist, file.path(out_dir, "table", "distance_Endo_Epi_within_donor.tsv"))

# median distance per cell type and test (donors are different between EC and FIBRO)
df_dist %>%
  group_by(CellType) %>%
  summarise(n_donors = n(), median_distance = median(distance), .groups = "drop") %>%
  print()
print(wilcox.test(distance ~ CellType, data = df_dist))

ggplot(df_dist, aes(x = CellType, y = distance, col = CellType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  scale_color_manual(values = col_hm$CellType) +
  labs(y = "Endoneurial - Epineurial distance (vst, within donor)") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "plot", "boxplot_distance_Endo_Epi.pdf"), width = 4, height = 5)

# DE model ----------------------------------------------------------------
# the design must not be rank-deficient
stopifnot(qr(mat_design)$rank == ncol(mat_design))

dds <- dds_base
design(dds) <- design_block
dds <- DESeq(dds)
resultsNames(dds)
saveRDS(dds, file.path(out_dir, "object", "dds.rds"))

# coefficient names (treatment coding, reference = EC / Endoneurial)
cf_comp <- "Compartment_Epineurial_vs_Endoneurial"
cf_int <- "CellTypeFIBRO.CompartmentEpineurial"
stopifnot(all(c(cf_comp, cf_int) %in% resultsNames(dds)))

# contrasts as numeric vectors over resultsNames(dds)
contrasts_mat <- limma::makeContrasts(
  # effect in EC (reference level of CellType)
  Compartment_EC = Compartment_Epineurial_vs_Endoneurial,
  # effect in FIBRO
  Compartment_FIBRO = Compartment_Epineurial_vs_Endoneurial + CellTypeFIBRO.CompartmentEpineurial,
  # effect averaged over EC and FIBRO = effect after accounting for CellType
  Compartment_avg = Compartment_Epineurial_vs_Endoneurial + 0.5 * CellTypeFIBRO.CompartmentEpineurial,
  # does the compartment effect differ between FIBRO and EC?
  Interaction_FIBRO_vs_EC = CellTypeFIBRO.CompartmentEpineurial,
  levels = resultsNames(dds))

# one numeric vector per contrast, named by contrast
contrasts_list <- contrasts_mat %>%
  as.data.frame() %>%
  as.list()

res <- map(contrasts_list, function(x){
  results(dds, contrast = x, alpha = padj_thr)
})
map(res, summary)

res_shr <- map(res, function(x){
  lfcShrink(dds, res = x, type = "ashr")
})

# one long table: contrast x gene
df_all <- res %>%
  map(function(x){
    data.frame(x) %>%
      rownames_to_column("symbol") %>%
      arrange(pvalue)
  }) %>%
  bind_rows(.id = "contrast_name")

df_all_shr <- res_shr %>%
  map(function(x){
    data.frame(x) %>%
      rownames_to_column("symbol") %>%
      arrange(pvalue)
  }) %>%
  bind_rows(.id = "contrast_name")

write_tsv(df_all, file.path(out_dir, "table", "DE.tsv"))
write_tsv(df_all_shr, file.path(out_dir, "table", "DE_shr.tsv"))

pdf(file.path(out_dir, "plot", "dispersion.pdf"), width = 5, height = 5)
plotDispEsts(dds)
dev.off()

# summary of significant genes
df_all_shr %>%
  group_by(contrast_name) %>%
  summarise(n_padj = sum(padj < padj_thr, na.rm = TRUE),
            n_padj_lfc = sum(padj < padj_thr & abs(log2FoldChange) > lfc_thr, na.rm = TRUE),
            n_up = sum(padj < padj_thr & log2FoldChange > lfc_thr, na.rm = TRUE),
            n_down = sum(padj < padj_thr & log2FoldChange < -lfc_thr, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# wide gene table: lfc + padj for each contrast, plus a class from the averaged compartment
# effect and the interaction (significant = padj < padj_thr and |LFC| > lfc_thr)
df_wide <- df_all_shr %>%
  dplyr::select(symbol, contrast_name, log2FoldChange, padj) %>%
  pivot_wider(names_from = contrast_name, values_from = c(log2FoldChange, padj),
              names_glue = "{contrast_name}__{.value}") %>%
  # padj = NA (e.g. independent filtering) falls through to FALSE
  mutate(sig_comp = case_when(Compartment_avg__padj < padj_thr & abs(Compartment_avg__log2FoldChange) > lfc_thr ~ TRUE,
                              TRUE ~ FALSE),
         sig_int = case_when(Interaction_FIBRO_vs_EC__padj < padj_thr & abs(Interaction_FIBRO_vs_EC__log2FoldChange) > lfc_thr ~ TRUE,
                             TRUE ~ FALSE),
         class = case_when(sig_comp & sig_int ~ "compartment+interaction",
                           sig_comp ~ "compartment_only",
                           sig_int ~ "interaction_only",
                           TRUE ~ "ns")) %>%
  dplyr::select(-sig_comp, -sig_int)
write_tsv(df_wide, file.path(out_dir, "table", "DE_gene_classification.tsv"))
dplyr::count(df_wide, class)

# EC vs FIBRO compartment effect ------------------------------------------
# genes DE (padj < padj_thr and |LFC| > lfc_thr) in EC and in FIBRO, split by direction, compared with the averaged contrast.
# Compartment_avg is by construction (Compartment_EC + Compartment_FIBRO) / 2, so the LFCs are
# identical; the DEG sets are not, because "avg" tests the mean = 0 while the intersection
# requires both to be significant. Unshrunk LFC are used (shrinkage is not linear).
df_dir <- df_all %>%
  filter(contrast_name %in% c("Compartment_EC", "Compartment_FIBRO", "Compartment_avg")) %>%
  dplyr::select(symbol, contrast_name, log2FoldChange, padj) %>%
  pivot_wider(names_from = contrast_name, values_from = c(log2FoldChange, padj)) %>%
  # padj = NA falls through to FALSE
  mutate(EC_up = case_when(padj_Compartment_EC < padj_thr & log2FoldChange_Compartment_EC > lfc_thr ~ TRUE,
                           TRUE ~ FALSE),
         EC_down = case_when(padj_Compartment_EC < padj_thr & log2FoldChange_Compartment_EC < -lfc_thr ~ TRUE,
                             TRUE ~ FALSE),
         FIBRO_up = case_when(padj_Compartment_FIBRO < padj_thr & log2FoldChange_Compartment_FIBRO > lfc_thr ~ TRUE,
                              TRUE ~ FALSE),
         FIBRO_down = case_when(padj_Compartment_FIBRO < padj_thr & log2FoldChange_Compartment_FIBRO < -lfc_thr ~ TRUE,
                                TRUE ~ FALSE),
         avg_up = case_when(padj_Compartment_avg < padj_thr & log2FoldChange_Compartment_avg > lfc_thr ~ TRUE,
                            TRUE ~ FALSE),
         avg_down = case_when(padj_Compartment_avg < padj_thr & log2FoldChange_Compartment_avg < -lfc_thr ~ TRUE,
                              TRUE ~ FALSE),
         overlap = case_when(EC_up & FIBRO_up ~ "up_both",
                             EC_down & FIBRO_down ~ "down_both",
                             (EC_up | EC_down) & (FIBRO_up | FIBRO_down) ~ "discordant",
                             EC_up | EC_down ~ "EC_only",
                             FIBRO_up | FIBRO_down ~ "FIBRO_only",
                             TRUE ~ "ns"),
         avg_sig = avg_up | avg_down)
write_tsv(df_dir, file.path(out_dir, "table", "DE_compartment_EC_vs_FIBRO_overlap.tsv"))

# 1. the LFC of the average contrast is the mean of the two: the difference must be ~0
df_dir %>%
  summarise(max_abs_diff_lfc = max(abs(log2FoldChange_Compartment_avg -
                                         (log2FoldChange_Compartment_EC + log2FoldChange_Compartment_FIBRO) / 2),
                                   na.rm = TRUE)) %>%
  print()

# 2. how do the DEG in both cell types relate to the DEG of the averaged contrast?
dplyr::count(df_dir, overlap, avg_sig) %>% print()

# 3. LFC in EC vs LFC in FIBRO, coloured by the overlap class
df_dir %>%
  filter(!is.na(padj_Compartment_EC), !is.na(padj_Compartment_FIBRO)) %>%
  arrange(overlap != "ns") %>%
  ggplot(aes(x = log2FoldChange_Compartment_EC, y = log2FoldChange_Compartment_FIBRO, col = overlap)) +
  geom_point(alpha = 0.4) +
  geom_abline(col = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, col = "gray") + geom_vline(xintercept = 0, col = "gray") +
  scale_color_manual(values = c(ns = "gray80", EC_only = "firebrick", FIBRO_only = "steelblue",
                                up_both = "darkgreen", down_both = "purple", discordant = "orange")) +
  labs(x = "LFC Epi vs Endo, EC", y = "LFC Epi vs Endo, FIBRO") +
  theme_bw()
ggsave(file.path(out_dir, "plot", "scatter_LFC_Compartment_EC_vs_FIBRO.pdf"), width = 7, height = 6)

# 4. upset of the compartment DEG sets (EC and FIBRO) split by direction, without the averaged
# contrast: common genes = intersection of same-direction sets, set size = number of genes
sets_dir <- c("EC_up", "EC_down", "FIBRO_up", "FIBRO_down")
df_dir %>%
  filter(if_any(all_of(sets_dir), function(x) x)) %>%
  as.data.frame() %>%
  ComplexUpset::upset(intersect = sets_dir, name = "Epi vs Endo DEG", width_ratio = 0.4,
                      set_sizes = ComplexUpset::upset_set_size() +
                        geom_text(aes(label = after_stat(count)), hjust = -0.1, colour = "white", stat = "count"))
ggsave(file.path(out_dir, "plot", "upset_DEG_Compartment_direction.pdf"), width = 10, height = 6)

# diagnostics -------------------------------------------------------------
# p-value histograms (genes with baseMean > 1)
df_all %>%
  filter(baseMean > 1) %>%
  ggplot(aes(x = pvalue)) + geom_histogram(breaks = 0:20 / 20) +
  facet_wrap(~contrast_name) +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(out_dir, "plot", "histogram_pvalue.pdf"), width = 8, height = 6)

# volcano
# top n_lab significant genes per contrast are labelled
n_lab <- 15

list(volcano = df_all, volcano_shr = df_all_shr) %>%
  map(function(df) {
    df_volcano <- df %>%
      filter(!is.na(padj)) %>%
      mutate(sig = padj < padj_thr & abs(log2FoldChange) > lfc_thr)

    ggplot(df_volcano, aes(x = log2FoldChange, y = -log10(padj), col = sig)) +
      geom_point(aes(alpha = sig)) +
      scale_alpha_manual(values = c(`FALSE` = 0.1, `TRUE` = 0.5)) +
      scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), col = "red", linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_thr), col = "red", linetype = "dashed") +
      geom_text_repel(data = df_volcano %>% filter(sig) %>% group_by(contrast_name) %>%
                        arrange(padj) %>% dplyr::slice(1:n_lab),
                      aes(label = symbol), segment.alpha = 0.4, show.legend = FALSE) +
      facet_wrap(~contrast_name) +
      theme_bw() +
      theme(strip.background = element_blank(), legend.position = "none")
  }) %>%
  iwalk(function(p, tag) {
    ggsave(file.path(out_dir, "plot", sprintf("%s.pdf", tag)), p, width = 10, height = 10)
  })

# MA plot
df_all %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj < padj_thr & abs(log2FoldChange) > lfc_thr) %>%
  arrange(sig) %>%
  ggplot(aes(baseMean, log2FoldChange, col = sig)) +
  geom_point(alpha = 0.2) +
  geom_hline(yintercept = c(-lfc_thr, lfc_thr), col = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c(`FALSE` = "gray", `TRUE` = "red")) +
  facet_wrap(~contrast_name) +
  theme_bw() + theme(strip.background = element_blank(), legend.position = "none")
ggsave(file.path(out_dir, "plot", "MA.pdf"), width = 10, height = 8)

# overlap of the significant gene sets across contrasts
# one row per gene, one logical column per significant set (ComplexUpset input)
df_upset <- df_all_shr %>%
  filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr) %>%
  mutate(present = TRUE) %>%
  distinct(symbol, contrast_name, present) %>%
  pivot_wider(names_from = contrast_name, values_from = present, values_fill = FALSE) %>%
  as.data.frame()
sets_upset <- setdiff(colnames(df_upset), "symbol")

if (length(sets_upset) > 1) {
  ComplexUpset::upset(df_upset, intersect = sets_upset,
                      name = "contrast", width_ratio = 0.25)
  ggsave(file.path(out_dir, "plot", "upset_DEG.pdf"), width = 12, height = 6)
}

# heatmaps ----------------------------------------------------------------
mat_vst <- assay(vds)
col_order <- colData %>% arrange(CellType, Compartment, Sample_Name) %>% rownames()

ann_hm <- HeatmapAnnotation(CellType = colData[col_order, "CellType"],
                            Compartment = colData[col_order, "Compartment"],
                            col = col_hm)

# top genes per contrast (heatmaps with < 2 genes are skipped)
# the cell type specific contrasts are drawn on the samples of that cell type only
df_top <- df_all_shr %>%
  filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr) %>%
  group_by(contrast_name) %>%
  arrange(padj) %>%
  dplyr::slice(1:max_hm_genes) %>%
  summarise(genes = list(symbol), .groups = "drop") %>%
  filter(lengths(genes) >= 2) %>%
  mutate(celltype = case_when(contrast_name == "Compartment_EC" ~ "EC",
                              contrast_name == "Compartment_FIBRO" ~ "FIBRO",
                              TRUE ~ "all"),
         title = paste(contrast_name, "(Epi vs Endo, adj. for Sample_Name)"),
         file = file.path(out_dir, "plot", sprintf("heatmap_DEG_%s.pdf", contrast_name)))

pwalk(df_top, function(genes, celltype, title, file, ...) {
  # samples to show, in the usual order
  cols <- col_order[celltype == "all" | colData[col_order, "CellType"] == celltype]

  # z-score over the shown samples only
  m <- t(scale(t(mat_vst[genes, cols, drop = FALSE])))
  ht <- Heatmap(m, name = "z-score vst", column_title = title,
                top_annotation = HeatmapAnnotation(CellType = colData[cols, "CellType"],
                                                   Compartment = colData[cols, "Compartment"],
                                                   col = col_hm),
                column_split = colData[cols, "CellType"],
                cluster_columns = FALSE, show_row_dend = FALSE,
                show_column_names = FALSE,
                show_row_names = length(genes) <= 80,
                use_raster = TRUE)
  pdf(file, width = 8, height = 12)
  draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")
  dev.off()
})

# genes DE in both EC and FIBRO with the same direction (see "EC vs FIBRO compartment effect")
df_conc <- df_dir %>%
  filter(overlap %in% c("up_both", "down_both")) %>%
  arrange(overlap)

m_conc <- t(scale(t(mat_vst[df_conc$symbol, col_order, drop = FALSE])))
ht_conc <- Heatmap(m_conc, name = "z-score vst",
                   column_title = "Epi vs Endo DEG common to EC and FIBRO (same direction)",
                   top_annotation = ann_hm,
                   column_split = colData[col_order, "CellType"],
                   row_split = df_conc$overlap,
                   cluster_columns = FALSE, show_row_dend = FALSE,
                   show_column_names = FALSE, show_row_names = FALSE,
                   use_raster = TRUE)
pdf(file.path(out_dir, "plot", "heatmap_DEG_Compartment_common_EC_FIBRO.pdf"), width = 8, height = 12)
draw(ht_conc, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

# same genes, scaled separately in EC and in FIBRO (z-score over the samples of that cell type
# only), and the two heatmaps joined side by side. Rows are sorted by the averaged LFC within
# up_both / down_both so that the same row is the same gene in both heatmaps.
df_conc_sorted <- df_conc %>%
  arrange(overlap, log2FoldChange_Compartment_avg)

cols_ec <- col_order[colData[col_order, "CellType"] == "EC"]
cols_fibro <- col_order[colData[col_order, "CellType"] == "FIBRO"]

m_conc_ec <- t(scale(t(mat_vst[df_conc_sorted$symbol, cols_ec, drop = FALSE])))
m_conc_fibro <- t(scale(t(mat_vst[df_conc_sorted$symbol, cols_fibro, drop = FALSE])))

# same colour scale for the two heatmaps
col_z <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

ht_conc_ec <- Heatmap(m_conc_ec, name = "z-score vst", col = col_z,
                      column_title = "EC",
                      top_annotation = HeatmapAnnotation(Compartment = colData[cols_ec, "Compartment"],
                                                         col = col_hm["Compartment"]),
                      row_split = df_conc_sorted$overlap,
                      cluster_rows = T,
                      cluster_columns = T,
                      show_column_names = FALSE,
                      show_row_names = FALSE,
                      use_raster = TRUE)
ht_conc_fibro <- Heatmap(m_conc_fibro, name = "z-score vst FIBRO", col = col_z,
                         column_title = "FIBRO",
                         top_annotation = HeatmapAnnotation(Compartment = colData[cols_fibro, "Compartment"],
                                                            col = col_hm["Compartment"],
                                                            show_legend = FALSE),
                         cluster_rows = T,
                         cluster_columns = T,
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         show_heatmap_legend = FALSE,
                         use_raster = TRUE)

pdf(file.path(out_dir, "plot", "heatmap_DEG_Compartment_common_EC_FIBRO_scaledByCellType.pdf"), width = 10, height = 12)
draw(ht_conc_ec + ht_conc_fibro, column_title = "Epi vs Endo DEG common to EC and FIBRO (same direction)",
     heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

# session -----------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "table", "sessionInfo.txt"))
