# 03_fano_factor.R
# Description:
# Analysis of gene expression variability using Fano factor across EGF and HRG conditions, focusing on cell-cycle genes.

# Load libraries
library(Seurat)
library(tidyverse)
library(msigdbr)
library(dorothea)
library(ggtext)

set.seed(1)

# Define paths
data_path <- "data/downsampled_seurat.rds"
output_dir <- "results/fano/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
seurat_obj <- readRDS(data_path)

DefaultAssay(seurat_obj) <- "RNA"

expr_data <- as.matrix(GetAssayData(seurat_obj, layer = "data"))

# Get cell-cycle genes
go_cc <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  filter(grepl("CELL_CYCLE", gs_name, ignore.case = TRUE)) %>%
  select((gs_name, gene_symbol)) %>%
  distinct()

genes_cc <- intersect(unique(go_cc$gene_symbol), rownames(expr_data))

# Fano factor function
calculate_fano <- function(x) {
  m <- mean(x, na.rm = TRUE)
  v <- var(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  v / m
}

# Compute Fano per condition
Idents(seurat_obj) <- "orig.ident"

conditions <- levels(Idents(seurat_obj))
conditions <- grep("^(EGF|HRG)_(4|8|16)h$", conditions, value = TRUE)

fano_df <- map_dfr(conditions, function(cond) {
  cells <- WhichCells(seurat_obj, idents = cond)
  expr <- expr_data[genes_cc, cells, drop = FALSE]
  
  fano <- apply(expr, 1, calculate_fano)
  
  tibble(
    gene = names(fano),
    Fano = as.numeric(fano),
    condition = cond
  )
}) %>%
  filter(is.finite(Fano))

# Add metadata
fano_df <- fano_df %>%
  mutate(
    ligand = ifelse(str_detect(condition, "^EGF"), "EGF", "HRG"),
    time_h = as.numeric(str_extract(condition, "\\d+"))
  )

# Summary
boot_ci <- function(x, B = 1000, conf = 0.95) {
  meds <- replicate(B, median(sample(x, replace = TRUE)))
  tibble(
    med = median(x),
    lo = quantile(meds, 0.025),
    hi = quantile(meds, 0.975)
  )
}

fano_summary <- fano_df %>%
  group_by(ligand, time_h) %>%
  summarise(boot_ci(Fano), .groups = "drop")

# Plot global trend
p_global <- ggplot(fano_summary,
                   aes(x = time_h, y = med, color = ligand)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  scale_color_manual(values = c(EGF = "blue", HRG = "red")) +
  theme_classic()

ggsave(file.path(output_dir, "fano_global.png"),
       p_global, width = 6, height = 5)

# Gene-level comparison
safe_log <- function(x, eps = 1e-6) log10(x + eps)

plot_fano_scatter <- function(timepoint) {
  df <- fano_df %>%
    filter(time_h == timepoint) %>%
    pivot_wider(names_from = ligand, values_from = Fano) %>%
    filter(EGF > 0, HRG > 0) %>%
    mutate(delta = safe_log(HRG) - safe_log(EGF))
  
  ggplot(df, aes(x = safe_log(EGF), y = safe_log(HRG), color = delta)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_gradient2(low = "blue", high = "red") +
    theme_classic()
}

p_16h <- plot_fano_scatter(16)

ggsave(file.path(output_dir, "fano_scatter_16h.png"),
       p_16h, width = 6, height = 5)

# Distribution
delta_all <- fano_df %>%
  pivot_wider(names_from = ligand, values_from = Fano) %>%
  mutate(delta = safe_log(HRG) - safe_log(EGF))

p_dist <- ggplot(delta_all, aes(x = delta)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = 0, linetype = "dashed")

ggsave(file.path(output_dir, "fano_distribution.png"),
       p_dist)

# TF-target variability
dorothea_db <- dorothea_hs %>%
  filter(confidence %in% c("A","B","C"))

myc_targets <- dorothea_db %>% filter(tf == "MYC") %>% pull(target)
fos_targets <- dorothea_db %>% filter(tf == "FOS") %>% pull(target)

tf_fano <- fano_df %>%
  filter(gene %in% c(myc_targets, fos_targets)) %>%
  mutate(TF = case_when(
    gene %in% myc_targets ~ "MYC",
    gene %in% fos_targets ~ "FOS"
  ))

p_tf <- ggplot(tf_fano, aes(x = factor(time_h), y = Fano, fill = ligand)) +
  geom_boxplot() +
  facet_wrap(~TF) +
  scale_fill_manual(values = c(EGF = "blue", HRG = "red"))

ggsave(file.path(output_dir, "fano_tf_targets.png"),
       p_tf, width = 8, height = 5)

# Save results
saveRDS(fano_df, file = file.path(output_dir, "fano_values.rds"))
