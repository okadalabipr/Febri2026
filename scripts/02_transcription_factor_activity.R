# 02_transcription_factor_activity.R
# Description:
# Transcription factor activity inference using decoupleR with DoRothEA regulons, followed by visualization and statistical comparison between EGF and HRG conditions.

# Load libraries
library(Seurat)
library(decoupleR)
library(dorothea)
library(tidyverse)
library(pheatmap)
library(viridisLite)
library(rstatix)
library(OmnipathR)

set.seed(1)

# Define paths
data_path <- "data/stimulated_1.rds"
output_dir <- "results/tf_activity/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
seurat_obj <- readRDS(data_path)

# Extract expression matrix
expression_matrix <- GetAssayData(seurat_obj, layer = "data", assay = "RNA")

# Load DoRothEA regulon
dorothea_db <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# TF activity inference (ULM)
acts <- run_ulm(
  mat = expression_matrix,
  net = dorothea_db,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  minsize = 5
)

# Store TF activity in Seurat object
tf_mat <- acts %>%
  pivot_wider(id_cols = "source",
              names_from = "condition",
              values_from = "score") %>%
  column_to_rownames("source")

seurat_obj[["tfsulm"]] <- CreateAssayObject(tf_mat)
DefaultAssay(seurat_obj) <- "tfsulm"

seurat_obj <- ScaleData(seurat_obj)

# Aggregate TF activity per condition
df <- t(as.matrix(seurat_obj@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(condition = seurat_obj$orig.ident) %>%
  pivot_longer(cols = -condition,
               names_to = "TF",
               values_to = "score") %>%
  group_by(condition, TF) %>%
  summarise(mean = mean(score), .groups = "drop")

# Select TFs of interest
selected_tfs <- c("E2F1","E2F2","E2F3","MYC","FOXM1","RB1","TFDP1",
                  "FOSL1","FOSB","FOS","JUN","SP1","CREB1","STAT3",
                  "SMAD3","NFKB1")

df_filtered <- df %>%
  filter(TF %in% selected_tfs)

# Heatmap visualization
tf_matrix <- df_filtered %>%
  pivot_wider(id_cols = condition,
              names_from = TF,
              values_from = mean) %>%
  column_to_rownames("condition") %>%
  as.matrix()

# Order conditions
order <- c("Control_0h","EGF_4h","EGF_8h","EGF_16h",
           "HRG_4h","HRG_8h","HRG_16h")

tf_matrix <- tf_matrix[order, ]
tf_matrix <- t(tf_matrix)

pheatmap(
  tf_matrix,
  color = viridis(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE
)

# Statistical comparison (EGF vs HRG)
tf_scores <- t(as.matrix(seurat_obj@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(condition = seurat_obj$orig.ident) %>%
  rownames_to_column("cell")

tf_scores_long <- tf_scores %>%
  pivot_longer(cols = "SP1", names_to = "TF", values_to = "score")

tf_scores_long <- tf_scores_long %>%
  mutate(
    ligand = ifelse(grepl("^EGF", condition), "EGF", "HRG"),
    time_h = as.numeric(sub(".*_(\\d+)h", "\\1", condition))
  )

# Statistical test
p_df <- tf_scores_long %>%
  group_by(time_h) %>%
  wilcox_test(score ~ ligand) %>%
  adjust_pvalue(method = "BH")

# Plot
p <- ggplot(tf_scores_long,
            aes(x = factor(time_h), y = score, color = ligand)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.3) +
  scale_color_manual(values = c(EGF = "blue", HRG = "red")) +
  theme_bw()

ggsave(file.path(output_dir, "SP1_activity.png"),
       p, width = 6, height = 5)

# Summary statistics
n_fos_targets <- dorothea_db %>%
  filter(tf == "FOS") %>%
  distinct(target) %>%
  nrow()

n_myc_targets <- dorothea_db %>%
  filter(tf == "MYC") %>%
  distinct(target) %>%
  nrow()

# Save object
saveRDS(seurat_obj,
        file = file.path(output_dir, "tf_activity_seurat.rds"))