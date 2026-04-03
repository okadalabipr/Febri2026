# 04_network_analysis.R
# Description:
# Visualization of signaling network with gene variability (Fano factor) mapped onto nodes across conditions.

# Load libraries
library(igraph)
library(ggraph)
library(tidyverse)

# Define paths
fano_path <- "results/fano/fano_values.rds"
output_dir <- "results/network/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Fano data
fano_df <- readRDS(fano_path)

# Define signaling network
edges <- tribble(
  ~from, ~to,
  "RTK", "MAPK3",
  "MAPK3", "AKT1",
  "AKT1", "MYC",
  "AKT1", "FOS",
  "FOS", "CCND1",
  "MYC", "CCND1",
  "CCND1", "CDK4",
  "CCND1", "CDK6",
  "CDK4", "RB1",
  "RB1", "E2F1",
  "E2F1", "CCNE1", 
  "E2F1", "CCNA1",
  "CCNE1", "CDK2"
)

g <- graph_from_data_frame(edges, directed = TRUE)

# Select timepoint
timepoint_use <- 8 #4, 16

# Prepare node attributes
safe_log <- function(x, eps = 1e-6) log10(x + eps)

fano_net <- fano_df %>%
  filter(time_h == timepoint_use) %>%
  select(gene, ligand, Fano) %>%
  pivot_wider(names_from = ligand, values_from = Fano)

node_df <- data.frame(name = V(g)$name) %>%
  left_join(fano_net, by = c("name" = "gene")) %>%
  mutate(
    EGF_log = safe_log(EGF),
    HRG_log = safe_log(HRG),
    delta = HRG_log - EGF_log
  )

# Assign values to graph
V(g)$EGF <- node_df$EGF_log
V(g)$HRG <- node_df$HRG_log
V(g)$delta <- node_df$delta

# Global color scaling
all_vals <- c(node_df$EGF_log, node_df$HRG_log)

global_min <- min(all_vals, na.rm = TRUE)
global_max <- max(all_vals, na.rm = TRUE)
global_mid <- median(all_vals, na.rm = TRUE)

# Plot function
plot_network <- function(graph, value_col, title = "", is_delta = FALSE) {
  
  p <- ggraph(graph, layout = "tree") +
    geom_edge_link(
      arrow = arrow(length = unit(4, "mm")),
      end_cap = circle(3, "mm"),
      alpha = 0.6
    ) +
    geom_node_label(
      aes(label = name, fill = .data[[value_col]]),
      color = "black",
      size = 5
    ) +
    theme_void() +
    ggtitle(title)
  
  if (is_delta) {
    p <- p +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        name = "Δlog10(Fano)",
        na.value = "grey"
      )
  } else {
    p <- p +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = global_mid,
        limits = c(global_min, global_max),
        name = "log10(Fano)",
        na.value = "white"
      )
  }
  
  return(p)
}

# Generate plots
p_egf <- plot_network(g, "EGF", paste0("EGF (", timepoint_use, "h)"))
p_hrg <- plot_network(g, "HRG", paste0("HRG (", timepoint_use, "h)"))
p_delta <- plot_network(g, "delta", "ΔFano (HRG - EGF)", is_delta = TRUE)

# Save outputs
ggsave(file.path(output_dir, "network_EGF.png"),
       p_egf, width = 6, height = 5)

ggsave(file.path(output_dir, "network_HRG.png"),
       p_hrg, width = 6, height = 5)

# Combined plot
library(patchwork)

p_combined <- p_egf | p_hrg

ggsave(file.path(output_dir, "network_combined.png"),
       p_combined, width = 10, height = 5)
