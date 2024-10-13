library(Seurat)
library(tidyverse)

seu <- readRDS("rds/manual_annoted.rds")
head(seu);table(seu$celltype)
DimPlot(seu, group.by = "celltype", label = T)
DimPlot(seu, group.by = "seurat_clusters", label = T)

Idents(seu) <- seu$celltype
deg_20 <- FindMarkers(seu, ident.1 = 20, group.by = "seurat_clusters")
deg_21 <- FindMarkers(seu, ident.1 = 21, group.by = "seurat_clusters")
deg_22 <- FindMarkers(seu, ident.1 = 22, group.by = "seurat_clusters")
write.csv(deg_22,"deg_result/unknwon_22.csv")

seu$celltype <- ifelse(seu$seurat_clusters %in% 25, "Macrophages", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 20, "Immune cells", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 22, "Astrocytes", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 21, "Fibroblasts", seu$celltype)


# Load the necessary libraries
library(Seurat)
library(ggplot2)
library(ggsci)  # for ggsci color palettes

# Reorder the cell types: Microglia first, Neuron last
seu$celltype <- factor(seu$celltype, levels = c("Microglia", "Oligodendrocytes", "Astrocytes", 
                                                "Endothelial Cells", "Ependymal Cells", "VSMCs", 
                                                "Macrophages", "Fibroblasts", "Immune cells", "Neuron"))

# Specify the palette you are using
chosen_palette <- "npg"  # Change this to the palette you're using, e.g., "jco", "d3", etc.

# Create the UMAP plot with ggsci palette and enhanced aesthetics
DimPlot(seu, group.by = "celltype", label = TRUE) + 
  # Use a palette from ggsci, adjust based on the number of cell types
  scale_color_npg() +  # Choose from palettes like "npg", "d3", "jco" depending on your preference
  
  # Set theme for smaller axis labels and better overall appearance
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),       # Remove axis titles
    axis.text = element_blank(), # Smaller axis text
    axis.ticks = element_blank(),       # Remove ticks
    legend.position = "right",         # Position the legend at the bottom
    legend.title = element_blank(),     # Remove legend title
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  
  # Adjust the label size, font, and positioning
  # geom_text(aes(label = celltype), size = 5, fontface = "bold", vjust = 1.5, hjust = 1.5) +
  
  # Ensure that labels don't overlap, good for large numbers of clusters
  guides(label = guide_legend(nrow = 1, byrow = TRUE)) +

  # Set the title and improve readability, including the palette name
  labs(title = paste("UMAP of Mouse Brain Cells - Palette:", chosen_palette))



