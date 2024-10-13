# Example of running FindMarkers to compare two groups (group1 vs group2)
table(seu$group)
markers <- FindMarkers(seu, ident.1 = "AD + WT", ident.2 = "Control", group.by = "group")
# Load necessary libraries
library(ggplot2)

# Assume the output of FindMarkers is stored in 'markers' dataframe
# Add a column for gene names (you can customize this to match your dataset)
markers$gene <- rownames(markers)

# Create a new column to label significant points and highlight the genes of interest
markers$highlight <- ifelse(markers$gene %in% c("Mertk"), "highlight", "normal")

# Create a new column to categorize based on adjusted p-value and log fold change
markers$significance <- ifelse(markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25, "Significant", "Not Significant")

# Plot the volcano plot
ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  
  # Highlight specific genes with different color and size
  geom_point(data = subset(markers, highlight == "highlight"), aes(color = gene), size = 4) +
  
  # Add labels for the highlighted genes
  geom_text_repel(data = subset(markers, highlight == "highlight"), size = 5, fontface = "bold", 
                  box.padding = 0.5, point.padding = 0.3, segment.color = "black") +
  
  # Set colors: adjust as needed to match the reference image
  scale_color_manual(values = c("grey", "red", "blue")) +
  
  # Customize plot aesthetics to match the reference image
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  
  # Add title and axis labels
  labs(title = "Volcano Plot: Group1 vs Group2", 
       x = "log2 Fold Change", 
       y = "-log10 Adjusted P-value") +
  
  # Optional: add vertical and horizontal lines to indicate thresholds
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# Loss: Seurat's VlnPlot() function, by default, does not allow split violins for more than two groups
# VlnPlot(seu, "Mertk", group.by = "group", split.by = "celltype", split.plot = TRUE)
# Load necessary libraries
library(ggplot2)
library(ggsci)
library(patchwork)

# Create a sample data frame for plotting
df <- data.frame(
  x = factor(rep(1:3, each = 3)),
  y = c(2, 4, 6, 2, 4, 6, 2, 4, 6),
  group = factor(rep(1:3, 3))
)

# Define a list of ggsci palettes to display (first 3 colors for each palette)
palettes <- list(
  "npg" = scale_fill_npg(),
  "d3" = scale_fill_d3(),
  "jco" = scale_fill_jco(),
  "lancet" = scale_fill_lancet(),
  "uchicago" = scale_fill_uchicago(),
  "startrek" = scale_fill_startrek(),
  "tron" = scale_fill_tron(),
  "futurama" = scale_fill_futurama()
)

# Generate individual plots for each palette
plots <- lapply(names(palettes), function(palette_name) {
  ggplot(df, aes(x = x, y = y, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    palettes[[palette_name]] +
    ggtitle(palette_name) +
    theme_minimal() +
    theme(legend.position = "none")
})

# Combine the plots into a single layout
combined_plots <- wrap_plots(plots, ncol = 4)

# Show the combined plot
combined_plots

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggsci)
library(patchwork)

# Define the cell types you're interested in
celltypes <- c("Microglia", "Neuron", "Immune cells", "Astrocytes", "Oligodendrocytes", "Endothelial Cells")

# Create individual violin plots for each cell type
vln_plots <- list()
for (celltype in celltypes) {
  p <- VlnPlot(seu, features = "Mertk", group.by = "group", 
               split.by = "group", pt.size = 0.1, split.plot = FALSE) +
    
    # Set the plot title as the celltype
    ggtitle(celltype) +
    
    # Leave x-axis label blank
    xlab("") +
    
    # Customize the plot appearance
    theme_minimal(base_line_size = 0) +  # Remove base grid lines
    theme(
      panel.grid = element_blank(),      # Remove all grid lines
      axis.line = element_line(),        # Keep axis lines
      axis.ticks = element_line(),       # Keep ticks
      legend.position = "none"           # Remove legend
    ) +
    
    # Apply the d3 color palette to the 'group' factor and color the points with the same palette
    scale_fill_d3() +
    
    # Modify the jitter to match the fill color
    geom_jitter(aes(color = group), size = 0.1, width = 0.1) +
    scale_color_d3()  # Use the same d3 color palette for the points
  
  if (T) {print(p)}
  vln_plots[[celltype]] <- p
}
# Combine the plots using patchwork (ncol = 2)
combined_plot <- wrap_plots(vln_plots, ncol = 2)

# Show the final combined plot
print(combined_plot)




# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggsci)
library(patchwork)

# Define the cell types you're interested in
celltypes <- c("Microglia", "Neuron", "Immune cells", "Astrocytes", "Oligodendrocytes", "Endothelial Cells")

# Create individual violin plots for each cell type
vln_plots <- list()
for (celltype in celltypes) {
  
  # Extract the data for plotting
  data_to_plot <- FetchData(seu, vars = c("Mertk", "group", "celltype"))
  
  # Filter for the current celltype
  data_filtered <- subset(data_to_plot, celltype == celltype)
  
  # Create the violin plot using ggplot2
  p <- ggplot(data_filtered, aes(x = group, y = Mertk, fill = group)) +
    geom_violin(scale = "width", trim = FALSE) +
    
    # Add jitter points matching the violin colors
    geom_jitter(aes(color = group), size = 0.01, alpha = 0.1, width = 0.03) +
    
    # Set the plot title as the celltype
    ggtitle(celltype) +
    
    # Leave x-axis label blank
    xlab("") +
    
    # Customize the plot appearance
    theme_minimal(base_line_size = 0) +  # Remove base grid lines
    theme(
      panel.grid = element_blank(),      # Remove all grid lines
      axis.line = element_line(),        # Keep axis lines
      axis.ticks = element_line(),       # Keep ticks
      legend.position = "none"           # Remove legend
    ) +
    
    # Apply the d3 color palette to both the violins and the jitter points
    scale_fill_d3() +
    scale_color_d3()
  p
  vln_plots[[celltype]] <- p
}

# Combine the plots using patchwork (ncol = 2)
combined_plot <- wrap_plots(vln_plots, ncol = 2)

# Show the final combined plot
print(combined_plot)

pdf("visual/mertk_vln.pdf" ,width =8, height = 10)
print(combined_plot)
dev.off()

