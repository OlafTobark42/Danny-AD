fig_dir <- "visual/"


# SpatialDimPlot ----------------------------------------------------------


dir.create(paste0(fig_dir, 'human_spatial_featureplots'))

selected_samples <- c("Dec_13_2021_Human3", "Dec_13_2021_Human4", "Dec_13_2021_Human5", "Dec_13_2021_Human6")
seurat_subset <- subset(seurat_obj, Sample %in% selected_samples)


barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
slice_obj <- subset(seurat_obj, cells = barcodes)

######################################################
# dimplot
######################################################

plot_df <- seurat_obj@meta.data

plot_list <- list()
for(cur_sample in selected_samples){
  print(cur_sample)
  cur_df <- plot_df %>% subset(Sample == cur_sample)
  cur_dx <- unique(cur_df$Diagnosis)
  plot_list[[cur_sample]] <- cur_df  %>%
    ggplot(aes(x=imagerow, y=imagecol, color=annotation)) +
    geom_point(size=0.5) +
    # umap_theme +
    ggtitle(cur_dx) +
    theme(plot.title=element_text(hjust=0.5))+
}

patch <- wrap_plots(plot_list, ncol=2) +
  plot_layout(guides='collect')

pdf(paste0(fig_dir, 'rep_human_annotations.pdf'), width=6, height=6)
print(patch)
dev.off()



# SpatialFeaturePlot ------------------------------------------------------

######################################################
# featureplot
######################################################
plot_genes <- c("CX3CR1","P2RY12","TMEM119","MERTK")
plot_genes %in% rownames(seurat_obj)

barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
slice_obj <- subset(seurat_obj, cells = barcodes)

# test
cur_gene <- "MERTK"
print(cur_gene)
p <- SpatialFeaturePlot(
  object=slice_obj,
  feature=cur_gene, images = slice,
  #samples_to_plot = selected_samples,
  # sample_labels = c("Diagnosis"),  # , "Sex", 'Age'
  # ncol = 10,
  # raster=TRUE,
  # plot_max = 'q99',
  # plot_min = 'q0',
  # colfunc = rainbow,
  # rev_colors=TRUE,
  # dpi=400
  )+theme_minimal()+
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()
  )

p

for(cur_gene in plot_genes){
  print(cur_gene)
  p <- SampleFeaturePlot(
    slice_obj,
    feature=cur_gene,
    #samples_to_plot = selected_samples,
    sample_labels = c("Diagnosis"),  # , "Sex", 'Age'
    ncol = 10,
    raster=TRUE,
    plot_max = 'q99',
    plot_min = 'q0',
    colfunc = rainbow,
    rev_colors=TRUE,
    dpi=400
  )
  
  pdf(paste0(fig_dir, 'human_spatial_featureplots/', cur_gene, '_featureplot.pdf'), width=20, height=8)
  print(p)
  dev.off()
  
  p <- VlnPlot(seurat_human, features = cur_gene, group.by = 'annotation', split.by = 'Diagnosis', pt.size=0) + xlab('')
  
  pdf(paste0(fig_dir, 'human_spatial_featureplots/', cur_gene, '_vlnplot.pdf'), width=7, height=4)
  print(p)
  dev.off()
  
}



plot_df <- seurat_human@meta.data
#plot_df <- plot_df %>% subset(Sample %in% selected_samples)
plot_df <- plot_df %>% subset(Sample == "Dec_13_2021_Human3")
vertices <- BayesSpace:::.make_hex_spots(plot_df, 'annotation')

splot <- ggplot(
  data = vertices,
  aes_(x = ~x.vertex, y = ~y.vertex,group = ~spot, fill = ~fill)
) +
  geom_polygon(size=0) +
  labs(fill = cur_gene) +
  #  scale_fill_gradientn(colors=viridis(256)) +
  coord_equal() +
  theme_void()


pdf(paste0(fig_dir, 'rep_human_feat_hex.pdf'), width=6, height=6)
print(splot)
dev.off()



