##作图
library(ggplot2)
dir.create("visual/fig1")
library(RColorBrewer)
col_set <- brewer.pal(n=length(unique(scRNA$celltype)),name = "Set3")
# Fig 1a1
umap.singler <- DimPlot(scRNA, reduction = "umap", group.by='celltype',
                        repel=T, label=T, label.size=5,
                        cols = col_set)+labs(title = "")
ggsave("visual/fig1/singleR.pdf", umap.singler, width = 16, height = 9)
# Fig 1a2
umap.orig <- DimPlot(scRNA, reduction = "umap", group.by='orig.ident',
                     repel=F, label=F, label.size=5)+labs(title = "")
ggsave("visual/fig1/orig.pdf", umap.orig, width = 16, height = 9)
#Fig 1a3
umap.seur <- DimPlot(scRNA, reduction = "umap", group.by='seurat_clusters',
                     repel=T, label=T, label.size=5)+labs(title = "")
ggsave("visual/fig1/seur.pdf", umap.seur, width = 16, height = 9)
#Fig 1a4
umap.singler.split <- DimPlot(scRNA, reduction = "umap", group.by='celltype',
                              split.by = "orig.ident", ncol = 2, 
                              repel=T, label=T, label.size=5)+labs(title = "")
ggsave("visual/fig1/singler-split.pdf", umap.singler.split, width = 16, height = 9)


## Fig1b proportion
## Cell Proportion
cellprop <- table(scRNA$orig.ident,scRNA$celltype)
cellprop <- proportions(cellprop,1)
cellpro <- as.data.frame(cellprop)
fig.prop <- list()
dir.create("visual/fig1")
for (sample.orig in unique(cellpro$Var1)) {
  sample_prop <- cellpro[cellpro$Var1==sample.orig,]
  sample_prop$ymax=cumsum(sample_prop$Freq)
  sample_prop$ymin=c(0,head(sample_prop$ymax,n=-1))
  #设置标签
  lab=paste0(sample_prop$Var2,'\n',round(sample_prop$Freq*100,1),"%")
  sample_prop$lab=lab
  sample_prop
  
  fig.prop[[sample.orig]] <- ggplot(sample_prop,aes(ymax=ymax,ymin=ymin,
                                                    xmax=4,xmin=2))+
    geom_rect(aes(fill=Var2))+
    theme_void()+
    xlim(0,4)+
    coord_polar(theta="y")+
    scale_fill_manual(values = col_set) +
    geom_text(aes(y= (ymax + ymin)/2,x=3,label=lab),size=5,
              check_overlap = F)+
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, vjust = -6))+
    ggtitle(sample.orig)
  ggsave(paste0("visual/fig1/",sample.orig,".pdf"),
         fig.prop[[sample.orig]],
         width = 16, height = 9)
  
}

prop.plot <- fig.prop[[1]]

for (i in 2:length(fig.prop)) {
  prop.plot <- prop.plot+fig.prop[[i]]
  
}

ggsave("visual/fig1/patch.pdf", prop.plot, width = 16, height = 9)

