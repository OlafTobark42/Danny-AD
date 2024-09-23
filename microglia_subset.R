rm(list = ls())
getwd()
scRNA <- readRDS("rds/scRNA_annoted.rds")

library(Seurat)
Idents(scRNA)
table(scRNA$celltype)

# 如果idents是Seurat clusters，先将Idents修改为celltype
Idents(scRNA) <- "celltype"
Idents(scRNA)
# 提取小胶质细胞亚群
microglia_cells <- subset(scRNA, idents = "Microglia")

# 标准预处理流程
microglia_cells <- NormalizeData(microglia_cells)
microglia_cells <- FindVariableFeatures(microglia_cells)
microglia_cells <- ScaleData(microglia_cells)
microglia_cells <- RunPCA(microglia_cells)

# 使用UMAP降维并进行聚类，调整resolution
microglia_cells <- RunUMAP(microglia_cells, dims = 1:20)
microglia_cells <- FindNeighbors(microglia_cells, dims = 1:20)
microglia_cells <- FindClusters(microglia_cells, resolution = 0.2)  # 分辨率可以根据需求调整

table(microglia_cells$seurat_clusters)
# 使用ROGUE评估分辨率
library(Seurat)
library(tidyverse)
library(ROGUE)

# 提取表达矩阵
expr <- LayerData(microglia_cells, layer = "counts") %>% as.matrix()

# 过滤矩阵，去除在少于10个细胞中表达的基因
expr <- matr.filter(expr, min.cells = 10)

# 计算样本的 Shannon Entropy
ent.res <- SE_fun(expr)

# 绘制 Shannon Entropy 图
SEplot(ent.res)

# 计算 ROGUE 值
rogue.res <- rogue(expr, samples = "orig.ident",
                   labels = microglia_cells$seurat_clusters, platform = "UMI", span = 0.6)

# 查看 ROGUE 值
print(rogue.res)

# 绘制 ROGUE boxplot
rogue.boxplot(rogue.res)


# 可视化clustree
library(clustree)
clustree(microglia_cells)

library(ggplot2)
# UMAP图展示小胶质细胞亚群
DimPlot(microglia_cells, reduction = "umap", label = TRUE) 
  + ggtitle("Microglia Subclusters")

library(ggrastr); library(Nebulosa)
ggrastr::rasterize(Nebulosa::plot_density(microglia_cells,
                       c("Mertk","Pparg"),
                       size = 0.2), dpi = 300)
FeaturePlot(microglia_cells, c("Mertk","Pparg"))

# 按样本 `orig.ident` 分组的 UMAP 图
DimPlot(microglia_cells, reduction = "umap", ncol = 3,
        split.by = "group", label = T) + ggtitle("UMAP of Microglia Subclusters by Group")


# 找出每个小胶质细胞亚群的特异性marker
microglia_markers <- FindAllMarkers(microglia_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 保存marker到CSV文件
write.csv(microglia_markers, file = "microglia_markers.csv")

# 找WT vs AD是否有Mertk
Idents(microglia_cells) <- "group"
head(microglia_cells)

wt_vs_ad_deg <- FindMarkers(microglia_cells, 
                            ident.1 = "AD + WT",
                            ident.2 = "Control")


wt_vs_ad_deg
dir.create("deg_result")
write.csv(wt_vs_ad_deg, "deg_result/wt_vs_ad_deg.csv")

wt_vs_ad_deg["Ppar",]

# KO vs AD
ad_vs_ko <- FindMarkers(microglia_cells, 
                        ident.1 = "AD + Mertk-/-",
                        ident.2 = "AD + WT")
ad_vs_ko["Ppar",]

# 上下调富集
# 加载所需的R包
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠基因注释包
library(ggplot2)
library(enrichplot)  # 可视化包
library(DOSE)
# 未筛选P值
if (F) {# 读取差异表达数据
  deg_data <- read.csv("deg_result/wt_vs_ad_deg.csv")
  
  # 分别提取上调和下调的基因列表
  top100_upregulated <- deg_data %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    top_n(100, avg_log2FC)
  
  top100_downregulated <- deg_data %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    top_n(30, desc(avg_log2FC))
  # Error :这个错误的原因是你在 top_n() 函数中用了 desc(avg_log2FC)，但 top_n() 已经隐含了一个排序方向。
  # 在 top_n() 中使用 desc() 会引发冲突，因为 top_n() 本身默认会对数值型数据从小到大排序。因此，修复方法是去掉 desc()，直接传递列名即可。
  top100_downregulated <- deg_data %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    top_n(100, avg_log2FC)  # 注意这里不需要用desc
  
  
  # 提取基因列表 (注意你的基因列的名字可能需要调整)
  upregulated_genes <- top100_upregulated$Unnamed..0
  downregulated_genes <- top100_downregulated$Unnamed..0
}
# 先只筛选P<0.05而非adj P
if (T) {
  # 加载所需的R包
  library(dplyr)
  
  # 读取差异表达数据
  deg_data <- read.csv("deg_result/wt_vs_ad_deg.csv")
  
  # 首先筛选 p 值小于 0.05 的基因
  significant_deg <- deg_data %>%
    filter(p_val < 0.05)
  
  # 分别提取上调和下调的基因列表，先按照普通p值筛选，再进行log2FC排序
  top100_upregulated <- significant_deg %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    top_n(100, avg_log2FC)
  
  top100_downregulated <- significant_deg %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    top_n(100, avg_log2FC)
  
  # 提取基因列表 (注意基因列的名字可能需要调整)
  upregulated_genes <- top100_upregulated$X
  downregulated_genes <- top100_downregulated$X
}

## Step 2: 进行GO富集分析（上下调基因分别分析）
# 上调基因的GO富集分析
up_go_enrich <- enrichGO(gene          = upregulated_genes,
                         OrgDb         = org.Mm.eg.db,
                         keyType       = "SYMBOL",
                         ont           = "ALL",  # 使用BP, MF, CC
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

# 下调基因的GO富集分析
down_go_enrich <- enrichGO(gene          = downregulated_genes,
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)

## Step 3: 富集结果可视化（区分上下调）
# # 设置上调和下调基因的富集分析结果输出格式
up_go_enrich_df <- as.data.frame(up_go_enrich)
down_go_enrich_df <- as.data.frame(down_go_enrich)

# 绘制上调基因的富集气泡图
dotplot(up_go_enrich, split = "ONTOLOGY", showCategory = 10) + 
  facet_grid(ONTOLOGY~., scale="free") +
  ggtitle("Upregulated Gene GO Enrichment") +
  theme_minimal(base_size = 15)

# 绘制下调基因的富集气泡图
dotplot(down_go_enrich, split = "ONTOLOGY", showCategory = 10) + 
  facet_grid(ONTOLOGY~., scale="free") +
  ggtitle("Downregulated Gene GO Enrichment") +
  theme_minimal(base_size = 15)


## Step 4: 进行KEGG通路富集分析
# Step 1: 使用 bitr() 将基因符号转换为 Entrez ID
# 加载所需的R包
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠基因注释包

# 将上调基因转换为Entrez ID
upregulated_genes_entrez <- bitr(upregulated_genes, 
                                 fromType = "SYMBOL", 
                                 toType = "ENTREZID", 
                                 OrgDb = org.Mm.eg.db)

# 将下调基因转换为Entrez ID
downregulated_genes_entrez <- bitr(downregulated_genes, 
                                   fromType = "SYMBOL", 
                                   toType = "ENTREZID", 
                                   OrgDb = org.Mm.eg.db)

# 上调基因的KEGG富集分析
up_kegg_enrich <- enrichKEGG(gene         = upregulated_genes_entrez$ENTREZID,
                             organism     = 'mmu',  # 小鼠物种代码
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05)

# 下调基因的KEGG富集分析
down_kegg_enrich <- enrichKEGG(gene         = downregulated_genes_entrez$ENTREZID,
                               organism     = 'mmu',
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

# 检查结果并可视化
if (!is.null(up_kegg_enrich)) {
  dotplot(up_kegg_enrich, showCategory = 10) + 
    ggtitle("Upregulated Gene KEGG Enrichment") +
    theme_minimal(base_size = 15)
}

if (!is.null(down_kegg_enrich)) {
  dotplot(down_kegg_enrich, showCategory = 10) + 
    ggtitle("Downregulated Gene KEGG Enrichment") +
    theme_minimal(base_size = 15)
}

# 提取通路名字
# 如果 up_kegg_enrich 是有效的 enrichResult 对象，提取通路名称
if (!is.null(up_kegg_enrich)) {
  # 将 KEGG 富集结果转换为数据框
  up_kegg_enrich_df <- as.data.frame(up_kegg_enrich)
  
  # 提取通路名字
  kegg_pathways <- up_kegg_enrich_df$Description
  
  # 打印前几条通路名字
  print(kegg_pathways)
}

## Step 2: 读取数据并制作火山图
# 加载所需的包
library(ggplot2)
library(ggrepel)
library(dplyr)

# 读取差异表达数据
# deg_data <- read.csv("/path/to/your/file/wt_vs_ad_deg.csv")

# 添加一列 -log10(p值) 以便于绘制火山图
deg_data <- deg_data %>%
  mutate(logP = -log10(p_val))

# 设置基因颜色分类规则：上调为红色，下调为蓝色，其他为灰色
deg_data <- deg_data %>%
  mutate(color = ifelse(X %in% c("Mertk", "Pparg"), "highlight",
                        ifelse(avg_log2FC > 1 & p_val < 0.05, "upregulated",
                               ifelse(avg_log2FC < -1 & p_val < 0.05, "downregulated", "notsig"))))

# 绘制火山图
volcano_plot <- ggplot(deg_data, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = color), size = 2) +  # 设置点的大小为2
  scale_color_manual(values = c("highlight" = "red", "upregulated" = "red", "downregulated" = "blue", "notsig" = "grey")) +
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none") +  # 隐藏图例
  labs(title = "Volcano Plot: WT vs AD Microglia", x = "Log2 Fold Change", y = "-log10(p-value)") +
  
  # 添加 Pparg 的标签，直接高亮显示名字
  geom_text_repel(data = subset(deg_data, X == "Pparg"),
                  aes(label = X), size = 5, box.padding = 0.5, point.padding = 0.5, 
                  segment.color = 'black') +
  
  # 为 Mertk 添加线和名字
  geom_point(data = subset(deg_data, X == "Mertk"), aes(x = avg_log2FC, y = logP), color = "red", size = 3) +
  geom_text_repel(data = subset(deg_data, X == "Mertk"),
                  aes(label = X), size = 5, box.padding = 0.5, point.padding = 0.5,
                  segment.color = 'black', arrow = arrow(length = unit(0.02, "npc")))  # 使用箭头线连接

# 显示火山图
print(volcano_plot)




