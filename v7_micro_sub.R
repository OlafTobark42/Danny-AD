# microlia subset
library(Seurat)
microglia_cells <- subset(seu, celltype %in% "Microglia")
microglia_cells <- readRDS("rds/micro_sub.rds")
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

seu <- microglia_cells
saveRDS(microglia_cells, file = "rds/micro_sub.rds")
library(ggplot2)
# UMAP图展示小胶质细胞亚群
DimPlot(microglia_cells, reduction = "umap", label = TRUE) 
+ ggtitle("Microglia Subclusters")

library(ggrastr); library(Nebulosa)
ggrastr::rasterize(Nebulosa::plot_density(microglia_cells,
                                          c("Mertk","Pparg"),
                                          size = 0.2), dpi = 300)
FeaturePlot(microglia_cells, c("Mertk"),max.cutoff = 5,split.by = "group")


dam_sub <- subset(microglia_cells, seurat_clusters %in% 3)
table(dam_sub$seurat_clusters)
VlnPlot(dam_sub,"Mertk",group.by = "group")
VlnPlot(dam_sub,"Cxcl9",group.by = "group")
Idents(dam_sub) <- dam_sub$group
table(dam_sub$group)
dam_deg <- FindMarkers(dam_sub, ident.1 = "AD + WT", ident.2 = "Control")
library(tidyverse)
dam_top_deg <- top_n(dam_deg,n = 10)

# 按样本 `orig.ident` 分组的 UMAP 图
DimPlot(microglia_cells, reduction = "umap", ncol = 3,
        split.by = "group", label = T) + ggtitle("UMAP of Microglia Subclusters by Group")


# 找出每个小胶质细胞亚群的特异性marker
microglia_markers <- FindAllMarkers(microglia_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Load necessary libraries
library(Seurat)
library(pheatmap)
library(stringr)
# List of marker genes (example: from Homeostatic, DAM, LDAM categories)
marker_genes <- str_to_title(c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", 
                  "CD9", "ITGAX", "CLEC7A", "CD63", "SPP1", "LPL", "TREM2", "APOE",
                  "NAMPT", "ACSL1", "DPYD", "CD163"))

# Check if the genes exist in the dataset
available_genes <- rownames(seu[["RNA"]]@layers)

# Check which of your marker genes are not in the dataset
setdiff(marker_genes, available_genes)


# Create heatmap using Seurat's DoHeatmap function
heatmap_plot <- DoHeatmap(seu, features = marker_genes, size = 4, group.by = "seurat_clusters") +
  scale_fill_gradientn(colors = c("navy", "blue", "green", "yellow", "red"))

# Prepare expression matrix for heatmap (markers by cluster)
data_matrix <- FetchData(seu, vars = marker_genes)
head(data_matrix)
# Create heatmap using pheatmap with fixed cluster widths and viridis palette
pheatmap::pheatmap(data_matrix_transposed, show_rownames = T, show_colnames = F,
                   cluster_cols = F,  # Clustering columns
                   cluster_rows = F,  # Clustering rows
                   color = viridis::viridis(50),  # Viridis color palette
                   cellwidth = 15,       # Fixed cell width
                   cellheight = 15,      # Fixed cell height
                   scale = "row")        # Z-score scaling across rows (genes)

# Assuming your data is called data_matrix
# Transpose the matrix to have genes as rows and cells as columns
data_matrix_transposed <- t(data_matrix) %>% as.data.frame()

# Now use pheatmap
library(pheatmap)

# Plot the heatmap with genes as rows
pheatmap(data_matrix_transposed)


# 加载必要的库
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠基因注释
library(enrichplot)
library(ggplot2)

# 定义上下调基因列表
# 上调基因列表
upregulated_genes <- c('Cxcl9', 'Gpnmb', 'Apold1', 'Gpr157', 'Gtf3c1', 'Thbs1', 
                       'Apoc4', 'Cdk5rap2', 'Ganc', 'Dcstamp', 'Sufu', 'Ppp6r2', 
                       'Apoo', 'Egln3', '4930557K07Rik', 'Zfp512', 'Cd99l2', 
                       'Ms4a7', 'Pianp', 'Nusap1', 'Sirt1', 'Pygo2', 'Iqgap2', 
                       'Spsb3', 'Nr6a1os', 'Hook2', 'B3glct', 'Pik3r2', 'Mrs2', 
                       'Bcl2l12', 'Trp53inp1', 'Clybl', 'Srp54a', 'Ube2c', 
                       'Msr1', 'Trappc13', '6030458C11Rik', 'Mfn1', 'Metrn', 
                       'Mdc1', 'Fastkd2', 'Eif4enif1', 'Hibch', 'Pak4', 
                       'Rad51c', 'Sin3a', 'Nedd9', 'Cdc42ep4', 'Igf1', 'Casp9')

downregulated_genes <- c('Draxin', '4930447C04Rik', 'Zyg11a', 'Lamc3', 'Gm10655', 
                         'Insc', 'Ptcra', 'Gm46367', 'Dzank1', 'Bcas1os1', 
                         'Qrich2', 'Rbp4', 'Tvp23bos', 'Gm50333', 'Pbld1', 
                         '0610039K10Rik', 'Mreg', 'Rhox5', 'Gm19325', 'Lins1')

# 转换基因符号为ENTREZ ID
upregulated_entrez <- bitr(upregulated_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
downregulated_entrez <- bitr(downregulated_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# GO富集分析
up_go <- enrichGO(gene = upregulated_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
down_go <- enrichGO(gene = downregulated_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
up_go
# KEGG富集分析
up_kegg <- enrichKEGG(gene = upregulated_entrez$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)
down_kegg <- enrichKEGG(gene = downregulated_entrez$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

# 绘制GO富集结果
dotplot(up_go, showCategory=10) + ggtitle("Upregulated Genes - GO Enrichment")
dotplot(down_go, showCategory=10) + ggtitle("Downregulated Genes - GO Enrichment")

# 绘制KEGG富集结果
dotplot(up_kegg, showCategory=10) + ggtitle("Upregulated Genes - KEGG Enrichment")
dotplot(down_kegg, showCategory=10) + ggtitle("Downregulated Genes - KEGG Enrichment")

Idents(microglia_cells) <- microglia_cells$seurat_clusters
m3_marker <- FindMarkers(microglia_cells, ident.1 = 3)

# 加载必要的R包
library(clusterProfiler)
library(org.Mm.eg.db)  # 小鼠基因注释
library(enrichplot)
library(ggplot2)

# 假设 m3_marker 数据框已经存在
# head(m3_marker) 的输出应该是类似这样的：
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# Lgals3bp     0   3.071960 0.814 0.200         0
# Cst7         0   4.254675 0.632 0.038         0
# Ifi27l2a     0   3.888764 0.629 0.063         0
# Bcl2a1d      0   2.526854 0.747 0.217         0
# Gpr65        0   2.790305 0.627 0.099         0
# Lyz2         0   2.926828 0.900 0.378         0

# 筛选最上调和最下调的100个基因
top_upregulated <- m3_marker[order(m3_marker$avg_log2FC, decreasing = TRUE), ][1:100, ]
top_downregulated <- m3_marker[order(m3_marker$avg_log2FC, decreasing = FALSE), ][1:100, ]

# 提取基因名称列
upregulated_genes <- rownames(top_upregulated)
downregulated_genes <- rownames(top_downregulated)
upregulated_genes <- c(upregulated_genes,downregulated_genes)
# 转换基因符号为ENTREZ ID
upregulated_entrez <- bitr(upregulated_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
downregulated_entrez <- bitr(downregulated_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# GO富集分析 (上调的基因)
up_go <- enrichGO(gene = upregulated_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
up_go
# GO富集分析 (下调的基因)
down_go <- enrichGO(gene = downregulated_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# KEGG富集分析 (上调的基因)
up_kegg <- enrichKEGG(gene = upregulated_entrez$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

# KEGG富集分析 (下调的基因)
down_kegg <- enrichKEGG(gene = downregulated_entrez$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05)

# 绘制GO富集结果 (上调)
dotplot(up_go, showCategory=10) + ggtitle("Upregulated Genes - GO Enrichment")

# 绘制GO富集结果 (下调)
dotplot(down_go, showCategory=10) + ggtitle("Downregulated Genes - GO Enrichment")

# 绘制KEGG富集结果 (上调)
dotplot(up_kegg, showCategory=10) + ggtitle("Upregulated Genes - KEGG Enrichment")

# 绘制KEGG富集结果 (下调)
dotplot(down_kegg, showCategory=10) + ggtitle("Downregulated Genes - KEGG Enrichment")


# 加载必要的R包
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSVA)
library(GSEABase)
library(ggplot2)
table(Idents(dam_sub))
m3_ad_ko <- FindMarkers(dam_sub, ident.1 = "AD + Mertk-/-",
                        ident.2 = "AD + WT")
m3_marker <- m3_ad_ko
# 提取 foldchange 最大的 250 个基因
top_upregulated <- m3_marker %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 1000) %>%
  rownames_to_column(var = "gene")

# 提取 foldchange 最小的 250 个基因
top_downregulated <- m3_marker %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 1000) %>%
  rownames_to_column(var = "gene")

# 合并上调和下调的基因
selected_genes <- bind_rows(top_upregulated, top_downregulated)

# 创建gene_list，基于 avg_log2FC 值
gene_list <- selected_genes$avg_log2FC
names(gene_list) <- selected_genes$gene

# 按降序排列基因列表
gene_list <- sort(gene_list, decreasing = TRUE)
# 假设你已经有一个 m3_marker 数据框
# 基因表达数据应该以 avg_log2FC 排列
# 先将基因名称和foldchange值保存到向量中
gene_list <- m3_marker$avg_log2FC
names(gene_list) <- rownames(m3_marker)

# 按降序排列基因列表
gene_list <- sort(gene_list, decreasing = TRUE)
library(tidyverse)
# gene_list <- top_n(gene_list, n = 500)
# GSEA分析准备：导入吞噬功能相关的基因集
# 你可以自定义一个基因集，或者从 MSigDB 下载包含吞噬功能的基因集（例如GO的phagocytosis相关基因集）

# 示例：手动定义一个小胶质细胞吞噬相关基因集
phagocytosis_genes <- c('Apoe', 'Cd68', 'Itgam', 'Trem2', 'Lyz2', 'Csf1r', 'Fcgr1', 'Fcgr2b')

microglia_cells <- AddModuleScore(microglia_cells,
                                  features = phagocytosis_genes,
                                  name = "phago")

head(microglia_cells)
# Home1
mydata <- FetchData(microglia_cells, vars = c("umap_1", "umap_2", "phago1"))

# 计算颜色梯度的上下限
Home1_min <- quantile(mydata$phago1, 0.01)  # 设置5%作为最低点
Home1_max <- quantile(mydata$phago1, 0.95)  # 设置95%作为最高点
library(viridis)
# 调整颜色梯度的上下限并添加梯度线阈值
ggplot(mydata, aes(x = umap_1, y = umap_2, colour = phago1)) +
  geom_point(size = 1) +  # 调整点的大小
  scale_color_viridis(option = "viridis") +  # 限制颜色梯度会导致灰色点太多 , limits = c(Home1_min, Home1_max)
  theme_minimal() +  # 使用简约主题
  theme(
    legend.position = "none",  # 完全移除 legend
    panel.background = element_rect(fill = "black"),  # 设置黑色背景
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),  # 移除小网格线
    plot.background = element_rect(fill = "black"),  # 设置黑色绘图区背景
    axis.title = element_text(color = "white"),  # 坐标轴标题颜色
    axis.text = element_blank()  # 移除坐标轴刻度
  ) +
  geom_density_2d(data = subset(mydata, phago1 > Home1_max * 0.6), color = "white")   # 添加白色梯度线（仅在高亮区域）



# 准备一个GSEA geneSet对象
phagocytosis_set <- list(phagocytosis = phagocytosis_genes)
gsea_set <- GSEABase::GeneSet(phagocytosis_genes, setName = "Phagocytosis")

# 运行GSEA分析
gsea_results <- GSEA(gene_list, TERM2GENE = data.frame(term="Phagocytosis", gene=phagocytosis_genes), 
                     pvalueCutoff = 0.05)

# 打印GSEA结果
print(gsea_results)

# 绘制GSEA结果图
gseaplot2(gsea_results, geneSetID = "Phagocytosis", title = "GSEA - Phagocytosis")

# GSVA分析
# 如果你有不同样本的表达矩阵，可以进行GSVA
# 假设expr_matrix是表达矩阵，行是基因，列是样本
# rownames(expr_matrix) <- 基因名
# colnames(expr_matrix) <- 样本名

# 示例数据: 你可以根据自己的数据格式提供表达矩阵
# GSVA分析
gsva_results <- gsva(as.matrix(expr_matrix), phagocytosis_set, method = "gsva")

# 查看GSVA结果
head(gsva_results)

# 可视化GSVA结果
boxplot(gsva_results[1, ] ~ group, main="GSVA - Phagocytosis Activity", ylab="GSVA score")




