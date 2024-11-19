rm(list = ls())
getwd()
seu <- readRDS("rds/seu_annoted.rds")

library(Seurat)
Idents(seu)
table(seu$celltype)

# 提取小胶质细胞亚群
microglia_cells <- subset(seu, subset = celltype == "Microglia")  # GPT 的高级写法，无需修改Idents

# 标准预处理流程
microglia_cells <- NormalizeData(microglia_cells)
microglia_cells <- FindVariableFeatures(microglia_cells)
microglia_cells <- ScaleData(microglia_cells)
microglia_cells <- RunPCA(microglia_cells)

# 使用UMAP降维并进行聚类，调整resolution
microglia_cells <- RunUMAP(microglia_cells, dims = 1:20)
microglia_cells <- FindNeighbors(microglia_cells, dims = 1:20)
microglia_cells <- FindClusters(microglia_cells, resolution = 0.4)  # 分辨率可以根据需求调整

table(microglia_cells$seurat_clusters)
# 使用ROGUE评估分辨率
library(Seurat)
library(tidyverse)
# devtools::install_github("PaulingLiu/ROGUE")
library(ROGUE)

# 提取表达矩阵
expr <- Seurat::GetAssayData(microglia_cells) %>% as.matrix()

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

# 让第四群红色
micro_umap_pal <- ggsci::pal_npg()(6)
micro_umap_pal <- c("#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF","#E64B35FF","#8491B4FF")
library(ggplot2)
# UMAP图展示小胶质细胞亚群
micro_umap <- DimPlot(microglia_cells, reduction = "umap", split.by = "group",
                      pt.size = 0.6, alpha = 0.5, cols = micro_umap_pal,
                      label = TRUE, label.size = 6) + 
  # ggsci::scale_color_npg() +
  ggtitle("Microglia Subclusters"); micro_umap

ggsave("visual/micro_umap.pdf", micro_umap, width = 9, height = 6)


library(ggrastr); library(Nebulosa)
ggrastr::rasterize(Nebulosa::plot_density(microglia_cells,
                       c("Mertk","Pparg"),
                       size = 0.2), dpi = 300)
FeaturePlot(microglia_cells, c("Mertk","Pparg"))

# 按样本 `orig.ident` 分组的 UMAP 图
DimPlot(microglia_cells, reduction = "umap", ncol = 3,
        split.by = "group", label = T) + ggtitle("UMAP of Microglia Subclusters by Group")
saveRDS(microglia_cells, "rds/new_micro.rds")

# 找出每个小胶质细胞亚群的特异性marker
micro_all_markers <- FindAllMarkers(microglia_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 保存marker到CSV文件
write.csv(micro_all_markers, file = "rds/microglia_markers.csv")


# 亚群特异性marker展示 -----------------------------------------------------------
# 除了已经用过的 **DotPlot** 和 **热图 (Heatmap)** 外，还可以尝试以下可视化方式展示小胶质细胞 6 个亚群的 Marker 表达情况：
library(Seurat)
library(ggridges)
library(ggplot2)
library(dplyr)

# 筛选每个 cluster 中的 top 3 marker
top_markers <- micro_all_markers %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# 获取 top 3 基因列表
top_genes <- top_markers$gene

# 提取表达矩阵并转换为长格式
expr_data <- FetchData(microglia_cells, vars = c(top_genes, "seurat_clusters")) %>%
  pivot_longer(cols = top_genes, names_to = "gene", values_to = "expression") %>%
  mutate(seurat_clusters = factor(seurat_clusters))

# 绘制山脊图
ridge_plot <- ggplot(expr_data, aes(x = expression, y = gene, fill = seurat_clusters)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
  facet_wrap(~ seurat_clusters, scales = "free_x", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  labs(x = "Expression Level", y = "Gene Name", fill = "Cluster")

# 保存结果
ggsave(ridge_plot, filename = "visual/ridge_plot_fixed.pdf", width = 12, height = 10)


### 2. **Feature Plot (特征图)**

feature_plot <- FeaturePlot(microglia_cells, features = top_genes,
                            cols = c("lightgrey", "blue"), reduction = "umap")
ggsave(feature_plot, filename = "feature_plot.pdf", width = 12, height = 10)


### 3. **Violin Plot (小提琴图)**
# 合并小提琴图
violin_plot <- VlnPlot(
  microglia_cells,
  features = top_genes,
  group.by = "seurat_clusters",
  stack = TRUE,  # 关键参数，确保所有基因的图合并
  pt.size = 0      # 如果点太多，可以设置为 0 以隐藏
)+
  # scale_x_continuous(limits = c(0.01, NA)) +  # 设置 Y 轴起点略高于 0
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

ggsave(violin_plot, filename = "visual/violin_plot.pdf", width = 10, height = 8)


### 4. **Bubble Chart (气泡图)**
# - 用于展示基因在不同亚群中的表达量和检测比例，类似 dotplot，但可以增强视觉效果。
# - 示例代码：
# ```R
library(ggplot2)
library(dplyr)

# 确保基因名按 cluster 顺序排列
bubble_data <- top_markers %>%
  group_by(cluster, gene) %>%
  summarize(avg_exp = mean(avg_log2FC), pct_exp = mean(pct.1 * 100)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = unique(gene[order(cluster)])),  # 按 cluster 排列基因
         cluster = factor(cluster, levels = sort(unique(cluster))))  # 按 cluster 顺序排序

# 改进配色和样式
bubble_chart <- ggplot(bubble_data, aes(x = cluster, y = gene, size = pct_exp, color = avg_exp)) +
  geom_point(shape = 20) +  # 修改透明度和边框
  viridis::scale_color_viridis(option = "C", direction = -1)+
  # scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(bubble_data$avg_exp)) +  # 改为蓝白红渐变
  scale_size(range = c(3, 12)) +  # 调整气泡大小范围
  theme_minimal() +  # 极简主题
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),  # 横轴标签
    axis.text.y = element_text(size = 12, face = "bold"),  # 纵轴标签
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题
    panel.grid = element_blank(),   
    # panel.grid.major = element_line(color = "grey90", linetype = "dotted"),  # 柔和网格线
    legend.position = "right"  # 图例位置
  ) +
  labs(
    x = "Cluster",
    y = "Marker Gene",
    size = "Percent Expressed",
    color = "Average Expression"
  ) +
  guides(
    size = guide_legend(override.aes = list(color = "black", alpha = 1))  # 气泡大小图例保持一致颜色
  )

# 保存气泡图
ggsave(bubble_chart, filename = "bubble_chart_improved_ordered.pdf", width = 7, height = 8)



### 5. **Stacked Bar Plot (堆叠条形图)**
# - 展示不同亚群 marker 基因的总体分布情况。
# - 示例代码：
# ```R
marker_counts <- top_markers %>%
  count(cluster, gene) %>%
  mutate(prop = n / sum(n))

bar_plot <- ggplot(marker_counts, aes(x = cluster, y = prop, fill = gene)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(x = "Cluster", y = "Proportion", fill = "Gene")
ggsave(bar_plot, filename = "visual/stacked_bar_plot.pdf", width = 10, height = 8)


# 亚群分化方向 ------------------------------------------------------------------
# 加载所需包
library(CytoTRACE2)
library(Seurat)

# 提取表达矩阵
expression_data <- as.matrix(GetAssayData(microglia_cells, slot = "counts"))

# 确保行名和列名正确
if (is.null(rownames(expression_data)) | is.null(colnames(expression_data))) {
  stop("表达矩阵必须具有行名（基因名）和列名（细胞名）")
}

# 生成注释数据
annotation <- data.frame(
  celltype = microglia_cells$seurat_clusters,  # 假设你的 Seurat 对象中有 celltype 列
  group = microglia_cells$group        # 假设你的 Seurat 对象中有 group 列
)
rownames(annotation) <- colnames(microglia_cells)

# 运行 CytoTRACE 2
cytotrace2_result <- cytotrace2(microglia_cells, is_seurat = T)

# 可视化 CytoTRACE 2 结果
plots <- plotData(
  cytotrace2_result = cytotrace2_result, 
  annotation = annotation, is_seurat = T
)

# 查看生成的 UMAP 图
print(plots$Phenotype_UMAP)

# 获取 Seurat 的 UMAP 坐标
umap_coordinates <- Embeddings(microglia_cells, reduction = "umap")

# 将坐标传递给 CytoTRACE 可视化
plots <- plotData(
  cytotrace2_result = cytotrace2_result,
  annotation = annotation,
  umap = umap_coordinates  # 使用 Seurat 的 UMAP 坐标
)

# 将 CytoTRACE 结果添加到 Seurat 对象
microglia_cells$CytoTRACE2_Potency <- cytotrace2_result$CytoTRACE2_Potency
microglia_cells$CytoTRACE2_Relative <- cytotrace2_result$CytoTRACE2_Relative

# 在 Seurat 的 UMAP 上展示 CytoTRACE2 分数
FeaturePlot(
  object = microglia_cells, 
  features = c("CytoTRACE2_Relative"), 
  reduction = "umap"
)

micro_cytotrave <- FeaturePlot(microglia_cells, "CytoTRACE2_Relative",pt.size = 1.5) + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)
micro_cytotrave
# 如果需要保存图像：
ggsave(micro_cytotrave, filename = "visual/micro_cytotrave.pdf", width = 10, height = 8)

saveRDS(cytotrace2_result,"rds/micro_cytotrace.rds")

# Mertk Density -----------------------------------------------------------


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




