# 筛选神经元
Idents(seu) <- seu$celltype
neurons <- subset(seu, idents  = "Neuron")
table(neurons$celltype)


# 设置分组变量
neurons$group <- seurat_obj$group_variable # 确保这个是你分组的信息

# 查看PS相关基因的表达
PS_genes <- c("Anxa5", "BAX", "Casp3", "Casp8", "Tim4")
FeaturePlot(neurons, features = PS_genes, split.by = "group")
gene <- "Anxa5"
# 筛选掉该基因表达值为0的细胞
neurons_filtered <- subset(neurons, subset = eval(parse(text = gene_of_interest)) > 0)
neurons_filtered <- neurons[, neurons@assays$RNA$counts[gene,] > 0]
# 绘制小提琴图，隐藏单细胞点
VlnPlot(neurons_filtered, features = gene, group.by = "group", pt.size = 0)

VlnPlot(neurons, features = PS_genes, group.by = "group")

# 加载必要的R包
library(Seurat)
library(dplyr)
library(ggplot2)

# 定义感兴趣的基因
PS_genes <- c("Anxa5", "Casp3", "Casp8")
neurons <- JoinLayers(neurons)
# 循环遍历每个基因，并绘制小提琴图
for (gene in PS_genes) {
  # 筛选掉表达值为0的细胞
  neurons_filtered <- subset(neurons, subset = eval(parse(text = gene)) > 0)
  
  # 绘制小提琴图，隐藏单细胞点
  p <- VlnPlot(neurons_filtered, features = gene, group.by = "group", pt.size = 0) +
    ggtitle(paste("Expression of", gene)) +
    theme_minimal()
  
  print(p)
}



# test --------------------------------------------------------------------

# 确保基因名的大小写和Seurat对象一致
# 查看Seurat对象中的基因名
head(rownames(neurons))

# 如果需要转换基因名，可以这样做
PS_genes <- c("Anxa5", "Bax", "Casp3", "Casp8", "Tim4") # 请根据Seurat对象中的实际基因名调整

# 定义一个函数来绘制小提琴图，去除表达值为0的细胞
plot_violin_filtered <- function(seurat_obj, gene) {
  # 检查基因是否在Seurat对象中
  if (gene %in% rownames(seurat_obj)) {
    # 筛选掉表达值为0的细胞
    seurat_filtered <- subset(seurat_obj, subset = eval(parse(text = gene)) > 0)
    
    # 绘制小提琴图
    p <- VlnPlot(seurat_filtered, features = gene, group.by = "group", pt.size = 0) +
      ggtitle(paste("Expression of", gene)) +
      theme_minimal()
    
    print(p)
  } else {
    message(paste("Gene", gene, "not found in the Seurat object."))
  }
}

# 循环遍历每个基因，并绘制小提琴图
for (gene in PS_genes) {
  plot_violin_filtered(neurons, gene)
}




# pvalue ------------------------------------------------------------------

# 加载必要的R包
library(Seurat)
library(dplyr)

# 定义感兴趣的基因
PS_genes <- c("Anxa5", "Casp3", "Casp8")

# 获取group的前两个组别
groups <- unique(neurons_filtered$group)
group1 <- groups[1]
group2 <- groups[2]

# 创建一个空的数据框来存储结果
results <- data.frame(Gene = character(),
                      Group1 = character(),
                      Group2 = character(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# 对每个基因进行组间差异检验
for (gene in PS_genes) {
  if (gene %in% rownames(neurons_filtered@assays$RNA$counts)) {
    # 提取每组中该基因的表达值
    expression_group1 <- neurons_filtered@assays$RNA$counts[gene, neurons_filtered$group == group1]
    expression_group2 <- neurons_filtered@assays$RNA$counts[gene, neurons_filtered$group == group2]
    
    # 进行Wilcoxon秩和检验（非参数检验）
    test_result <- wilcox.test(expression_group1, expression_group2)
    
    # 保存结果
    results <- rbind(results, data.frame(Gene = gene,
                                         Group1 = group1,
                                         Group2 = group2,
                                         p_value = test_result$p.value))
  }
}

# 显示结果
results

