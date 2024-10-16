library(Seurat)
library(ggplot2)
seu <- readRDS("~/Downloads/scAD/rds/done_annotated.rds")
head(seu)
DimPlot(seu, label = T, group.by = "celltype")

micro_sub <- readRDS("~/Downloads/scAD/rds/micro_sub.rds")

# 确保Idents已经设置为seurat_clusters
Idents(micro_sub) <- micro_sub$seurat_clusters

# 为microglia子集添加subtype列，按cluster编号进行命名
micro_sub$subtype <- paste0("Mg", as.numeric(Idents(micro_sub)))

# 为neuron子集添加subtype列
neuron_sub <- subset(seu, celltype == "Neuron")
neuron_sub$subtype <- "Neuron"

# 检查是否添加成功
table(micro_sub$subtype)
table(neuron_sub$subtype)

# 使用merge函数合并两个对象
seu_merged <- merge(micro_sub, neuron_sub)

# 检查合并后的meta data是否包含subtype列
head(seu_merged@meta.data)
table(seu_merged$subtype)

# 对合并后的对象进行标准化处理
seu_merged <- NormalizeData(seu_merged)
seu_merged <- FindVariableFeatures(seu_merged)
seu_merged <- ScaleData(seu_merged)

# 如果有必要，可以再次进行PCA和UMAP降维
seu_merged <- RunPCA(seu_merged)
seu_merged <- RunUMAP(seu_merged, dims = 1:20)

# 查看聚类效果，确保正确处理
DimPlot(seu_merged, label = TRUE, group.by = "subtype")

# 创建CellChat对象
library(CellChat)
seu_merged <- JoinLayers(seu_merged)
# 设置meta数据，group.by参数使用subtype
cellchat <- createCellChat(object = seu_merged, meta = seu_merged@meta.data, group.by = 'subtype')

# 使用小鼠的CellChat数据库
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

# 识别过表达基因和细胞间互作
cellchat <- subsetData(cellchat)  # 提取必要的数据
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)  # 投影到蛋白质相互作用网络

# 计算通讯概率
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

# 聚合网络并分析通讯
cellchat <- aggregateNet(cellchat)
print(cellchat);head(cellchat)
# 可视化通讯
netVisual_circle(cellchat, weight.scale = TRUE)

# 对比不同组的通讯
compareInteractions(cellchat, show.legend = FALSE)

# 根据group拆分数据
sham <- subset(seu_merged, group == "Control")
mcao <- subset(seu_merged, group == "AD + WT")
ko <- subset(seu_merged, group == "AD + Mertk-/-")

# 创建CellChat对象，专注于中性粒细胞和内皮细胞的互作
library(CellChat)
cellchat_sham <- createCellChat(object = sham, meta = sham@meta.data, group.by = 'subtype')
cellchat_mcao <- createCellChat(object = mcao, meta = mcao@meta.data, group.by = 'subtype', assay = 'RNA')
cellchat_ko <- createCellChat(object = ko, meta = ko@meta.data, group.by = 'subtype', assay = 'RNA')

# 使用小鼠数据库
CellChatDB <- CellChatDB.mouse
head(CellChatDB$interaction)
cellchat_sham@DB <- CellChatDB
cellchat_mcao@DB <- CellChatDB
cellchat_ko@DB <- CellChatDB
# 重复之前的分析步骤，分别对每个组进行通讯分析
cellchat_sham <- subsetData(cellchat_sham)
cellchat_mcao <- subsetData(cellchat_mcao)
cellchat_ko <- subsetData(cellchat_ko)

# 按组进行分析
cellchat_sham <- identifyOverExpressedGenes(cellchat_sham)
cellchat_mcao <- identifyOverExpressedGenes(cellchat_mcao)
cellchat_ko <- identifyOverExpressedGenes(cellchat_ko)

cellchat_sham <- computeCommunProb(cellchat_sham)
cellchat_mcao <- computeCommunProb(cellchat_mcao)
cellchat_ko <- computeCommunProb(cellchat_ko)

# 保存每个组的分析结果
save(cellchat_sham, cellchat_mcao, cellchat_ko, file = "rds/cellchat_results_micro_neuron.rda")


