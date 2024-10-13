rm(list = ls())
library(Seurat)
# 设置工作目录并创建cellchat文件夹
if (!dir.exists("cellchat")) dir.create("cellchat")
setwd("~/Downloads/scAD")
getwd()

seu <- readRDS("rds/manual_annoted.rds")
head(seu); colnames(seu@meta.data)
table(seu$celltype)
seu$celltype <- ifelse(is.na(seu$celltype),"Unknown",seu$celltype)
# Step 13: Cellchat -------------------------------------------------------
Idents(seu) <- seu$group
# 提取MCAO和Sham组
sham <- subset(seu, ident = "Control")
mcao <- subset(seu, ident = "AD + WT")
ko <- subset(seu, ident = "AD + Mertk-/-")

# 创建CellChat对象，专注于中性粒细胞和内皮细胞的互作
library(CellChat)
cellchat_sham <- createCellChat(object = sham, meta = sham@meta.data, group.by = 'celltype', assay = 'RNA')
cellchat_mcao <- createCellChat(object = mcao, meta = mcao@meta.data, group.by = 'celltype', assay = 'RNA')
cellchat_ko <- createCellChat(object = ko, meta = ko@meta.data, group.by = 'celltype', assay = 'RNA')

# 使用小鼠数据库
CellChatDB <- CellChatDB.mouse
head(CellChatDB$interaction)
cellchat_sham@DB <- CellChatDB
cellchat_mcao@DB <- CellChatDB
cellchat_ko@DB <- CellChatDB

library(future)
library(parallel)
detectCores()
plan(strategy = 'multisession', workers = 6)
options(future.globals.maxSize = 4 * 1024^3)  # 设置为2GB 否则超过500会报错“Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the 4 globals exported for future expression (‘FUN()’) is 523.14 MiB.. 
# This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize'). 
# The three largest globals are ‘FUN’ (499.75 MiB of class ‘function’), ‘data2’ (22.60 MiB of class ‘numeric’) and ‘data1’ (804.92 KiB of class ‘numeric’)”
# plan(strategy = 'multisession', workers = 48)  # 继续使用multisession并调整核心数


# 识别基因和互作
cellchat_sham <- subsetData(cellchat_sham)
Sys.time()
cellchat_sham <- identifyOverExpressedGenes(cellchat_sham)
print("identifyOverExpressedGenes_done");Sys.time()
cellchat_sham <- identifyOverExpressedInteractions(cellchat_sham)
cellchat_sham <- projectData(cellchat_sham, PPI.mouse)
Sys.time()
cellchat_sham <- computeCommunProb(cellchat_sham)
print("computeCommunProb_done");Sys.time()

cellchat_sham <- computeCommunProbPathway(cellchat_sham)
cellchat_sham <- aggregateNet(cellchat_sham)
levels(cellchat_sham@idents)            #查看细胞顺序
# vertex.receiver = c(3, 6)          #指定靶细胞的索引 内皮3 中性粒 8
cellchat_sham@netP$pathways       


cellchat_mcao <- subsetData(cellchat_mcao)
cellchat_mcao <- identifyOverExpressedGenes(cellchat_mcao)
cellchat_mcao <- identifyOverExpressedInteractions(cellchat_mcao)
cellchat_mcao <- projectData(cellchat_mcao, PPI.mouse)
cellchat_mcao <- computeCommunProb(cellchat_mcao)
cellchat_mcao <- computeCommunProbPathway(cellchat_mcao)
cellchat_mcao <- aggregateNet(cellchat_mcao)
levels(cellchat_mcao@idents)            #查看细胞顺序
# vertex.receiver = c(3, 6)          #指定靶细胞的索引 内皮3 中性粒 8
cellchat_mcao@netP$pathways   

cellchat_ko <- subsetData(cellchat_ko)
cellchat_ko <- identifyOverExpressedGenes(cellchat_ko)
cellchat_ko <- identifyOverExpressedInteractions(cellchat_ko)
cellchat_ko <- projectData(cellchat_ko, PPI.mouse)
cellchat_ko <- computeCommunProb(cellchat_ko)
cellchat_ko <- computeCommunProbPathway(cellchat_ko)
cellchat_ko <- aggregateNet(cellchat_ko)
levels(cellchat_ko@idents)            #查看细胞顺序
# vertex.receiver = c(3, 6)          #指定靶细胞的索引 内皮3 中性粒 8
cellchat_ko@netP$pathways   

# 保存数据
save(cellchat_sham, cellchat_mcao, cellchat_ko, file = "rds/cellchat_results_3subsets.rda")

# 合并对象
object.list <- list(Sham = cellchat_sham, MCAO = cellchat_mcao, KO = cellchat_ko)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 比较两个组的总体通讯数目和强度
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2,3))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2,3), measure = "weight")
p1 + p2 + ggsci::scale_fill_npg()

# 查看细胞群体的通讯变化（中性粒细胞和内皮细胞的互作）
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,2,3),
                          title.name = "Differential number of intractions - MCAO vs KO")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2,3),
                          title.name = "Differential intraction strength - MCAO vs KO")

# 气泡图展示中性粒细胞和内皮细胞的互作
cellchat@meta$group <- cellchat@meta$samples
cellchat@meta <- cellchat@meta[cellchat@meta$group != "KO", ]
cellchat@meta$group <- factor(cellchat@meta$group)
table(cellchat@meta$group)
levels(cellchat@meta$group)
cellchat@meta$group <- droplevels(cellchat@meta$group)
valid_rows <- cellchat@meta$group %in% c("Sham", "MCAO")
cellchat@meta <- cellchat@meta[valid_rows, ]

netAnalysis_contribution(cellchat, signaling = "Neutrophil-Endothelial")

netVisual_bubble(cellchat,
                 sources.use = "Neutrophil",
                 targets.use = "Endothelial",
                 comparison = c("Sham", "MCAO"),
                 angle.x = 45)




# 信号通路比较
# pathways.show <- c("CCL")  # 替换为感兴趣的通路
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "OrRd")



# 1.总体比较：通讯数目与通讯强度差异
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1+ ggsci::scale_fill_npg() + p2 + ggsci::scale_fill_npg()

# 2.细胞亚群水平的通讯差异
### 2.1 细胞通讯差异网络图
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# 展示了样本间所有细胞亚群中的配受体对数目差异(左)及通讯概率差异(右)；
# 外周圆不同颜色代表不同细胞亚群，大小表示亚群的配受体对数目，圈越大，细胞间配受体对数目比值越大.
# 蓝线表示对照组通讯较强，红线表示比较组通讯较强。线越粗表示通讯变化程度越强。

### 2.2 细胞通讯差异热图
p3 <- netVisual_heatmap(cellchat,comparison = c(2,3))
p4 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(2,3))
p3 + p4

# 3.信号通路水平的通讯差异
### 3.1 组间富集信号通路差异条形图
##基于信息流或互作数对信号通路进行排序
p5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) #堆叠
p6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) #不堆叠
p5 + p6
# 当比较（LS）组与对照组（NL）的通路概率总和的比值＜0.95且差异pval＜0.05（秩和检验）时，
# 则该通路在对照组中的通讯强度显著增加（纵坐标为红色）；当比较组与对照组的通路概率总和的比值＞1.05且差异pval＜0.05时，
# 则该通路在比较组中的通讯强度显著增加（纵坐标为蓝色）；
# 纵坐标为黑色表示该通路在两组间没有差异。左侧为比例图，右侧为实际数值比对图。

### 3.2 传出信号通路水平热图
# complexheatmap is needed
library(CellChat)
library(ComplexHeatmap)
library(patchwork)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways,
                       object.list[[i+1]]@netP$pathways)
pathway.union

ccs <- netAnalysis_computeCentrality(cc_sham)
ccm <- netAnalysis_computeCentrality(cc_mcao)
save(ccs,ccm, file = "rds/cellchat-divided.rda")
object.list <- list(Sham = ccs,
                    MCAO = ccm)
# plot
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5,
                                        height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5,
                                        height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

### 3.3 传入信号通路水平热图
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

### 3.4 总体信号通路水平热图
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "all", #总体
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 5, height = 6,
                                        color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "all", #总体
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 5, height = 6,
                                        color.heatmap = "OrRd")
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))

# 纵坐标为信号通路，横坐标为细胞亚群，热图颜色代表信号强度，颜色越深，通讯越强。上侧和右侧的柱子是纵轴和横轴强度的累积。
# 左图NL组，右图LS组。

# 4.配受体对水平通讯差异
### 4.1 总配受体对概率差异气泡图
levels(cellchat@idents$joint) #查看细胞亚群
cellchat@group
levels(cellchat@meta$group)

netVisual_bubble(cellchat,
                 sources.use = "Neutrophil",
                 targets.use = c("Endothelial","Microglia"),
                 comparison = c(1,2),
                 angle.x = 45)

### 4.2 区分上下调配受体对 !!!! Important!!!

p7 <- netVisual_bubble(cellchat,
                       sources.use = 4,
                       targets.use = c(6),
                       comparison = c(1,2),
                       max.dataset = 2,
                       title.name = "Increased signaling in MCAO",
                       angle.x = 45,
                       remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息
p8 <- netVisual_bubble(cellchat,
                       sources.use = 4,
                       targets.use = c(6),
                       comparison = c(1,2),
                       max.dataset = 1,
                       title.name = "Decreased signaling in MCAO",
                       angle.x = 45,
                       remove.isolate = T) #Decreased为对照组通讯概率更强的配受体对信息
p7 + p8
# 在气泡图中，横坐标为细胞对，颜色区分样本；纵坐标为配受体。
# 气泡的大小表示p值，p值越小气泡越大。颜色表示通讯概率的大小。

# 5.单个/特定信号通路水平差异可视化
#使用网络图：
pathways.show <- c("CCL") #选择目标信号通路
weight.max <- getMaxWeight(object.list,
                           slot.name = c("netP"),
                           attribute = pathways.show) 
#控制不同数据集的边权重
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],
                      signaling = pathways.show,
                      layout = "circle",
                      edge.weight.max = weight.max[1],
                      edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

#使用热图：
pathways.show <- c("CCL")
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]],
                               signaling = pathways.show,
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ht7 <- ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# plot export
if (!dir.exists("cellchat")) dir.create("cellchat")
ls() %>% cat(sep = ",")
for (ht in list(ht1,ht2,ht3,ht4,ht5,ht6,ht7)) {
  pdf(file = paste0("cellchat/",ht,".pdf"),
      width = 10, height = 6.18)
  ht
  dev.off()
}

save.image("rds/done.rda")

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_heatmap(cellchat, signaling = "Neutrophil-Endothelial", comparison = c("Sham", "MCAO"))
str(cellchat)
netVisual_aggregate(cellchat,
                    sources.use = "Neutrophil",
                    targets.use = "Endothelial",
                    signaling = "",
                    comparison = c("Sham", "MCAO"))

communication_data <- subsetCommunication(cellchat, 
                                          sources.use = "Neutrophil", 
                                          targets.use = "Endothelial")
head(communication_data)
library(ggplot2)

# 绘制气泡图，x轴为配体，y轴为受体，气泡大小代表通讯强度
ggplot(communication_data, aes(x = ligand, y = receptor, size = prob, color = prob)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Neutrophil-Endothelial Communication",
       x = "Ligand",
       y = "Receptor",
       size = "Communication Strength",
       color = "Communication Strength") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 将Sham和MCAO通讯数据合并
communication_sham <- communication_data$Sham
communication_mcao <- communication_data$MCAO

# 添加组别标识
communication_sham$group <- "Sham"
communication_mcao$group <- "MCAO"

# 合并数据
communication_combined <- rbind(communication_sham, communication_mcao)

# 按interaction_name分组，计算MCAO和Sham的通讯强度差异
library(dplyr)
communication_diff <- communication_combined %>%
  group_by(interaction_name) %>%
  summarise(diff = diff(prob))
library(ggplot2)

# 绘制气泡图
ggplot(communication_diff, aes(x = interaction_name, y = diff, size = abs(diff), color = diff)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Sham vs MCAO: Neutrophil-Endothelial Communication",
       x = "Interaction",
       y = "Difference in Communication Strength (MCAO - Sham)",
       size = "Abs(Difference)",
       color = "Difference") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 图片导出
# 获取绘图历史
plots <- recordPlot()
if (!dir.exists("visual")) {
  dir.create("visual")
}
for (i in seq_along(plots)) {
  pdf(file = paste0("visual/plot", i, ".pdf"), width = 8, height = 6)  # 设置合适的比例（宽度8，高度6）
  replayPlot(plots[[i]])  # 重现保存的图像
  dev.off()  # 关闭当前图形设备
}

