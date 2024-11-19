#------可视化---------

load(file = 'meta_seu_obj_tran.RData')

# ① 各群marker气泡图

# 假设 cell_group 是存储细胞群体的列名
cell_groups <- unique(seu$celltype)  # 获取细胞群体
levels(cell_groups)
# [1] "Microglia"            "Neuron"               "Ependymal Cells"      "Astrocytes"          
# [5] "Oligodendrocytes"     "VSMCs"                "Choroid Plexus Cells" "Endothelial Cells"   
# [9] "Macrophages"          "Fibroblasts"          "Immune Cells"        
p2 <- FeaturePlot(seu,
                  features = c('Map2'),raster = F)
# 选择 3 个标记物为每个细胞群体（示例）
markers <- c("Tmem119","P2ry12","Cx3cr1",  # 小胶质 标记物
             'Snap25',"Gad1", "Syt1",  # 神经元 标记物
             'Ccdc153','Cfap126','Fam183b',  # 室管膜 标记物
             "Aqp4",'S100b','Aldh1l1',  # 星形胶质 标记物
             'Mog','Mag','Olig2',  # 少突胶质 标记物
             
             'Vtn','Lhfp','Myh11',  # 周细胞 标记物
             'Folr1', 'Ecrg4','Pcp4', # 脉络丛 标记物
             
             "Cldn5", "Ebf1",'Pglyrp1',  # 内皮细胞 标记物
        
            
             "Cd68", "Cd74","Lyz2", # Macro
             
             'Dcn','Col1a1','Col12a1',  # 成纤维 标记物
             
             'Cd3d',"Flt3", "Klrd1"  # T细胞 标记物
             
             )

# 使用 DotPlot 函数绘制气泡图，并设置彩虹色颜色条，同时调整表达值范围
p <- DotPlot(seu, features = markers, group.by = "celltype") + 
  scale_size(range = c(0.7, 10)) +  # 调整气泡的大小范围!!!!
  scale_color_gradientn(
    colors = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#FEE08B", "#FDAE61", "#D53E4F"), 
    limits = c(-1, 2.5)) +  # 设置表达值范围
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 旋转X轴标签以提高可读性
  scale_x_discrete(expand = c(0, 0)) + # 确保横轴没有额外的留白
  scale_x_discrete(expand = c(0.019, 0))

p

# 提取 y 轴的数据
data <- ggplot_build(p)
y_max <- max(data$data[[1]]$y)  # 获取气泡图的 y 最大值


# 创建细胞名称颜色条的背景
cell_colors <- colorRampPalette(
  c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#FEE08B", "#FDAE61", "#D53E4F"))(length(cell_groups))
# "#5E4FA2" "#4371B2" "#3C93B8" "#5BB6A9" "#81CCA4" "#ABDDA4" "#DCDE94" "#FDD582" "#FDB769" "#EC8159" "#D53E4F"
cell_colors <- cols
cell_groups_ordered <- levels(cell_groups)
# 为每个细胞覆盖三个 marker 的颜色条
p1 <- p + 
  geom_tile(data = data.frame(x = seq(2, 3 * length(cell_groups_ordered), by = 3),  # 确保颜色条在三个 marker 居中显示
                              y = rep(y_max + 1.5, length(cell_groups_ordered)),  # 向上移动颜色条，避免重叠
                              fill = cell_colors),
            aes(x = x, y = y, fill = fill), width = 3, height = 1) +  # 设置更高的颜色条
  scale_fill_identity() +  # 使用自定义颜色
  geom_text(data = data.frame(x = seq(2, 3 * length(cell_groups_ordered), by = 3),  # seq(1.5, ...) 调整为 seq(2, ...)：颜色条和标签的 x 起点整体向右平移了 0.5
                              y = rep(y_max + 1.5, length(cell_groups_ordered)),  # 与颜色条一起向上移动
                              label = cell_groups_ordered),
            aes(x = x, y = y, label = label), color = "white", size = 4, vjust = 0.5) +  # 在颜色条上添加细胞名称
  theme_classic() +  # 使用经典主题，保留横线
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.line.x = element_line(size = 0.9),  # 确保 x 轴横线显示
        axis.line.y = element_line(size=0.9),  # 确保 y 轴横线显示
        axis.ticks.length = unit(0.15, "cm"),  # 设置刻度线长度
        axis.ticks = element_line(size = 0.9),  # 加粗刻度线
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold",colour = 'black'),  # 加粗 x 轴字体
        axis.text.y = element_text(size = 12, face = "bold",colour = 'black'),  # 加粗 y 轴字体
        panel.grid = element_blank(),  # 移除所有内部虚线
        plot.margin = unit(c(0.2, 0.2, 0.5, 0.2), "cm"))  # 调整上下左右边距

p1

pdf(file = 'visual/cluster_marker.pdf',width = 18,height = 8)
p1
dev.off()
