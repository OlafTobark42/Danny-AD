
# 第一步：筛选满足条件的 slice -------------------------------------------------------


# 只保留 Control 和 earlyAD 的 slice
target_slices <- slice_info$slice[slice_info$diagnosis %in% c("Control", "earlyAD")]
target_slices


# 第二步：逐个绘图 & 保存单独 PDF -----------------------------------------------------
library(Seurat)
library(ggplot2)
library(viridis)
dir.create("visual")

# 创建一个用于合并的空列表
plot_list <- list()

# 循环遍历 slice 绘图并保存
for (slice in target_slices) {
  barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
  slice_obj <- subset(seurat_obj, cells = barcodes)
  
  p <- SpatialFeaturePlot(slice_obj, features = "MERTK", 
                          images = slice, pt.size.factor = 1.6) +
    ggtitle(paste("MERTK -", slice, "-", slice_info$diagnosis[slice_info$slice == slice]))
  
  # 单独保存
  ggsave(filename = paste0("visual/MERTK_", 
                           slice_info$diagnosis[slice_info$slice == slice],
                           "_",slice, ".pdf"), 
         plot = p, width = 8, height = 6)
  
  # 收集用于合并
  plot_list[[slice]] <- p
}


# 第三步：拼成一个总 PDF，每行一个组 -----------------------------------------------------
library(patchwork)

# 按 Diagnosis 分组拼图
control_slices <- slice_info$slice[slice_info$diagnosis == "Control"]
earlyad_slices <- slice_info$slice[slice_info$diagnosis == "earlyAD"]

# 提取并组合成每一行
row1 <- wrap_plots(plot_list[control_slices], nrow = 1)
row2 <- wrap_plots(plot_list[earlyad_slices], nrow = 1)

# 拼成一个图（上下两行）
final_plot <- row1 / row2 + plot_layout(heights = c(1, 1))

# 保存成大 PDF
ggsave("visual/MERTK_all_Control_vs_earlyAD.pdf", plot = final_plot, width = 20, height = 16)


# 修改绘图收集代码（保留两套版本）： -------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

# 初始化列表
plot_list_with_title <- list()
plot_list_no_title <- list()

# 提取目标 slice
target_slices <- slice_info$slice[slice_info$diagnosis %in% c("Control", "earlyAD")]

for (slice in target_slices) {
  barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
  slice_obj <- subset(seurat_obj, cells = barcodes)
  
  # 获取诊断信息
  diagnosis <- slice_info$diagnosis[slice_info$slice == slice]
  
  # 带标题的图（用于单独保存）
  p1 <- SpatialFeaturePlot(slice_obj, features = "MERTK", images = slice, pt.size.factor = 1.6) +
    ggtitle(paste("MERTK -", slice, "-", diagnosis))
  
  # 无标题图（用于拼图）
  p2 <- SpatialFeaturePlot(slice_obj, features = "MERTK", images = slice, pt.size.factor = 1.6) +
    theme(plot.title = element_blank())
  
  # 保存单图
  ggsave(filename = paste0("visual/MERTK_", 
                           slice_info$diagnosis[slice_info$slice == slice],
                           "_",slice, ".pdf"),  
         plot = p1, width = 6, height = 5)
  
  # 加入拼图列表
  plot_list_with_title[[slice]] <- p1
  plot_list_no_title[[slice]] <- p2
}


# 拼图（无标题版）+ 添加分组标签 --------------------------------------------------------
# 分组提取
control_slices <- slice_info$slice[slice_info$diagnosis == "Control"]
earlyad_slices <- slice_info$slice[slice_info$diagnosis == "earlyAD"]

# 各行拼接
row1 <- wrap_plots(plot_list_no_title[control_slices], nrow = 1)
row2 <- wrap_plots(plot_list_no_title[earlyad_slices], nrow = 1)

# 拼接整体图 + 加大标题
final_plot <- row1 / row2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Spatial expression of MERTK in Control and earlyAD",
    tag_levels = list(c("Control", "earlyAD")),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# 保存为总 PDF
ggsave("visual/MERTK_merged_Control_vs_earlyAD.pdf", plot = final_plot, width = 20, height = 10)


# 修改后的拼图代码（统一 legend）： ----------------------------------------------------

library(patchwork)

# 提取无标题版本图
row1 <- wrap_plots(plot_list_no_title[control_slices], nrow = 1)
row2 <- wrap_plots(plot_list_no_title[earlyad_slices], nrow = 1)

# 拼接总图并共享图例
final_plot <- (row1 / row2) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "Spatial expression of MERTK",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
final_plot
# 用 patchwork 标签指定每一行
final_plot <- wrap_plots(
  list(row1, row2),
  ncol = 1,
  guides = "collect",
  heights = c(1, 1),
  tag_level = "new"
) +
  plot_annotation(
    title = "Spatial expression of MERTK",
    tag_levels = list(c("Control", "earlyAD")),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
ggsave("MERTK_merged_Control_earlyAD_shared_legend.pdf", plot = final_plot, width = 20, height = 10)



# 统一比例尺 -------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

# 全局表达范围，统一颜色 scale（默认 log-normalized）
mertk_values <- FetchData(seurat_obj, vars = "MERTK")
min_mertk <- min(mertk_values, na.rm = TRUE)
max_mertk <- max(mertk_values, na.rm = TRUE)

# 图列表
plot_list_for_merge <- list()

# 筛选目标 slices
target_slices <- slice_info$slice[slice_info$diagnosis %in% c("Control", "earlyAD")]

for (slice in target_slices) {
  barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
  slice_obj <- subset(seurat_obj, cells = barcodes)
  diagnosis <- slice_info$diagnosis[slice_info$slice == slice]
  
  # 绘图，统一颜色、去图例，加标题（仅 slice 名）
  p <- SpatialFeaturePlot(
    slice_obj, features = "MERTK", images = slice,
    min.cutoff = min_mertk, max.cutoff = max_mertk,
    pt.size.factor = 1.6
  ) +
    ggtitle(slice) +  # 仅加 slice ID
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "none"  # 去掉图例
    ) +
    viridis::scale_fill_viridis()
  p
  # 单图保存也保持统一色阶（可选）
  # 保存单图
  ggsave(filename = paste0("visual/MERTK_", 
                           slice_info$diagnosis[slice_info$slice == slice],
                           "_",slice, ".pdf"),  
         plot = p, width = 6, height = 5)
  # 收集用于合并
  plot_list_for_merge[[slice]] <- p
}



# 分组列表
control_slices <- slice_info$slice[slice_info$diagnosis == "Control"]
earlyad_slices <- slice_info$slice[slice_info$diagnosis == "earlyAD"]

# 分行合并图
row1 <- wrap_plots(plot_list_for_merge[control_slices], nrow = 1)
row2 <- wrap_plots(plot_list_for_merge[earlyad_slices], nrow = 1)

# 合并为总图
final_plot <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Spatial expression of MERTK (color scale unified across all slices)",
    tag_levels = list(c("Control", "earlyAD")),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# 保存
ggsave("MERTK_merged_Control_earlyAD_clean.pdf", plot = final_plot, width = 20, height = 10)


# 调整early AD位置 ------------------------------------------------------------

# 确保 diagnosis 信息是 factor，并排序（可选但推荐）
# slice_info$diagnosis <- factor(slice_info$diagnosis, levels = c("Control", "earlyAD"))

# 提前构建顺序明确的 slice 列表
control_slices <- slice_info$slice[slice_info$diagnosis == "Control"]
earlyad_slices <- slice_info$slice[slice_info$diagnosis == "earlyAD"]

# 取图
plots_control <- plot_list_for_merge[as.character(control_slices)]
plots_earlyad <- plot_list_for_merge[as.character(earlyad_slices)]

# 两行拼图
row1 <- wrap_plots(plots_control, nrow = 1)
row2 <- wrap_plots(plots_earlyad, nrow = 1)

# 拼成一个完整图
final_plot <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Spatial expression of MERTK (color scale unified across all slices)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# 保存
ggsave("MERTK_merged_Control_earlyAD_fixedOrder.pdf", plot = final_plot, width = 22, height = 10)

