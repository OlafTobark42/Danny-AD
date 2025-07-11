seurat_obj <- readRDS("GSE233208_Human_visium_ADDS_seurat_processed.rds")
#  步骤一：构建每个 slice 对应的 Sample 和 Diagnosis 映射表 ------------------------------

# 初始化一个空的记录列表
slice_info <- data.frame(
  slice = character(),
  sample = character(),
  diagnosis = character(),
  stringsAsFactors = FALSE
)

# 遍历每个 slice
for (slice in names(seurat_obj@images)) {
  print(paste0("--- Process Slice: ", slice, " ---"))
  # 取出该切片的barcode（spot名）
  barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
  print(head(barcodes))
  # 查找这些 barcode 所属的 sample 和 diagnosis
  sample_names <- unique(seurat_obj$Sample[barcodes])
  
  diagnoses <- unique(seurat_obj$`Diagnosis...13`[barcodes])
  
  print(paste0(diagnoses, " belongs to ", slice))
  # 警告处理：如果不唯一，打印
  if (length(sample_names) != 1 | length(diagnoses) != 1) {
    warning(paste("Slice", slice, "has multiple sample or diagnosis values"))
  }
  
  # 存入表格
  slice_info <- rbind(slice_info, data.frame(
    slice = slice,
    sample = sample_names[1],
    diagnosis = diagnoses[1]
  ))
}

# 查看
head(slice_info)
table(slice_info$diagnosis)

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# 设置要绘图的基因
plot_genes <- c("CX3CR1", "P2RY12", "TMEM119", "MERTK")

# 获取目标切片名
slice_info <- slice_info[slice_info$diagnosis %in% c("Control", "earlyAD"), ]
target_slices <- slice_info$slice

# 筛选目标元数据（spot注释）
meta_df <- seurat_obj@meta.data
meta_df <- meta_df[meta_df$slice %in% target_slices, ]

# 创建保存目录
dir.create("individual", showWarnings = FALSE)

for (gene in plot_genes) {
  # 从 Seurat 提取表达值，仅限于目标 barcode
  gene_expr <- FetchData(seurat_obj, vars = gene)
  # gene_expr <- gene_expr[rownames(meta_df), , drop = FALSE]  # 顺序严格对齐
  
  # 统一色阶范围
  vmin <- min(gene_expr[[gene]], na.rm = TRUE)
  vmax <- max(gene_expr[[gene]], na.rm = TRUE)
  
  # 合并表达与元数据（行名即 barcode，自动对齐）
  full_df <- cbind(meta_df, expr = gene_expr[[gene]])
  
  for (slice in target_slices) {
    df_sub <- full_df[full_df$slice == slice, ]
    
    p <- ggplot(df_sub, aes(x = imagerow, y = imagecol, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_viridis(option = "A", limits = c(vmin, vmax)) +
      scale_y_reverse() +
      coord_fixed() +
      labs(title = paste0(gene, " - ", slice)) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none"
      )
    
    ggsave(filename = file.path("individual", paste0(gene, "_", slice, ".pdf")),
           plot = p, width = 4, height = 4)
  }
}



# old below ---------------------------------------------------------------


# 实现目标：干净的切片区域表达图
# 
# 你可以参考以下 改造版代码框架，专门针对 MERTK 这样某个 marker 的表达，
# 在没有背景图的基础上，只画组织区域的表达强度，并用统一配色（如 viridis）：
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
seurat_obj$MERTK <- NULL
# 提取 MERTK 表达量，合并到 meta.data 中
seurat_obj$MERTK <- FetchData(seurat_obj, vars = "MERTK")[, 1]

# 筛选目标切片
target_slices <- slice_info$slice[slice_info$diagnosis %in% c("Control", "earlyAD")]

# 在进入循环前，直接从 seurat_obj@meta.data 里筛选出对应的行
plot_df <- seurat_obj@meta.data[seurat_obj$slice %in% target_slices, ]

# 所需列
plot_df <- plot_df %>%
  mutate(barcode = rownames(plot_df)) %>%
  dplyr::select(slice,Sample, Diagnosis, imagerow, imagecol, MERTK)

# check if really only Control and earlyAD
unique(plot_df$Diagnosis);unique(plot_df$slice)

# 统一 color scale
color_max <- max(plot_df$MERTK, na.rm = TRUE)

# 创建 plot_list
plot_list_clean <- list()

for (cur_slice in target_slices) {
  cur_df <- plot_df[plot_df$slice == cur_slice, ]
  cur_diagnosis <- unique(cur_df$Diagnosis)
  
  p <- ggplot(cur_df, aes(x = imagerow, y = imagecol, color = MERTK)) +
    geom_point(size = 0.5) +
    scale_color_viridis(limits = c(0, color_max)) +
    coord_fixed() +
    scale_y_reverse() +  # 符合图像上下方向
    labs(title = cur_slice) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
  
  plot_list_clean[[cur_slice]] <- p
}

# p
# 获取诊断信息
slice_info <- plot_df %>%
  group_by(slice) %>%
  summarise(diagnosis = first(Diagnosis)) %>%
  arrange(diagnosis)

# 分组
control_samples <- slice_info$slice[slice_info$diagnosis == "Control"]
earlyad_samples <- slice_info$slice[slice_info$diagnosis == "earlyAD"]

row1 <- wrap_plots(plot_list_clean[control_samples], nrow = 1)
row2 <- wrap_plots(plot_list_clean[earlyad_samples], nrow = 1)

ggsave("MERTK_clean_nobg_control.pdf", row1, width = 20, height = 10)
ggsave("MERTK_clean_nobg_earlyAD.pdf", row2, width = 15, height = 10)
ggsave("MERTK_clean_nobg_control.png", row1, width = 20, height = 10)
ggsave("MERTK_clean_nobg_earlyAD.png", row2, width = 15, height = 10)

final_plot_clean <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Spatial expression of MERTK (clean, background-free)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave("MERTK_clean_nobg_control_earlyAD.pdf", final_plot_clean, width = 18, height = 10)



# 多个基因循环 ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# -------------------- 输入设置 -------------------- #
plot_genes <- c("CX3CR1", "P2RY12", "TMEM119", "MERTK")
slice_info <- slice_info[slice_info$diagnosis %in% c("Control", "earlyAD"), ]
target_slices <- slice_info$slice
plot_df <- seurat_obj@meta.data[seurat_obj$slice %in% target_slices, ]
plot_df$gene <- NULL  # 保证干净环境

# 先准备 plot_df：加上 barcode，提前过滤目标 slice
plot_df <- seurat_obj@meta.data[seurat_obj$slice %in% target_slices, ]
plot_df$barcode <- rownames(plot_df)

# -------------------- individual 绘图保存 -------------------- #
for (gene in plot_genes) {
  # 当前基因表达矩阵
  gene_expr <- FetchData(seurat_obj, vars = gene)
  gene_expr$barcode <- rownames(gene_expr)
  
  # 统一配色范围（放在这儿 OK）
  vmin <- min(gene_expr[[gene]], na.rm = TRUE)
  vmax <- max(gene_expr[[gene]], na.rm = TRUE)
  
  for (slice in target_slices) {
    df_sub <- plot_df[plot_df$slice == slice, ]
    
    # 用 barcode 连接表达量
    df_plot <- df_sub %>%
      left_join(gene_expr, by = "barcode")
    
    p <- ggplot(df_plot, aes(x = imagerow, y = imagecol, color = .data[[gene]])) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "viridis", limits = c(vmin, vmax)) +
      scale_y_reverse() +
      coord_fixed() +
      ggtitle(paste0(gene, " - ", slice)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 10))
    # p
    ggsave(filename = file.path("individual", paste0(gene, "_", slice, ".pdf")),
           plot = p, width = 4, height = 4)
  }
}


# -------------------- a. 每个基因的整合图 -------------------- #
for (gene in plot_genes) {
  plots <- list()
  for (slice in target_slices) {
    df_sub <- plot_df[plot_df$slice == slice, ]
    df_sub$expr <- FetchData(seurat_obj, vars = gene)[rownames(df_sub), gene]
    
    plots[[slice]] <- ggplot(df_sub, aes(x = imagerow, y = imagecol, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "viridis", limits = c(vmin, vmax)) +
      scale_y_reverse() +
      coord_fixed() +
      ggtitle(slice) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 8), legend.position = "none")
  }
  
  wrap_plot <- wrap_plots(plots, nrow = 2) +
    plot_annotation(title = paste0("Spatial expression of ", gene),
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(filename = file.path("combined_by_gene", paste0(gene, "_combined.pdf")),
         plot = wrap_plot, width = 20, height = 8)
}

# -------------------- b. 每张切片的四基因合图 -------------------- #
for (slice in target_slices) {
  plots <- list()
  for (gene in plot_genes) {
    df_sub <- plot_df[plot_df$slice == slice, ]
    df_sub$expr <- FetchData(seurat_obj, vars = gene)[rownames(df_sub), gene]
    
    plots[[gene]] <- ggplot(df_sub, aes(x = imagerow, y = imagecol, color = expr)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "viridis", limits = c(vmin, vmax)) +
      scale_y_reverse() +
      coord_fixed() +
      ggtitle(gene) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 8), legend.position = "none")
  }
  
  wrap_plot <- wrap_plots(plots, nrow = 1) +
    plot_annotation(title = paste0("Slice ", slice, " expression of selected genes"),
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(filename = file.path("combined_by_slice", paste0(slice, "_genes.pdf")),
         plot = wrap_plot, width = 16, height = 4)
}

