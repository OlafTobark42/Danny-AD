
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


#  步骤二：添加 slice 名字到每个spot（cell）的 metadata 中 -------------------------------

# 创建barcode到slice的映射
barcode_to_slice <- lapply(names(seurat_obj@images), function(slice) {
  barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
  data.frame(barcode = barcodes, slice = slice, stringsAsFactors = FALSE)
})
# 先生成 barcode → slice 映射的命名向量
barcode_to_slice_vec <- setNames(barcode_to_slice$slice, barcode_to_slice$barcode)

# 然后赋值给 Seurat 对象
seurat_obj$slice <- barcode_to_slice_vec[rownames(seurat_obj)]


barcode_to_slice <- do.call(rbind, barcode_to_slice)

# 匹配到meta.data
seurat_obj$slice <- barcode_to_slice$slice[match(rownames(seurat_obj), barcode_to_slice$barcode)]

head(rownames(seurat_obj))
head(barcode_to_slice$barcode)


# 用 colnames(seurat_obj) 获取 barcode
shared_barcodes <- intersect(colnames(seurat_obj), barcode_to_slice$barcode)

# 检查是否有共通 barcode
length(shared_barcodes)  # 应该大于0，说明有正确匹配

# 重新构造命名向量，并赋值
barcode_slice_matched <- setNames(barcode_to_slice$slice, barcode_to_slice$barcode)[shared_barcodes]

# 添加到 metadata 中
seurat_obj$slice <- NA  # 初始化
seurat_obj$slice[shared_barcodes] <- barcode_slice_matched

table(seurat_obj$Sample, seurat_obj$slice)

# 步骤三：可选，画图时自动添加 sample / diagnosis 分组信息 ----------------------------------
SpatialDimPlot(seurat_obj, split.by = "sample")  # 或者 split.by = "Diagnosis...13"






