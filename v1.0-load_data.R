rm(list = ls())
library(Seurat)
library(tidyverse)
library(magrittr)
# library(readxl)

# load the seurat object
seurat_obj <- readRDS('GSE233208_Human_visium_ADDS_seurat_processed.rds')
str(seurat_obj)
dim(seurat_obj)

# load the patient meta-data table from Supplementary Table 1
library(readxl)
meta_df <- read_excel("SupplementaryTables.xlsx", 
                                  sheet = "Supp. Table 1 ")
# View(meta_df)
meta_df <- dplyr::rename(meta_df, Sample = Seq.ID)

# meta_df <- read.delim("adds_patient_data.txt", sep='\t', header=1) %>%
#   dplyr::rename( Sample = Seq.ID)

# identify the common columns between the seurat object and the meta_df,
# and remove them from the seurat object
common_cols <- intersect(colnames(seurat_obj@meta.data), colnames(meta_df))
common_cols <- common_cols[common_cols != 'Sample'] # keep this one for the join operation
seurat_obj@meta.data %<>% select(-all_of(common_cols))

# join the two tables based on the Sample name to create an updated meta-data table
updated_meta <- dplyr::left_join(
  seurat_obj@meta.data, 
  meta_df, 
  by = 'Sample'
)
seurat_obj@meta.data <- updated_meta
colnames(seurat_obj)
rownames(seurat_obj@meta.data) <- colnames(seurat_obj)
# check the updated Diagnosis column to see if the numbers in each group are correct
patient_meta <- seurat_obj@meta.data %>%
  select(c(Sample, Diagnosis)) %>% distinct()
table(patient_meta$Diagnosis)

# in Github the author say it should be
#  AD   AD_DS Control earlyAD 
# 10      10      10       9 

# but my output is 
#     AD   AD_DS Control earlyAD 
# 8       8      15       8 

# Let's check
table(meta_df$Diagnosis...13)  # is right
table(meta_df$Diagnosis...34)  # more weird,      
# AD   AD_DS Control 
# 18      14       4 

table(seurat_obj$Diagnosis)

# use dignosis...13 as the right column
unique(seurat_obj$Sample)
unique(meta_df$Sample)

names(seurat_obj)
seurat_obj@images
SpatialDimPlot(seurat_obj)

slice1 <- seurat_obj@images$slice1
str(slice1)
Seurat::readPNG(slice1)

"MERTK" %in% rownames(seurat_obj)
p <- SpatialFeaturePlot(seurat_obj, features = "MERTK", ncol = 5)
ggsave("mertk.pdf", p, width = 40, height = 32)

table(seurat_obj$SampleID)
table(seurat_obj$Slide)
table(seurat_obj$Sectioning)
table(seurat_obj$Sample)
table(seurat_obj$combined_id)

seu <- subset(seurat_obj, SampleID == "Dec_20_2021_Human8")

names(seurat_obj)
str(seurat_obj)

slice <- "slice1.6.1"
# 提取当前slice对应的barcode
seu <- seurat_obj@images[[slice]]
str(seu)


slice_barcodes <- colnames(seurat_obj)[WhichCells(seurat_obj, expression = seurat_obj@images[[slice]]@coordinates$tissue == 1)]
slice_barcodes
# 子集Seurat对象
slice_obj <- subset(seurat_obj, cells = slice_barcodes)


# seurat_obj@images[[slice]] 提取出来的是一个 VisiumV1 对象，
# 本质上只是保存了切片图像、spot坐标等信息，而不是一个完整的 Seurat 对象，
# 所以不能直接用 SpatialDimPlot() 或 SpatialFeaturePlot()。

# 你应该从 meta.data 中提取当前slice的 barcode，
# 然后用 subset() 得到该slice对应的 Seurat 子对象，这样才能配合 SpatialDimPlot()：


# 获取该切片对应的barcode（坐标的rownames就是barcode）
barcodes <- rownames(seurat_obj@images[[slice]]@coordinates)
barcodes

# 从原始对象中子集提取该切片对应的Seurat对象
slice_obj <- subset(seurat_obj, cells = barcodes)
slice_obj@active.ident
Idents(slice_obj) <- slice_obj$annotation
unique(slice_obj$Sample)

# 画图（你可以改成 SpatialFeaturePlot 等）
p <- SpatialDimPlot(slice_obj, images = slice, pt.size.factor = 1.6) +
  ggtitle(slice)

# 显示或保存
print(p)



