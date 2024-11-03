table(done_annotated$orig.ident, done_annotated$celltype)
head(done_annotated)
library(Seurat)
library(tidyverse)
# 查看每个样本中microglia的占比
microglia_prop <- done_annotated@meta.data %>%
  subset(celltype == "Microglia") %>%
  group_by(orig.ident) %>%
  summarise(microglia_count = n()) %>%
  mutate(total_cells = table(done_annotated$orig.ident)[orig.ident],
         proportion = microglia_count / total_cells)
