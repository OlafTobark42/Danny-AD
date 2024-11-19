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

# > microglia_prop
# # A tibble: 8 × 4
# orig.ident microglia_count total_cells proportion 
# <chr>                <int> <table[1d]> <table[1d]>
#   1 a1                    2613 8585        0.3043681  
# 2 a2                    2814 8555        0.3289305  
# 3 c1                    2904 9734        0.2983357  
# 4 c2                    2470 9332        0.2646807  
# 5 c3                    2272 7440        0.3053763  
# 6 k1                    2526 7492        0.3371596  
# 7 k2                    2208 8245        0.2677987  
# 8 k3                    2772 8493        0.3263864