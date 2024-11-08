# 加载必要的库
library(ggplot2)
library(dplyr)

# 假设 up_go_enrich_df 已经存在
# 按 p.adjust 升序排列并选择前几行
top_descriptions <- up_go_enrich_df %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust, Count) %>%
  head(50)

# 显示结果
print(top_descriptions$Description)

# 假设 up_go_enrich_df 已经创建
# 筛选感兴趣的 GO terms
selected_terms <- c("phagocytic vesicle", "protein tyrosine kinase", "cell adhesion molecule binding", 
                    "very-low-density lipoprotein", "protein-lipid complex",
                    "cytokine binding", "lipoprotein particle binding")

# 筛选这些特定的 terms 并加上排名靠前的其他 GO terms
# top_terms <- up_go_enrich_df %>%
#   filter(Description %in% selected_terms | rank(p.adjust) <= 10) %>%
#   arrange(p.adjust)

top_terms <- up_go_enrich_df %>%
    filter(Description %in% selected_terms) %>%
    arrange(p.adjust)
# 加载必要的库
library(ggplot2)

# 绘制挑选后的 GO 富集结果气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = rgb(224, 102, 99, maxColorValue = 255), high = "lightblue") +  # 红蓝配色，显著性高的为红色
  theme_minimal(base_size = 15) +
  labs(
    title = "Selected GO Terms Enrichment Bubble Plot",
    x = "-log10(P.adjust)",
    y = "GO Term"
  ) +
  theme(
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text.y = element_text(size = 12),  # 调整 y 轴标签字体大小
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_size(range = c(3, 10))  # 控制气泡大小范围


# 筛选 ----------------------------------------------------------------------
# 加载必要的包
library(clusterProfiler)
library(enrichplot)

# 选取包含感兴趣的 GO terms 的描述列表
selected_terms <- c("defense response to virus", "regulation of innate immune response", 
                    "activation of innate immune response", "tumor necrosis factor superfamily cytokine production", 
                    "positive regulation of innate immune response", "leukocyte migration", 
                    "myeloid leukocyte migration", "granulocyte migration", 
                    "positive regulation of inflammatory response", "innate immune response-activating signaling pathway",
                    "phagocytic vesicle")  # 包含 "phagocytic vesicle"

# 使用原始 enrichResult 对象进行筛选，保留感兴趣的 GO terms
top_terms_filtered <- up_go_enrich[up_go_enrich$Description %in% selected_terms, ]

# 绘制气泡图
dotplot(top_terms_filtered, showCategory = length(selected_terms), split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  ggtitle("Selected GO Terms Enrichment for Upregulated Genes") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # 去除背景网格线
    strip.text = element_text(size = 12),  # 调整分面标签字体大小
    plot.title = element_text(hjust = 0.5)
  )


# 加载必要的库
library(ggplot2)
library(dplyr)

# 假设 `top_terms_filtered` 已包含筛选后的 GO terms 数据
# 其中应包括 Description, p.adjust 和 Count 列
top_terms_filtered <- up_go_enrich_df %>%
  filter(Description %in% c("defense response to virus", "regulation of innate immune response", 
                            "activation of innate immune response", "tumor necrosis factor superfamily cytokine production", 
                            "positive regulation of innate immune response", "leukocyte migration", 
                            "myeloid leukocyte migration", "granulocyte migration", 
                            "positive regulation of inflammatory response", "innate immune response-activating signaling pathway",
                            "phagocytic vesicle")) %>%
  arrange(p.adjust)

# 绘制气泡图
ggplot(top_terms_filtered, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "#E06663", high = "blue") +  # 自定义配色，低 p 值为红色，高 p 值为蓝色
  theme_minimal(base_size = 15) +
  labs(
    title = "Selected GO Terms Enrichment Bubble Plot",
    x = "-log10(P.adjust)",
    y = "GO Term"
  ) +
  theme(
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text.y = element_text(size = 10),  # 调整 y 轴标签字体大小
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_size(range = c(3, 8))  # 控制气泡大小范围

