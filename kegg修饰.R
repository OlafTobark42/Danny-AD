library(stringr)

# 去掉后缀
up_kegg_enrich@result$Description <- str_remove(up_kegg_enrich@result$Description, " - Mus musculus \\(house mouse\\)")

# 绘制紧凑的 dotplot
kegg_up <- dotplot(up_kegg_enrich, showCategory = 10) + 
  ggtitle("Upregulated Gene KEGG Enrichment") +
  theme_minimal(base_size = 15) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  ) +
  
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +  # 控制文字换行
  theme(
    panel.grid = element_blank(),  # 去除背景网格线
    axis.text.y = element_text(size = 10), # 调整y轴文本的大小
    plot.title = element_text(hjust = 0.5), # 标题居中
    axis.title = element_text(size = 12)
  )

ggsave("visual/micro_kegg_ori.pdf", kegg_up,
       width = 7, height = 5)


dotplot(up_kegg_enrich, showCategory = 10) + 
  ggtitle("Upregulated Gene KEGG Enrichment") +
  theme_minimal(base_size = 15) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +  # 控制文字换行
  theme(
    axis.text.y = element_text(size = 10), # 调整y轴文本的大小
    plot.title = element_text(hjust = 0.5), # 标题居中
    axis.title = element_text(size = 12)
  ) +
  # coord_flip(clip = "off") + # 如果要更紧凑可以关闭坐标系裁剪
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.13))) # 调整 x 轴间距
 
  