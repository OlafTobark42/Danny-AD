# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Calculate the proportion of each cell type in each group
proportion_data <- seu@meta.data %>%
  group_by(RNA_snn_res.0.2, group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count))

# Plot the proportions
ggplot(proportion_data, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Group", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 按样本 ---------------------------------------------------------------------

# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Calculate the proportion of microglia in each sample
proportion_data <- seu@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(celltype == "Microglia") # Filter to keep only microglia proportions

# Add grouping information (AD or Aβ group) if stored in `group`
proportion_data <- proportion_data %>%
  left_join(seu@meta.data %>% select(orig.ident, group) %>% distinct(), by = "orig.ident")

# Perform statistical test (e.g., Wilcoxon test) to compare proportions between AD and Aβ groups
# Adjust `group` values based on your specific group labels
test_result <- wilcox.test(proportion ~ group, data = proportion_data, subset = (group %in% c("AD + WT", 
                                                                                              "Control")))

# Print test results
print(test_result)

# Plot microglia proportions across samples, grouped by AD and Aβ
ggplot(proportion_data, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(x = "Group", y = "Proportion of Microglia", title = "Microglia Proportion by Sample and Group") +
  theme_minimal()

