rm(list=ls(all=T))
# Load required libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(readr)
library(janitor)

# ==== Read data files =====
# From processed_data folder
fticr_data <- read_csv("processed_data/Processed_FTICR_60923_Data.csv")
names(fticr_data)[1] = 'molecular_formula'
mol_properties <- read_csv("processed_data/Processed_FTICR_60923_Mol.csv")
mol_properties = mol_properties[,-1]
vk_abundance <- read_csv("processed_data/Relative_Abundance_FTICR_60923_VK_Class1.csv")
summary_properties <- read_csv("processed_data/Summary_FTICR_60923_Properties.csv")
metadata = read_csv('processed_data/Metadata_Ricketts.csv') %>%
  clean_names()

# From data_from_EMSL_data_Central folder
coordinates <- read_csv("data_from_EMSL_data_Central/Coordinates.csv")%>%
  mutate(
    proposal_id = as.character(proposal_id),
    sampling_set = as.character(sampling_set)
  )

coordinates = merge(coordinates,metadata, by = c('proposal_id','sampling_set'))

coordinates = coordinates %>% filter(coordinates$rep == 1)

coordinates$sample_name = gsub("M_1","M", coordinates$sample_name)
coordinates$sample_name = gsub("P_1","P", coordinates$sample_name)

# ===== Data processing for visualization =====
sample_info <- tibble(
  sample_name = colnames(fticr_data)[colnames(fticr_data) != "molecular_formula"]
) %>%
  separate(sample_name, into = c("proposal_id", "sample_set", "depth"), 
           sep = "_", remove = FALSE) %>%
  mutate(
    proposal_id = as.character(proposal_id),
    sample_set = as.character(sample_set)
   ) %>%
   left_join(coordinates, by = c("proposal_id", "sample_set" = "sampling_set", "sample_name"))


# Prepare data for NMDS analysis
fticr_matrix <- fticr_data %>%
  column_to_rownames("molecular_formula") %>%
  t()

# Join molecular formula data with properties for VK plots
formula_presence <- fticr_data %>%
  pivot_longer(
    cols = -molecular_formula,
    names_to = "sample_name",
    values_to = "present"
  ) %>%
  filter(present == 1) %>%
  left_join(mol_properties, by = "molecular_formula")

# Color palette for van Krevelen classes
vk_colors <- c(
  "Amino Sugar" = "#FF9AA2",
  "Carbohydrate" = "#FFDAC1",
  "Cond Hydrocarbon" = "#E2F0CB",
  "Lignin" = "#B5EAD7",
  "Lipid" = "#C7CEEA",
  "Other" = "#9AA6C4",
  "Protein" = "#F5B5FC",
  "Tannin" = "#BDB2FF",
  "Unsat Hydrocarbon" = "#FFCAB0"
)

# ==== 1. NMDS Plot ====

set.seed(123) # For reproducibility
nmds_result <- metaMDS(fticr_matrix, distance = "jaccard", trymax = 100)

# Extract NMDS coordinates and join with metadata
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites")) %>%
  rownames_to_column("sample_name") %>%
  left_join(sample_info, by = "sample_name")

# Create NMDS plot
nmds_plot <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = cultivar, shape = site)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "NMDS Analysis of FTICR Data",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       color = "Cultivar",
       shape = "Site") +
  theme(legend.position = "right") +
  coord_fixed()

print(nmds_plot)
  ggsave("plots/nmds_plot.png", nmds_plot, width = 8, height = 6, dpi = 300)

  # ==== 2. Bar plot of relative abundance with state and cultivar ====
  # ==== 2. Bar plot with proper grouping like the image ====
  
  # Prepare data with proper grouping
  top_vk_abundance_grouped <- vk_abundance %>%
    filter(str_detect(sample_name, "_TOP$")) %>%
    pivot_longer(
      cols = -sample_name,
      names_to = "class",
      values_to = "abundance"
    ) %>%
    mutate(class = factor(class)) %>%
    left_join(sample_info, by = "sample_name") %>%
    # Create the grouping variable
    mutate(
      group = interaction(site, cultivar, sep = "\n"),
      sample_id = str_extract(sample_name, "\\d+_\\d+")
    ) %>%
    # Remove any NAs
    filter(!is.na(site), !is.na(cultivar))
  
  # Create the grouped bar plot
  grouped_bar_plot <- ggplot(top_vk_abundance_grouped, 
                             aes(x = sample_id, y = abundance, fill = class)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = vk_colors) +
    facet_wrap(~ group, scales = "free_x", nrow = 1) +
    labs(
      title = "Relative Abundance of Van Krevelen Classes in TOP Samples",
      x = "",
      y = "Percentage (%)",
      fill = "VK Class"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.spacing = unit(1, "lines")
    )
  
  print(grouped_bar_plot)
  ggsave("plots/vk_abundance_top_grouped_final.png", grouped_bar_plot, width = 16, height = 6, dpi = 300)
  
  # Same for BTM samples
  btm_vk_abundance_grouped <- vk_abundance %>%
    filter(str_detect(sample_name, "_BTM$")) %>%
    pivot_longer(
      cols = -sample_name,
      names_to = "class",
      values_to = "abundance"
    ) %>%
    mutate(class = factor(class)) %>%
    left_join(sample_info, by = "sample_name") %>%
    mutate(
      group = interaction(site, cultivar, sep = "\n"),
      sample_id = str_extract(sample_name, "\\d+_\\d+")
    ) %>%
    filter(!is.na(site), !is.na(cultivar))
  
  grouped_bar_plot_btm <- ggplot(btm_vk_abundance_grouped, 
                                 aes(x = sample_id, y = abundance, fill = class)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = vk_colors) +
    facet_wrap(~ group, scales = "free_x", nrow = 1) +
    labs(
      title = "Relative Abundance of Van Krevelen Classes in BTM Samples",
      x = "",
      y = "Percentage (%)",
      fill = "VK Class"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right",
      panel.spacing = unit(1, "lines")
    )
  
  print(grouped_bar_plot_btm)
  ggsave("plots/vk_abundance_btm_grouped_final.png", grouped_bar_plot_btm, width = 16, height = 6, dpi = 300)

  
# ==== 3. Box plots for VK classes ====
  # Option 1: Box plots faceted by VK class, grouped by site-cultivar
  boxplot_by_class <- ggplot(top_vk_abundance_grouped, 
                             aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    facet_wrap(~ class, scales = "free_y", ncol = 3) +
    labs(
      title = "Distribution of Van Krevelen Classes in TOP Samples",
      subtitle = "Grouped by Site and Cultivar",
      x = "Site - Cultivar",
      y = "Relative Abundance (%)",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
  
  print(boxplot_by_class)
  ggsave("plots/vk_boxplot_by_class_top.png", boxplot_by_class, width = 14, height = 10, dpi = 300)
  
  # Option 3: Single plot with all groups and classes
  boxplot_combined <- ggplot(top_vk_abundance_grouped, 
                             aes(x = class, y = abundance, fill = group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = "Distribution of Van Krevelen Classes in TOP Samples",
      subtitle = "Comparison across Site and Cultivar",
      x = "VK Class",
      y = "Relative Abundance (%)",
      fill = "Site - Cultivar"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right"
    )
  
  print(boxplot_combined)
  ggsave("plots/vk_boxplot_combined_top.png", boxplot_combined, width = 14, height = 6, dpi = 300)
  
  # Same for BTM samples

  boxplot_by_class <- ggplot(btm_vk_abundance_grouped, 
                             aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    facet_wrap(~ class, scales = "free_y", ncol = 3) +
    labs(
      title = "Distribution of Van Krevelen Classes in BTM Samples",
      subtitle = "Grouped by Site and Cultivar",
      x = "Site - Cultivar",
      y = "Relative Abundance (%)",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
  
  print(boxplot_by_class)
  ggsave("plots/vk_boxplot_by_class_btm.png", boxplot_by_class, width = 14, height = 10, dpi = 300)
  
  # Option 3: Single plot with all groups and classes
  boxplot_combined <- ggplot(btm_vk_abundance_grouped, 
                             aes(x = class, y = abundance, fill = group)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = "Distribution of Van Krevelen Classes in BTM Samples",
      subtitle = "Comparison across Site and Cultivar",
      x = "VK Class",
      y = "Relative Abundance (%)",
      fill = "Site - Cultivar"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right"
    )
  
  print(boxplot_combined)
  ggsave("plots/vk_boxplot_combined_btm.png", boxplot_combined, width = 14, height = 6, dpi = 300)
  
  
  #  Compare TOP vs BTM in one plot
  all_samples_boxplot <- bind_rows(
    top_vk_abundance_grouped %>% mutate(depth = "TOP"),
    btm_vk_abundance_grouped%>% mutate(depth = "BTM")
  )
  
  depth_comparison_boxplot <- ggplot(all_samples_boxplot, 
                                     aes(x = class, y = abundance, fill = depth)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    facet_wrap(~ group, nrow = 1) +
    scale_fill_manual(values = c("TOP" = "#E69F00", "BTM" = "#56B4E9")) +
    labs(
      title = "Van Krevelen Classes: TOP vs BTM Comparison",
      subtitle = "By Site and Cultivar",
      x = "VK Class",
      y = "Relative Abundance (%)",
      fill = "Depth"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "right"
    )
  
  print(depth_comparison_boxplot)
  ggsave("plots/vk_boxplot_depth_comparison.png", depth_comparison_boxplot, width = 18, height = 6, dpi = 300)  
  
  
  # ===== 4. Visualize AImod, NOSC, and GFE across samples =====
# Join sample info with summary properties using the same grouping approach
  sample_properties <- summary_properties %>%
    mutate(
      sample_name = paste(Proposal_ID, Sample_number, Depth, sep = "_"),
      Proposal_ID = as.character(Proposal_ID),
      Sample_number = as.character(Sample_number)
    ) %>%
    # Join with sample_info instead of coordinates to get site and cultivar
    left_join(sample_info, by = "sample_name") %>%
    # Create the same grouping variable we used for VK plots
    mutate(
      group = interaction(site, cultivar, sep = "\n"),
      sample_id = str_extract(sample_name, "\\d+_\\d+")
    ) %>%
    filter(!is.na(site), !is.na(cultivar))
  
  # Prepare data for visualization
  properties_long <- sample_properties %>%
    select(sample_name, sample_id, Mean_AImod, Mean_NOSC, Mean_GFE, site, cultivar, group, Depth) %>%
    pivot_longer(
      cols = c(Mean_AImod, Mean_NOSC, Mean_GFE),
      names_to = "property",
      values_to = "value"
    ) %>%
    mutate(
      property = factor(property,
                        levels = c("Mean_AImod", "Mean_NOSC", "Mean_GFE"),
                        labels = c("Mean Modified Aromaticity Index", "Mean NOSC", "Mean Gibbs Free Energy"))
    )
  
  # Boxplot visualization for molecular properties by site-cultivar groups
  properties_boxplot <- ggplot(properties_long,
                               aes(x = group, y = value, fill = Depth)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    #geom_jitter(width = 0.1, alpha = 0.5, size = 1.5, 
              #  position = position_jitterdodge(dodge.width = 0.8)) +
    facet_wrap(~property, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("TOP" = "#E69F00", "BTM" = "#56B4E9")) +
    theme_bw() +
    labs(
      title = "Distribution of Molecular Properties",
      subtitle = "By Site, Cultivar, and Depth",
      x = "Site - Cultivar",
      y = "Value",
      fill = "Depth"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  print(properties_boxplot)
  ggsave("plots/Mean_molecular_properties_grouped_boxplot.png", properties_boxplot, width = 12, height = 10, dpi = 300)
  
  # Properties by depth
  properties_by_depth <- ggplot(properties_long,
                                aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    facet_grid(property ~ Depth, scales = "free_y") +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    theme_bw() +
    labs(
      title = "Molecular Properties: TOP vs BTM Comparison",
      subtitle = "By Site and Cultivar",
      x = "Site - Cultivar",
      y = "Value",
      fill = "Group"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"  # Remove legend since x-axis shows the groups
    )
  
  print(properties_by_depth)
  ggsave("plots/Mean_molecular_properties_by_depth.png", properties_by_depth, width = 14, height = 10, dpi = 300)
  
  # Summary statistics table
  properties_summary <- properties_long %>%
    group_by(group, Depth, property) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    arrange(property, group, Depth)
  
  print(properties_summary)
  write_csv(properties_summary, "processed_data/molecular_properties_summary_by_group.csv")

# ==== 5. Unique and shared molecular formulas with groupings ====
  
# Analyze unique and shared formulas between TOP and BTM samples by group
  formula_group_analysis <- fticr_data %>%
    pivot_longer(
      cols = -molecular_formula,
      names_to = "sample_name",
      values_to = "present"
    ) %>%
    filter(present == 1) %>%
    # Join with sample_info to get grouping information
    left_join(sample_info, by = "sample_name") %>%
    filter(!is.na(site), !is.na(cultivar)) %>%
    mutate(
      depth = str_extract(sample_name, "(TOP|BTM)$"),
      group = interaction(site, cultivar, sep = " - "),
      sample_id = str_extract(sample_name, "\\d+_\\d+")
    ) %>%
    # Analyze within each group
    group_by(group, molecular_formula) %>%
    summarize(
      n_samples = n(),
      in_top = any(depth == "TOP"),
      in_btm = any(depth == "BTM"),
      category = case_when(
        in_top & in_btm ~ "Shared",
        in_top ~ "TOP only",
        in_btm ~ "BTM only"
      ),
      .groups = "drop"
    )
  
  # Join with molecular properties
  formula_group_properties <- formula_group_analysis %>%
    left_join(mol_properties, by = "molecular_formula")
  
  # Create summary of unique/shared formulas by group
  formula_group_summary <- formula_group_analysis %>%
    group_by(group, category) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(
      total = sum(n),
      percentage = n / total * 100
    ) %>%
    ungroup()
  
  # 1. Bar plot showing distribution by group
  formula_distribution_plot <- ggplot(formula_group_summary, 
                                      aes(x = group, y = n, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Distribution of Molecular Formulas by Depth and Group",
      subtitle = "Unique vs Shared Formulas between TOP and BTM samples",
      x = "Site - Cultivar",
      y = "Number of Formulas",
      fill = "Formula Category"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right"
    ) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), size = 3)
  
  print(formula_distribution_plot)
  ggsave("plots/formula_distribution_by_group.png", formula_distribution_plot, 
         width = 12, height = 6, dpi = 300)
  
  # 3. Analysis of formula properties by category and group
  formula_properties_analysis <- formula_group_properties %>%
    filter(!is.na(HC), !is.na(OC)) %>%
    select(group, category, molecular_formula, HC, OC, AImod, NOSC) %>%
    pivot_longer(
      cols = c(HC, OC,  AImod, NOSC),
      names_to = "property",
      values_to = "value"
    )
  
  
  # 4. Van Krevelen plot colored by formula category
  vk_plot_by_category <- ggplot(formula_group_properties %>% filter(!is.na(HC), !is.na(OC)), 
                                aes(x = OC, y = HC, color = category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    facet_wrap(~ group, nrow = 1) +
    scale_color_brewer(palette = "Set2") +
    labs(
      title = "Van Krevelen Diagram by Formula Category",
      subtitle = "Unique vs Shared Formulas between TOP and BTM samples",
      x = "O:C Ratio",
      y = "H:C Ratio",
      color = "Formula Category"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
  
  print(vk_plot_by_category)
  ggsave("plots/vk_plot_by_formula_category.png", vk_plot_by_category, 
         width = 16, height = 6, dpi = 300)
  
 # 5. Summary statistics table
  formula_stats_summary <- formula_group_summary %>%
    pivot_wider(names_from = category, values_from = c(n, percentage), 
                values_fill = 0) %>%
    arrange(group)
  
  print(formula_stats_summary)
  write_csv(formula_stats_summary, "processed_data/formula_category_summary_by_group.csv")
  
  # 6. Overall comparison across all groups (original analysis but with group info)
  overall_formula_analysis <- fticr_data %>%
    pivot_longer(
      cols = -molecular_formula,
      names_to = "sample_name",
      values_to = "present"
    ) %>%
    filter(present == 1) %>%
    mutate(depth = str_extract(sample_name, "(TOP|BTM)$")) %>%
    group_by(molecular_formula) %>%
    summarize(
      n_samples = n(),
      in_top = any(depth == "TOP"),
      in_btm = any(depth == "BTM"),
      category = case_when(
        in_top & in_btm ~ "Shared",
        in_top ~ "TOP only",
        in_btm ~ "BTM only"
      ),
      .groups = "drop"
    )
  
  overall_summary <- overall_formula_analysis %>%
    count(category) %>%
    mutate(percentage = n / sum(n) * 100)
  
  # Overall pie chart
  overall_pie <- ggplot(overall_summary, aes(x = "", y = n, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Set2") +
    geom_text(aes(label = paste0(category, "\n", n, " (", round(percentage, 1), "%)")), 
              position = position_stack(vjust = 0.5)) +
    theme_void() +
    labs(
      title = "Overall Distribution of Molecular Formulas",
      subtitle = "All samples: TOP vs BTM comparison",
      fill = "Category"
    ) +
    theme(legend.position = "bottom")
  
  print(overall_pie)
  ggsave("plots/overall_formula_distribution.png", overall_pie, 
         width = 8, height = 8, dpi = 300)  