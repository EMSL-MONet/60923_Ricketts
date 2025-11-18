rm(list=ls(all=T))
# Load required libraries
library(tidyverse)
library(vegan)
library(ggrepel)
library(patchwork)
library(readr)
library(janitor)

# ==== Read data files =====
geochem = read_csv('data_from_EMSL_data_central_11-17-25/processed_data/Soil_BioChemical_properties.csv')
geochem_top <- geochem %>%
  filter(str_detect(Sample_Name, "_TOP$"))

# From processed_data folder
fticr_data_full <- read_csv("processed_data/Processed_FTICR_60923_Data.csv")
names(fticr_data_full)[1] = 'molecular_formula'

# Filter to only TOP samples
top_columns <- c("molecular_formula", 
                 names(fticr_data_full)[str_detect(names(fticr_data_full), "_TOP$")])
fticr_data <- fticr_data_full %>%
  select(all_of(top_columns))

mol_properties <- read_csv("processed_data/Processed_FTICR_60923_Mol.csv")
mol_properties = mol_properties[,-1]

# Filter vk_abundance to only TOP samples
vk_abundance_full <- read_csv("processed_data/Relative_Abundance_FTICR_60923_VK_Class1.csv")
vk_abundance <- vk_abundance_full %>%
  filter(str_detect(sample_name, "_TOP$"))

# Filter summary_properties to only TOP samples
summary_properties_full <- read_csv("processed_data/Summary_FTICR_60923_Properties.csv")
summary_properties <- summary_properties_full %>%
  filter(str_detect(paste(Proposal_ID, Sample_number, Depth, sep = "_"), "_TOP$"))

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

# ===== Data processing for visualization (TOP samples only) =====
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

# Function to check which geochem variables are complete
check_complete_data <- function(x) {
  # Count missing values (NA, empty strings, common QC flags)
  qc_flags <- c("Failed_QC", "Below_Detection", "Below_Reporting_Limit", "Outlier", "", NA)
  missing_count <- sum(is.na(x) | x %in% qc_flags | trimws(as.character(x)) == "")
  return(missing_count == 0)
}

# identify variables with complete data
complete_vars <- geochem_top %>%
  # Remove non-variable columns (identifiers and flag columns)
  select(-c(Sample_Name, Proposal_ID, Sampling_Set, Core_Section)) %>%
  # Remove flag columns (they contain QC info, not environmental data)
  select(-contains("flag_")) %>%
  # Remove non-numeric extraction method columns
  select(-contains("extraction")) %>%
  select(-contains("_unit")) %>%
  # Keep only numeric columns
  select(where(is.numeric) | where(function(x) all(is.na(as.numeric(as.character(x))) == is.na(x)))) %>%
  # Check each column for completeness
  summarise(across(everything(), check_complete_data)) %>%
  # Get column names where all values are complete
  select(where(isTRUE)) %>%
  names()

# Prepare geochem data for envfit using compleye variables
geochem_for_envfit <- geochem_top %>%
  # Create sample_name to match your NMDS data format
  mutate(sample_name = paste(Proposal_ID, Sampling_Set, Core_Section, sep = "_")) %>%
  # Keep only complete variables
  select(sample_name, all_of(complete_vars)) %>%
  # Make sure sample names match exactly with your NMDS rownames
  filter(sample_name %in% rownames(fticr_matrix)) %>%
  arrange(match(sample_name, rownames(fticr_matrix)))

# Check if we have any variables left
if(length(complete_vars) == 0) {
  stop("No variables with complete data found. You may need to adjust the completeness criteria.")
}

# Extract just the environmental variables (remove sample_name)
env_vars <- geochem_for_envfit %>%
  select(-sample_name)

# Convert any character columns that should be numeric
env_vars <- env_vars %>%
  mutate(across(everything(), as.numeric))

# Remove any columns that became all NA after conversion
env_vars <- env_vars %>%
  select(where(~ !all(is.na(.))))

# Prepare data for NMDS analysis (TOP samples only)
fticr_matrix <- fticr_data %>%
  column_to_rownames("molecular_formula") %>%
  t()

# Join molecular formula data with properties for VK plots (TOP samples only)
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

# ==== 1. NMDS Plot (TOP samples only) ====
set.seed(123) # For reproducibility
nmds_result <- metaMDS(fticr_matrix, distance = "jaccard", trymax = 100)
# Extract NMDS coordinates and join with metadata
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites")) %>%
  rownames_to_column("sample_name") %>%
  left_join(sample_info, by = "sample_name")

# Create NMDS plot (TOP samples only)
nmds_plot <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, shape = cultivar, color = site)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = " ",
       subtitle = paste("Stress =", round(nmds_result$stress, 3)),
       shape = "Cultivar",
       color = "Site") +
  theme(legend.position = "right") +
  coord_fixed()
print(nmds_plot)
#ggsave("plots/nmds_plot_TOP_only.png", nmds_plot, width = 8, height = 6, dpi = 300)

# Run envfit analysis
envfit_result <- envfit(nmds_result, env_vars, permutations = 999, na.rm = TRUE)

# View the results
print(envfit_result)

# Extract significant variables (p < 0.05)
if(length(envfit_result$vectors$pvals) > 0) {
  sig_vars <- envfit_result$vectors$pvals < 0.05
  
  if(sum(sig_vars) > 0) {
    sig_envfit <- envfit_result$vectors$arrows[sig_vars, , drop = FALSE]
    sig_pvals <- envfit_result$vectors$pvals[sig_vars]
    sig_r2 <- envfit_result$vectors$r[sig_vars]
    
    # Create a dataframe for plotting arrows
    envfit_df <- data.frame(
      variable = rownames(sig_envfit),
      NMDS1 = sig_envfit[, 1],
      NMDS2 = sig_envfit[, 2],
      pval = sig_pvals,
      r2 = sig_r2
    ) %>%
      mutate(
        # Scale arrows for better visualization
        NMDS1_scaled = NMDS1 * 0.6,
        NMDS2_scaled = NMDS2 * 0.6,
        # Create significance labels
        significance = case_when(
          pval < 0.001 ~ "***",
          pval < 0.01 ~ "**", 
          pval < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
    
    # Enhanced NMDS plot with environmental vectors
    nmds_plot_envfit <- ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, shape = cultivar, color = site)) +
      geom_point(size = 3, alpha = 0.8) +
      # Add site ellipses
      stat_ellipse(aes(color = site), level = 0.95, linetype = "dashed", size = 1, alpha = 0.7) +
      # Add environmental vectors (arrows)
      geom_segment(data = envfit_df, 
                   aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
                   arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                   color = "black", size = 1, alpha = 0.8, inherit.aes = FALSE) +
      # Add variable labels
      geom_text_repel(data = envfit_df,
                      aes(x = NMDS1_scaled * 1.1, y = NMDS2_scaled * 1.1, 
                          label = paste0(variable, significance)),
                      color = "black", size = 3, inherit.aes = FALSE,
                      max.overlaps = 15) +
      theme_bw() +
      labs(subtitle = paste("Stress =", round(nmds_result$stress, 3), 
                            "| Arrows show significant soil variables"),
           shape = "Cultivar",
           color = "Site",
           caption = "* p<0.05, ** p<0.01, *** p<0.001") +
      theme(legend.position = "right") +
      coord_fixed()
    
    print(nmds_plot_envfit)
    ggsave("plots/nmds_plot_with_soil_envfit.png", nmds_plot_envfit, width = 12, height = 8, dpi = 300)
    
    # Create a summary table of significant variables
    envfit_summary <- envfit_df %>%
      arrange(pval) %>%
      mutate(
        r2_percent = round(r2^2 * 100, 1),
        pval_formatted = case_when(
          pval < 0.001 ~ "< 0.001",
          TRUE ~ as.character(round(pval, 3))
        )
      ) %>%
      select(variable, r2_percent, pval_formatted, significance)
    
    print("Significant soil variables driving site separation:")
    print(envfit_summary)
    
    # Save the summary
    write_csv(envfit_summary, "processed_data/envfit_soil_variables.csv")
    
  } else {
    cat("No significant environmental variables found (p < 0.05)\n")
  }
} else {
  cat("No environmental vectors could be fitted to the data\n")
}


# ==== 2. Bar plot of relative abundance for TOP samples ====
# Prepare data with proper grouping (TOP samples only)
top_vk_abundance_grouped <- vk_abundance %>%
  pivot_longer(
    cols = -sample_name,
    names_to = "class",
    values_to = "abundance"
  ) %>%
  mutate(class = factor(class)) %>%
  left_join(sample_info, by = "sample_name") %>%
  # Create variables for ordering within facets
  mutate(
    sample_id = str_extract(sample_name, "\\d+_\\d+"),
    # Create a combined variable for ordering: cultivar first, then sample ID
    cultivar_sample = paste(cultivar, sample_id, sep = "_")
  ) %>%
  # Remove any NAs
  filter(!is.na(site), !is.na(cultivar)) %>%
  # Order cultivars within each site
  arrange(site, cultivar, sample_id)

# Create factor levels to ensure proper ordering
top_vk_abundance_grouped <- top_vk_abundance_grouped %>%
  group_by(site) %>%
  mutate(cultivar_sample = factor(cultivar_sample, 
                                  levels = unique(cultivar_sample[order(cultivar, sample_id)]))) %>%
  ungroup()

# Create the grouped bar plot (TOP samples only) - Faceted by site
grouped_bar_plot <- ggplot(top_vk_abundance_grouped,
                           aes(x = cultivar_sample, y = abundance, fill = class)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = vk_colors) +
  facet_wrap(~ site, scales = "free_x", nrow = 1) +
  labs(
    title = "Relative Abundance of Van Krevelen Classes (TOP Samples Only)",
    subtitle = "Grouped by Site, with Cultivars Ordered Within Each Site",
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
ggsave("plots/vk_abundance_top_grouped_by_site.png", grouped_bar_plot, width = 16, height = 6, dpi = 300)

# ==== 3. Box plots for VK classes (TOP samples only) -  ====
# Box plots faceted by VK class, with cultivars grouped within sites
boxplot_by_class <- ggplot(top_vk_abundance_grouped,
                           aes(x = cultivar, y = abundance, fill = cultivar)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  facet_grid(class ~ site, scales = "free_y") +
  labs( x = " ",
    y = "Relative Abundance (%)",
    fill = "Cultivar"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  )
print(boxplot_by_class)
ggsave("plots/vk_boxplot_by_class_grouped_by_site.png", boxplot_by_class, width = 16, height = 12, dpi = 300)

# Option 1: Free y-axis scales (most likely what you want)
boxplot_combined_alt <- ggplot(top_vk_abundance_grouped,
                               aes(x = class, y = abundance, fill = cultivar)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~ site, nrow = 1, scales = "free_y") +  # Free y-axis scales
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(  x = " ",
       y = "Relative Abundance (%)",
       fill = "Cultivar"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  )

# Option 2: Add data points + better spacing
boxplot_combined_enhanced <- ggplot(top_vk_abundance_grouped,
                                    aes(x = class, y = abundance, fill = cultivar)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), 
               outlier.shape = NA) +  # Remove outlier points to avoid duplication
  geom_jitter(aes(color = cultivar), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
              alpha = 0.6, size = 1.5) +
  facet_wrap(~ site, nrow = 1, scales = "free_y") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  labs(x = " ",
       y = "Relative Abundance (%)",
       fill = "Cultivar"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  ) +
  guides(color = "none")  # Hide color legend since fill legend shows the same info

# Option 3: Vertical layout for better readability
boxplot_vertical <- ggplot(top_vk_abundance_grouped,
                           aes(x = class, y = abundance, fill = cultivar)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~ site, ncol = 1, scales = "free_y") +  # 2 columns instead of 1 row
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(x = "",
       y = "Relative Abundance (%)",
       fill = "Cultivar"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )


# Option 2 Enhanced: Add data points + better spacing + improved boxes
boxplot_combined_enhanced <- ggplot(top_vk_abundance_grouped,
                                    aes(x = class, y = abundance, fill = cultivar)) +
  geom_boxplot(alpha = 0.8, position = position_dodge(width = 0.8), 
               outlier.shape = NA, # Remove outlier points to avoid duplication
               size = 0.6, # Slightly thicker box lines
               fatten = 2) + # Thicker median line
  geom_jitter(aes(color = cultivar), 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15),
              alpha = 0.7, size = 1.2) +
  facet_wrap(~ site, ncol = 1, scales = "free_y") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_color_brewer(type = "qual", palette = "Set2") +
  labs( x = "",
       y = "Relative Abundance (%)",
       fill = "Cultivar"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 9, margin = margin(1, 0, 1, 0)), # Very small margin
   strip.background = element_rect(fill = "grey95", color = "grey70", size = 0.3), # Very light grey
    legend.position = "bottom",
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(10, 10, 5, 10), # Reduce bottom margin
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.minor = element_blank() # Remove minor grid lines
  ) +
  guides(color = "none") + # Hide color legend since fill legend shows the same info
  coord_cartesian(expand = FALSE) # Remove extra space around plot



# Print and save your preferred version
print(boxplot_combined_alt)  # Option 1: Free y-axis
ggsave("plots/vk_boxplot_combined_free_y.png", boxplot_combined_alt, width = 18, height = 6, dpi = 300)

print(boxplot_combined_enhanced)  # Option 2: With data points
ggsave("plots/vk_boxplot_enhanced.png", boxplot_combined_enhanced, width = 6, height = 18, dpi = 300)


# ===== 4. Visualize AImod, NOSC, and GFE across samples (TOP only) =====
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

# Prepare data for visualization (TOP samples only)
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

# Boxplot visualization for molecular properties by site-cultivar groups (TOP only)
properties_boxplot <- ggplot(properties_long,
                             aes(x = group, y = value, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_wrap(~property, scales = "free_y", ncol = 1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  theme_bw() +
  labs(
    title = "Distribution of Molecular Properties (TOP Samples Only)",
    subtitle = "By Site and Cultivar",
    x = "Site - Cultivar",
    y = "Value",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )
print(properties_boxplot)
ggsave("plots/Mean_molecular_properties_TOP_only_boxplot.png", properties_boxplot, width = 12, height = 10, dpi = 300)

# Properties comparison by site and cultivar separately
properties_by_site <- ggplot(properties_long,
                             aes(x = site, y = value, fill = cultivar)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.5, size = 1,
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1)) +
  facet_wrap(~property, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")) +
  theme_bw() +
  labs(x = "Site",
    y = "Value",
    fill = "Cultivar"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right"
  )

print(properties_by_site)
ggsave("plots/Mean_molecular_properties_by_site_TOP_only.png", properties_by_site, width = 12, height = 10, dpi = 300)

# Summary statistics table (TOP only)
properties_summary <- properties_long %>%
  group_by(group, property) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(property, group)
print(properties_summary)
write_csv(properties_summary, "processed_data/molecular_properties_summary_TOP_only.csv")

# ==== 5. Molecular formula analysis (TOP samples only) ====
# Analyze formulas present in different site-cultivar combinations
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
    group = interaction(site, cultivar, sep = " - "),
    sample_id = str_extract(sample_name, "\\d+_\\d+")
  )

# Count formulas by group
formula_count_by_group <- formula_group_analysis %>%
  group_by(group) %>%
  summarise(
    n_formulas = n_distinct(molecular_formula),
    n_samples = n_distinct(sample_name),
    .groups = "drop"
  )

# Analyze shared vs unique formulas between groups
formula_presence_matrix <- formula_group_analysis %>%
  select(molecular_formula, group) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = group, values_from = present, values_fill = 0)

# Calculate overlap between groups
group_names <- colnames(formula_presence_matrix)[-1]
overlap_analysis <- expand_grid(
  group1 = group_names,
  group2 = group_names
) %>%
  filter(group1 < group2) %>%
  rowwise() %>%
  mutate(
    shared = sum(formula_presence_matrix[[group1]] == 1 & formula_presence_matrix[[group2]] == 1),
    unique_group1 = sum(formula_presence_matrix[[group1]] == 1 & formula_presence_matrix[[group2]] == 0),
    unique_group2 = sum(formula_presence_matrix[[group1]] == 0 & formula_presence_matrix[[group2]] == 1),
    total_group1 = sum(formula_presence_matrix[[group1]] == 1),
    total_group2 = sum(formula_presence_matrix[[group2]] == 1),
    jaccard_similarity = shared / (total_group1 + total_group2 - shared)
  ) %>%
  ungroup()

print("Formula overlap analysis:")
print(overlap_analysis)
write_csv(overlap_analysis, "processed_data/formula_overlap_analysis_TOP_only.csv")

# Visualize formula counts by group
formula_count_plot <- ggplot(formula_count_by_group, aes(x = group, y = n_formulas, fill = group)) +
  geom_col(alpha = 0.8) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(
    title = "Number of Unique Molecular Formulas by Group (TOP Samples Only)",
    x = "Site - Cultivar",
    y = "Number of Formulas",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  ) +
  geom_text(aes(label = n_formulas), vjust = -0.5, size = 4)
print(formula_count_plot)
ggsave("plots/formula_count_by_group_TOP_only.png", formula_count_plot, width = 10, height = 6, dpi = 300)

