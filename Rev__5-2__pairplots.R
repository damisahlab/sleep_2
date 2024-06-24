#### Configuration ####
library(tidyverse)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
output_path_response <- 'Results/boxplots_response.svg'
output_path_type <- 'Results/boxplots_celltype.svg'

dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')

spikes_path <- c('Data/S01_Feb02_epochs_30s_spikes.csv',
                 'Data/S05_Jul11_epochs_30s_spikes.csv',
                 'Data/S05_Jul12_epochs_30s_spikes.csv',
                 'Data/S05_Jul13_epochs_30s_spikes.csv')

hypno_path <- c('Data/S01_Feb02_epochs_30s_hypno.csv',
                'Data/S05_Jul11_epochs_30s_hypno.csv',
                'Data/S05_Jul12_epochs_30s_hypno.csv',
                'Data/S05_Jul13_epochs_30s_hypno.csv')

cell_path <- 'Data/cell_types.csv'

micro_regions <- c('CLA', 'ACC', 'AMY')
micro_colors <- c('#E28DB8', '#A67A77', '#7BA387')

#### Munging ####
data <- data.frame()

# Load, merge, and concatenate data
for (i in 1:length(dates)) {
  
  hypno <- read.csv(hypno_path[i]) %>%
    mutate(stage = case_when(
      stage == 2 ~ 'NREM',
      stage == 3 ~ 'NREM',
      TRUE ~ 'WREM'))
  
  spikes <- read.csv(spikes_path[i]) %>%
    filter(unit_region %in% micro_regions) %>%
    mutate(unit_region = factor(unit_region, levels = micro_regions)) %>%
    select(epoch, unit_id, unit_laterality, unit_region, fr)
  
  loop_data <- merge(hypno, spikes, by = 'epoch')
  
  loop_data$unit_id_full <- paste0(loop_data$unit_id, '_', dates[i]) 
  
  cell <- read.csv(cell_path) %>%
    select(unit_id, response, cell_type) %>%
    rename(unit_id_full = unit_id)
  
  loop_data <- merge(loop_data, cell, by = 'unit_id_full') %>%
    mutate(unit_id = unit_id_full) %>%
    select(-unit_id_full)
  
  data <- rbind(data, loop_data)
}

# Z-score the firing rate
data <- 
  data %>%
  mutate(stage = factor(stage, levels = c('WREM', 'NREM'))) %>%
  group_by(unit_id) %>%
  mutate(zfr = (fr - mean(fr)) / sd(fr)) %>%
  ungroup()

#### Subsets ####
data$unit_region <- factor(data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
data$response <- factor(data$response, levels = c('Positive', 'None', 'Negative'))
data$cell_type <- factor(data$cell_type, levels = c('pyramidal', 'unknown', 'interneuron'))
data <- 
  data %>%
  mutate(cell_type = fct_recode(cell_type,
                                Pyramidal = "pyramidal",
                                Interneuron = "interneuron",
                                Unknown = "unknown")) 
  
# Group means of response subset
datar_sum <- 
  data %>%
  group_by(unit_id, unit_region, response, stage) %>%
  summarize(zfr = mean(zfr)) %>%
  ungroup()

# Group means of cell type subset
datac_sum <- 
  data %>%
  group_by(unit_id, unit_region, cell_type, stage) %>%
  summarize(zfr = mean(zfr)) %>%
  ungroup()

#### Statistics ####
# Initialize an empty dataframe for results
results_r <- data.frame(unit_region = character(),
                      response = character(),
                      p_value = numeric())

# Loop through each combination of unit_region and response
for (region in unique(datar_sum$unit_region)) {
  for (resp in unique(datar_sum$response)) {
    # Subset data for current combination
    subset_data <- datar_sum %>%
      filter(unit_region == region, response == resp)
    
    # Check if there are enough data points for a test
    if(nrow(subset_data) >= 2) {
      # Reshape data to wide format for paired test
      wide_data <- subset_data %>%
        pivot_wider(names_from = stage, values_from = zfr)
      
      # Perform Wilcoxon Signed-Rank Test and extract p-value
      test_result <- wilcox.test(wide_data$NREM, wide_data$WREM, paired = TRUE)
      p_value <- test_result$p.value
    } else {
      # Not enough data points for a test
      p_value <- NA
    }
    
    # Store the result
    results_r <- rbind(results_r, data.frame(unit_region = region,
                                             response = resp,
                                             p_value = p_value))
  }
}

# Initialize an empty dataframe for results
results_c <- data.frame(unit_region = character(),
                      cell_type = character(),
                      p_value = numeric())

# Loop through each combination of unit_region and cell_type
for (region in unique(datac_sum$unit_region)) {
  for (cell in unique(datac_sum$cell_type)) {
    # Subset data for current combination
    subset_data <- datac_sum %>%
      filter(unit_region == region, cell_type == cell)
    
    # Check if there are enough data points for a test
    if(nrow(subset_data) >= 2) {
      # Reshape data to wide format for paired test
      wide_data <- subset_data %>%
        pivot_wider(names_from = stage, values_from = zfr)
      
      # Perform Wilcoxon Signed-Rank Test and extract p-value
      test_result <- wilcox.test(wide_data$NREM, wide_data$WREM, paired = TRUE)
      p_value <- test_result$p.value
    } else {
      # Not enough data points for a test
      p_value <- NA
    }
    
    # Store the result
    results_c <- rbind(results_c, data.frame(unit_region = region,
                                             cell_type = cell,
                                             p_value = p_value))
  }
}

# Add a column for FDR-corrected p-values using Benjamini-Hochberg correction
results_r$fdr_p_value <- p.adjust(results_r$p_value, method = "BH")
datar_sum <- merge(datar_sum, results_r, by = c("unit_region", "response"))

# Add a column for FDR-corrected p-values using Benjamini-Hochberg correction
results_c$fdr_p_value <- p.adjust(results_c$p_value, method = "BH")
datac_sum <- merge(datac_sum, results_c, by = c("unit_region", "cell_type"))

#### Plotting ####
# Plot 1
datar_sum$label <- ifelse(datar_sum$fdr_p_value < 0.01,
                          sprintf("p[FDR] == %.1e", datar_sum$fdr_p_value),
                          sprintf("p[FDR] == %.2f", datar_sum$fdr_p_value))

p1 <- ggplot(datar_sum, aes(x = stage, y = zfr, color = unit_region, group = unit_id)) + 
  geom_line(alpha = 0.15) +
  geom_point(alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black", 
               aes(group = interaction(unit_region, response, stage)), fatten = 0.5) +
  geom_text(aes(label = label, x = 1.5, y = -0.8), parse = TRUE, size = 3,
            hjust = 0, vjust = 0, check_overlap = TRUE, color = "black") +
  labs(x = "Sleep Stage", y = "zFR") +
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_grid(unit_region ~ response) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill=NA, size=0.75),  
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    panel.grid.minor.y = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1))

p1

ggsave(output_path_response, plot = p1, width = 6, height = 6, units = "in")

# Plot 2
datac_sum$label <- ifelse(datac_sum$fdr_p_value < 0.01,
                          sprintf("p[FDR] == %.1e", datac_sum$fdr_p_value),
                          sprintf("p[FDR] == %.2f", datac_sum$fdr_p_value))

p2 <- ggplot(datac_sum, aes(x = stage, y = zfr, color = unit_region, group = unit_id)) + 
  geom_line(alpha = 0.15) +
  geom_point(alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black", 
               aes(group = interaction(unit_region, cell_type, stage)), fatten = 0.5) +
  geom_text(aes(label = label, x = 1.5, y = -0.8), parse = TRUE, size = 3,
            hjust = 0, vjust = 0, check_overlap = TRUE, color = "black") +
  labs(x = "Sleep Stage", y = "zFR") +
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_grid(unit_region ~ cell_type) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill=NA, size=0.75),  
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    panel.grid.minor.y = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1))

p2

ggsave(output_path_type, plot = p2, width = 6, height = 6, units = "in")

