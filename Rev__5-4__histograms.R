#### Configuration ####
library(tidyverse)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
output_path_hist <- 'Results/correlation_histograms.svg'

dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')
corr_path <- c('Data/S01_Feb02_correlation_10s_swa.csv',
               'Data/S05_Jul11_correlation_10s_swa.csv',
               'Data/S05_Jul12_correlation_10s_swa.csv',
               'Data/S05_Jul13_correlation_10s_swa.csv')

cell_path <- 'Data/cell_types.csv'

micro_regions <- c('CLA', 'ACC', 'AMY')
micro_colors <- c('#E28DB8', '#A67A77', '#7BA387')

# responses <- c('Positive', 'None', 'Negative')
# response_colors <- c('#D46A6A', 'gray', '#6A89D4')

#### Munging ####

data <- data.frame()

# Load, merge, and concatenate data
for (i in 1:length(dates)) {
  
  corr <- read.csv(corr_path[i]) %>%
    select(unit_id, rho)
  
  corr$unit_id_full <- paste0(corr$unit_id, '_', dates[i]) 
  
  cell <- read.csv(cell_path) %>%
    select(unit_id, unit_region, response, cell_type) %>%
    rename(unit_id_full = unit_id)
  
  corr <- merge(corr, cell, by = 'unit_id_full') %>%
    mutate(unit_id = unit_id_full) %>%
    select(-unit_id_full)
  
  data <- rbind(data, corr)
}

# Format columns
data$unit_region <- factor(data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
data$response <- factor(data$response, levels = c('Positive', 'None', 'Negative'))
data$cell_type <- factor(data$cell_type, levels = c('pyramidal', 'unknown', 'interneuron'))
data <- data %>%
  mutate(cell_type = fct_recode(cell_type,
                                Pyramidal = "pyramidal",
                                Interneuron = "interneuron",
                                Unknown = "unknown")) 

#### Plotting ####
p <- ggplot(data, aes(x = rho, fill = unit_region)) + 
  geom_histogram(alpha = 0.3) +
  geom_histogram(aes(color = unit_region), fill = NA) +  
  scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  scale_fill_manual(values = setNames(micro_colors, micro_regions)) +
  facet_grid(response ~ unit_region) +
  labs(x = "Spearman's Rho", y = "Unit-channel pairs (n)") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    panel.border = element_rect(colour = "black", fill=NA, size=0.75),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_path_hist, plot = p, width = 6, height = 6, units = "in")