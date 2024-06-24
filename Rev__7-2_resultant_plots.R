#### Configuration ####
library(tidyverse)
library(hexbin)
library(viridis)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
recording_dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')

input_path <- c('Cache/S01_Feb02_permuted_rayleigh.csv',
                'Cache/S05_Jul11_permuted_rayleigh.csv',
                'Cache/S05_Jul12_permuted_rayleigh.csv',
                'Cache/S05_Jul13_permuted_rayleigh.csv')

meta_path <- 'Data/cell_types.csv'

output_hist <- 'Results/sfc_histogram.svg'
output_heat <- 'Results/sfc_heatmap.svg'

micro_regions <- c('CLA', 'AMY', 'ACC')
micro_colors <- c('#E28DB8', '#7BA387', '#A67A77')

#### Munging ####
unit_meta <- read.csv(meta_path) %>%
  select(unit_id, unit_region, unit_roi, response) %>%
  rename(unit_laterality = unit_roi,
         unit_response = response)

data <- data.frame()

for (i in 1:length(recording_dates)) {
  
  temp <- read.csv(input_path[i]) %>%
    mutate(unit_id = paste(unit_id, recording_dates[i], sep = '_'),
           recording_date = recording_dates[i]) %>%
    merge(unit_meta, by = 'unit_id')
  
  data <- rbind(data, temp)
  
}

data$subset <- factor(data$subset, levels = c('DREM', 'NREM', 'WREM'))
data$unit_region <- factor(data$unit_region, levels = c('CLA', 'AMY', 'ACC'))

data_hist <- data %>% filter(p_corrected < 0.05)
data_heat <- data

#### Plotting ####
ggplot(data_hist, aes(x = mean_rl, fill = unit_region)) + 
  geom_histogram(alpha = 0.75) +
  scale_discrete_manual(values = micro_colors, aesthetics = c('fill')) + 
  scale_y_continuous(expand = c(0, 0)) + 
  facet_grid(subset ~ unit_region) + # , scales = 'free_y'
  labs(x = "Mean Resultant Length", y = "Pairs (n)") +
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

ggsave(file = output_hist, width = 6, height = 5, dpi = 300)

ggplot(data_heat, aes(x = mean_phase, y = mean_rl)) + 
  geom_hex() + 
  scale_fill_viridis(option = 'plasma') +
  #scale_x_continuous(expand = c(0, 0)) + 
  scale_x_continuous(breaks = c(seq(-pi, pi, pi/2)), 
                     labels = c('-π', '-π/2', '0', 'π/2', 'π'),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  facet_grid(subset ~ unit_region) + 
  labs(x = "Mean Phase", y = "Mean Resultant Length", fill = "Count") +
  theme_minimal() +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

ggsave(file = output_heat, width = 6, height = 5, dpi = 300)
