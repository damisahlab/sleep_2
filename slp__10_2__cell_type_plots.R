#### Configuration ----
library(tidyverse)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ----
metrics_path <- 'Data/cell_type_metrics.csv'

recording_dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')
response_path <- c('Data/S01_Feb02_unit_response.csv',
                   'Data/S05_Jul11_unit_response.csv',
                   'Data/S05_Jul12_unit_response.csv',
                   'Data/S05_Jul13_unit_response.csv')
                   
csv_output <- 'Data/cell_types.csv'
plot_output_1 <- 'Results/celltype_plot_1.svg'
plot_output_2 <- 'Results/celltype_plot_2.svg'
plot_output_3 <- 'Results/celltype_plot_3.svg'

micro_regions <- c('CLA', 'ACC', 'AMY')
micro_colors <- c('#E28DB8', '#A67A77', '#7BA387')

#### Load Data ----
response <- data.frame()
for(i in 1:length(recording_dates)){
  temp <- read.csv(response_path[i]) %>%
    mutate(unit_id_old = unit) %>%
    mutate(unit_id = paste(unit, recording_dates[i], sep = '_'),
           unit_laterality = str_split(unit_roi, " ", simplify = TRUE)[, 1],
           unit_region = str_split(unit_roi, " ", simplify = TRUE)[, 2]) %>%
    mutate(unit_region = factor(unit_region, levels = micro_regions),
           response_type = factor(response_type, levels = c('Positive', 'None', 'Negative'))) %>%
    select(unit_id, unit_id_old, unit_region, unit_laterality, response_type)
  
  temp$recording <- recording_dates[i]
  
  response <- rbind(response, temp)
}

df <- read.csv(metrics_path) %>%
  select(unit_id, firing_rate, trough_to_peak, fwhm, 
         burst_index, log_firing_rate) %>%
  merge(response, by = 'unit_id')

#### Classify cell type ----
df <- df %>%
  mutate(
    interneuron = as.integer((firing_rate > 10) + 
                               ifelse(is.na(trough_to_peak), 0, (trough_to_peak < 0.5)) + 
                               (burst_index < 0.2)),
    pyramidal = ifelse(is.na(trough_to_peak), 2 - (firing_rate > 15) - (burst_index < 0.2), 3 - interneuron)
  ) %>%
  mutate(
    cell_type = case_when(
      interneuron > pyramidal ~ "interneuron",
      interneuron < pyramidal ~ "pyramidal",
      TRUE ~ "unknown"
    )
  )

write.csv(df, csv_output, row.names = FALSE)

#### Summary ----

cross_tab <- df %>%
  group_by(unit_region, cell_type) %>%
  summarise(count = n()) %>%
  spread(key = cell_type, value = count, fill = 0)

cross_tab

#### Plotting ----
plot_1 <- 
  ggplot(df, aes(x = log_firing_rate, y = burst_index, 
                 color = unit_region, shape = cell_type)) +
  geom_point(size = 5, alpha = 0.75) + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(response_type ~ .) + 
  labs(x = 'Firing Rate (Log2 Hz)', y = 'Burst Index (au)') + 
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plot_1
ggsave(filename = plot_output_1, plot = plot_1, width = 9, height = 3)

plot_2 <- 
  # Remember that TTP latency is only calculated for negative-spiking units!
  ggplot(df, aes(x = log_firing_rate, y = trough_to_peak, 
                 color = unit_region, shape = cell_type)) +
  geom_point(size = 5, alpha = 0.75) + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(response_type ~ .) + 
  labs(x = 'Firing Rate (Log2 Hz)', y = 'Trough-to-Peak Latency (ms)') + 
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plot_2
ggsave(filename = plot_output_2, plot = plot_2, width = 9, height = 3)

plot_3 <- 
  ggplot(df, aes(x = log_firing_rate, y = fwhm, 
                 color = unit_region, shape = cell_type)) +
  geom_point(size = 5, alpha = 0.75) + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(response_type ~ .) + 
  labs(x = 'Firing Rate (Log2 Hz)', y = 'Full Width at Half Maximum (ms)') + 
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plot_3
ggsave(filename = plot_output_3, plot = plot_3, width = 9, height = 3)