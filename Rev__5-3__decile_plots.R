#### Configuration ####
library(tidyverse)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
output_path <- 'Results/response_decile_plot.svg'
  
dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')

spike_path <- c('Data/S01_Feb02_epochs_10s_spikes.csv',
                'Data/S05_Jul11_epochs_10s_spikes.csv',
                'Data/S05_Jul12_epochs_10s_spikes.csv',
                'Data/S05_Jul13_epochs_10s_spikes.csv')

swa_path <- c('Data/S01_Feb02_epochs_10s_swa.csv',
              'Data/S05_Jul11_epochs_10s_swa.csv',
              'Data/S05_Jul12_epochs_10s_swa.csv',
              'Data/S05_Jul13_epochs_10s_swa.csv')

channel_path <- c('Data/S01_electrodes.csv',
                  'Data/S05_electrodes.csv',
                  'Data/S05_electrodes.csv',
                  'Data/S05_electrodes.csv')

cell_path <- 'Data/cell_types.csv'

micro_regions <- c('CLA', 'ACC', 'AMY')
micro_colors <- c('#E28DB8', '#A67A77', '#7BA387')

#### Munging ####
# Load spike data
data <- data.frame() 

for (i in 1:length(dates)) {
  
  # Load spike data, add full unit_id (recording suffix)
  spike <- read.csv(spike_path[i]) %>%
    group_by(unit_id) %>%
    mutate(fr = scale(fr, center = TRUE, scale = TRUE)) %>%
    as.vector()
  
  spike$unit_id_full <- paste0(spike$unit_id, '_', dates[i])
  spike$recording <- dates[i]
  
  # Load cell type data
  cell <- read.csv(cell_path) %>%
    select(unit_id, response, cell_type) %>%
    rename(unit_id_full = unit_id)
  
  # Add cell type meta-data to spike data
  spike <- merge(spike, cell, by = 'unit_id_full') %>%
    mutate(unit_id = unit_id_full) %>%
    select(-unit_id_full)
  
  # Load macro channel meta-data
  chan <- read.csv(channel_path[i]) %>% 
    select(elec_label, hemisphere, lobe_1) %>%
    rename(channel = elec_label,
           chan_laterality = hemisphere,
           chan_lobe = lobe_1)
  
  # Load SWA data, merge with channel meta-data,
  # and then calculate global SWA (average first
  # by laterality/lobe, then globally)
  swa <- read.csv(swa_path[i]) %>%
    select(epoch, channel, zlog_power) %>%
    merge(chan, by = 'channel') %>%
    group_by(epoch, chan_laterality, chan_lobe) %>%
    summarise(power = mean(zlog_power, na.rm = TRUE)) %>%
    ungroup() %>% 
    group_by(epoch) %>%
    summarise(power = mean(power))
  
  loop_data <- merge(spike, swa, by = 'epoch') %>% 
    mutate(decile = ntile(power, 10))
  
  data <- rbind(data, loop_data)
}

# Format columns
data$unit_region <- factor(data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
data$response <- factor(data$response, levels = c('Positive', 'None', 'Negative'))
data$cell_type <- factor(data$cell_type, levels = c('pyramidal', 'unknown', 'interneuron'))

data <- data %>%
  mutate(cell_type = fct_recode(cell_type,
                                Pyramidal = "pyramidal",
                                Interneuron = "interneuron",
                                Unknown = "unknown")) %>%
  select(recording, unit_region, response, cell_type, decile, power, fr)

# Summarize by 'decile', 'unit_region', and 'response'
sum_data <- data %>%
  group_by(recording, decile, unit_region, response) %>%
  summarize(
    power = mean(power, na.rm = TRUE),
    se = sd(fr, na.rm = TRUE) / sqrt(n()),
    fr = mean(fr, na.rm = TRUE),
    .groups = 'drop'
  )

#### Statistics ####
# Create an empty dataframe to store R-squared values
rs_data <- data.frame(
  unit_region = character(),
  response = character(),
  R_squared = numeric()
)

# Iterate over each combination of unit_region and response
for (region in unique(sum_data$unit_region)) {
  for (resp in unique(sum_data$response)) {
    # Filter data for the current combination of unit_region and response
    filtered_data <- sum_data %>% 
      filter(unit_region == region, response == resp)
    
    # Fit a model
    model <- lm(fr ~ power, data = filtered_data)
    r_squared <- summary(model)$r.squared

    # Append the results to the rs_data dataframe
    rs_temp <- data.frame(unit_region = region, response = resp, R_squared = r_squared)
    rs_data <- rbind(rs_data, rs_temp)
  }
}

rs_data$unit_region <- factor(rs_data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
rs_data$response <- factor(rs_data$response, levels = c('Positive', 'None', 'Negative'))
print(rs_data)

#### Plotting ####
p <- ggplot(sum_data, aes(x = power, y = fr, color = unit_region)) + 
  geom_point(size = 1.5) +
  geom_linerange(aes(ymin = fr - se, ymax = fr + se)) + 
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5) +
  geom_text(data = rs_data, 
            aes(label = sprintf("R² = %.2f", R_squared), x = Inf, y = Inf), 
            hjust = 3.6, vjust = 1.5, 
            check_overlap = TRUE, color = "black") +
  scale_color_manual(values = setNames(micro_colors, micro_regions)) +
  facet_grid(unit_region ~ response) +
  labs(x = "Global SWA (z-log)", y = "zFR") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p

ggsave(output_path, plot = p, width = 8, height = 6, units = "in")
