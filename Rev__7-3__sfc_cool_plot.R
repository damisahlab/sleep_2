#### Configuration ####
library(tidyverse)
library(circular)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
selected_subset <- 'DREM'
p_threshold <- 0.05

recording_dates <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')

input_path <- c('Cache/S01_Feb02_permuted_rayleigh.csv',
                'Cache/S05_Jul11_permuted_rayleigh.csv',
                'Cache/S05_Jul12_permuted_rayleigh.csv',
                'Cache/S05_Jul13_permuted_rayleigh.csv')

channel_path <- c('Data/S01_electrodes.csv',
                  'Data/S05_electrodes.csv',
                  'Data/S05_electrodes.csv',
                  'Data/S05_electrodes.csv')

response_path <- 'Data/cell_types.csv'

output_path_1 <- 'Results/sfc_nspikes.svg'
output_path_2 <- 'Results/sfc_amy_norm.svg'
output_path_3 <- 'Results/sfc_coolplot.svg'

unit_color_vector <- c('CLA' = '#e8a2c5', 
                       'ACC' = '#805553', 
                       'AMY' = '#7da58a')

# NOTE: You may need to expand this for the other recordings
chan_region_order <- c('SFG', 'MFG', 'IFG', 'OFG', 'MCC', 'PCC', 'PM', 'S', 
                       'SPL', 'PNS', 'STG', 'MTG', 'TPL', 'AIN', 'PIN', 
                       'AMY', 'HPC')

#### Munging ####
unit_meta <- 
  read.csv(response_path) %>%
  select(unit_id, unit_region, unit_laterality, response_type) %>%
  rename(unit_response = response_type)

data <- data.frame()

for(i in 1:length(recording_dates)){
  chan_meta <- 
    read.csv(channel_path[i]) %>% 
    select(elec_label, hemisphere, roi_3) %>%
    rename(channel = elec_label,
           chan_laterality = hemisphere,
           chan_region = roi_3)
  
  temp <- 
    read.csv(input_path[i]) %>%
    mutate(unit_id = paste(unit_id, recording_dates[i], sep = '_'),
           recording_date = recording_dates[i]) %>%
    merge(chan_meta, by = 'channel') %>%
    merge(unit_meta, by = 'unit_id') %>%
    filter(p_corrected < p_threshold)
  
  data <- rbind(data, temp)
}

data$unit_region <- factor(data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
data$chan_region <- factor(data$chan_region, levels = rev(chan_region_order))
data$mean_phase <- circular(data$mean_phase, units = 'radians')

# Subset to DREM and WREM pairs
drem <- 
  data %>% 
  filter(subset == 'DREM') %>%
  #filter(mean_rl > 0.2) %>%
  filter(n_spikes > 5000) %>%
  select(unit_id, channel, subset, mean_phase, mean_rl,
         chan_region, unit_region, n_spikes)

wrem <- data %>% filter(subset == 'WREM')

# Remove rows where the unit_id-channel pair is present in WREM
excluded_pairs <- unique(with(wrem, paste(unit_id, channel, sep = "-")))
drem$unit_channel_pair <- with(drem, paste(unit_id, channel, sep = "-"))
drem <- subset(drem, !unit_channel_pair %in% excluded_pairs)

#### Plot 1 ####
ggplot(data, aes(x = log10(n_spikes), fill = unit_region)) + 
  geom_histogram(alpha = 0.75) +
  geom_vline(xintercept = 3, color = "black", linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept = 4, color = "black", linetype = 'dashed', size = 0.5) +
  scale_color_manual(values = unit_color_vector) +
  scale_fill_manual(values = unit_color_vector) +
  scale_y_continuous(expand = c(0, 0)) + 
  facet_grid(subset ~ unit_region) + # , scales = 'free_y'
  labs(title = 'Lines set to n = 1000, 10000', x = "Log10 Number of Spikes", y = "Pairs (n)") +
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

ggsave(file = output_path_1, width = 8, height = 6, dpi = 300)

#### Plot 2 ####
mean_amy_rl <- mean(drem[(drem$chan_region == 'AMY') & 
                           (drem$unit_region == 'AMY') & 
                           (drem$subset == 'DREM'), 'mean_rl'], na.rm = TRUE)


drem$norm_rl <- drem$mean_rl / mean_amy_rl

ggplot(drem, aes(x = mean_phase, y = norm_rl, color = unit_region)) +
  geom_point(alpha = 0.6) + 
  geom_hline(yintercept = 1, color = "black", linetype = 'solid', size = 0.5) +
  scale_x_continuous(breaks = c(seq(-pi, pi, pi/2)), 
                     labels = c('-π', '-π/2', '0', 'π/2', 'π'),
                     limits = c(-pi, pi),
                     expand = c(0.01, 0)) +
  scale_color_manual(values = unit_color_vector) +
  scale_fill_manual(values = unit_color_vector) +
  facet_wrap(unit_region ~ .) +
  labs(title = "Significant pairs in DREM (but not in WREM, n > 5000), normalized to its subset of AMY-AMY pairs",
       x = "Mean Phase", y = "Normalized MRL") + 
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.spacing = unit(0.2, 'lines'), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(colour = 'black', 
                                    fill = NA, 
                                    size = 1))

ggsave(file = output_path_2, width = 6, height = 3, dpi = 300)

#### Plot 3 ####
roi_angles <-
  drem %>%
  group_by(chan_region, unit_region) %>%
  summarize(mean_angle = mean.circular(mean_phase)[[1]]) %>%
  ungroup()

ggplot() +
  geom_point(data = drem, 
             aes(x = mean_phase, y = chan_region, color = unit_region),
             size = 3, alpha = 0.5) + 
  geom_point(data = roi_angles, 
             aes(x = mean_angle, y = chan_region, fill = unit_region), 
             color = 'black', shape = 23, size = 4) +
  scale_color_manual(values = unit_color_vector) +
  scale_fill_manual(values = unit_color_vector) +
  scale_x_continuous(breaks = c(seq(-pi, pi, pi/2)), 
                     labels = c('-π', '-π/2', '0', 'π/2', 'π'),
                     limits = c(-pi, pi),
                     expand = c(0.01, 0)) +
  labs(title = 'Pairs significant in DREM (but not in WREM) with min_spikes > 5000',
       x = "Mean Phase Angle (radians)", y = "Channel") + 
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.spacing = unit(0.2, 'lines'), 
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_rect(colour = 'black', 
                                    fill = NA, 
                                    size = 1))

ggsave(file = output_path_3, width = 6, height = 4, dpi = 300)
