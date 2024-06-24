# Load necessary libraries
library(ggplot2)
library(dplyr)
library(purrr)
library(ggbeeswarm)

# Parameters
setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

recording_lengths <- c(2, 9.68, 10.55, 10.4)
input_paths <- c('Data/S01_Feb02_spike_times.csv',
                 'Data/S05_Jul11_spike_times.csv',
                 'Data/S05_Jul13_spike_times.csv',
                 'Data/S05_Jul13_spike_times.csv')

output_path_cdf <- 'Results/spike_stability_cdf.svg'
output_path_swarm <- 'Results/spike_stability_swarm.svg'

# Define colors for unit regions
unit_region_colors <- c("CLA" = '#E28DB8',  # Pink
                        "AMY" = '#7BA387',  # Green
                        "ACC" = '#A67A77')  # Brown

# Import and bind the datasets
df_list <- map2(input_paths, recording_lengths, ~{
  df <- read.csv(.x)[, c("unit_laterality", "unit_region", "seconds", "unit_id")]
  df$recording_length <- .y
  return(df)
})

df <- bind_rows(df_list)

# Normalize seconds to percentage of total time
df <- df %>%
  mutate(total_time_sec = recording_length * 60 * 60,
         time_percent = seconds / total_time_sec * 100)

# Create 1% time bins and calculate cumulative spikes within each bin
df_binned <- df %>%
  mutate(time_bin = floor(time_percent)) %>%
  group_by(unit_region, unit_id, time_bin) %>%
  summarise(cumulative_spike_count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(unit_region, unit_id) %>%
  mutate(cumulative_spike_count = cumsum(cumulative_spike_count)) %>%
  ungroup()

# Normalize cumulative spike count to percentage
df_binned <- df_binned %>%
  group_by(unit_region, unit_id) %>%
  mutate(spike_percent = cumulative_spike_count / max(cumulative_spike_count) * 100) %>%
  ungroup()

# Plot the cumulative spike distribution function colored by unit_region
p <- ggplot(df_binned, aes(x=time_bin, y=spike_percent, group=unit_id, color=unit_region)) +
  geom_line(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="black") +
  scale_color_manual(values = unit_region_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Recording Time (%)', y = 'Spike Count (%)') +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Export the CDF plot as an SVG file
ggsave(output_path_cdf, plot = p, width = 8, height = 6, device = "svg")

# Calculate RMSE for each unit_id
rmse_df <- df_binned %>%
  group_by(unit_region, unit_id) %>%
  mutate(error = spike_percent - time_bin) %>%
  summarise(RMSE = sqrt(mean(error^2)), unit_region = first(unit_region)) %>%
  ungroup()

# Make a quasirandom (similar to swarm) plot of RMSE values for each unit_region
p_swarm <- ggplot(rmse_df, aes(x = unit_region, y = RMSE, color = unit_region, fill = unit_region)) +
  geom_quasirandom(shape = 21, size = 3, stroke = 1, alpha = 0.3) +  # Adjusted stroke for thicker outline
  guides(color = guide_legend(override.aes = list(alpha = 1))) +  # Override alpha for legend color
  scale_color_manual(values = unit_region_colors) +
  scale_fill_manual(values = unit_region_colors, guide = "none") +  # Match fill colors to color and remove fill legend
  labs(x = "Unit Region", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = 'none',
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_swarm

# Export the histogram plot as an SVG
ggsave(output_path_swarm, plot = p_swarm, width = 8, height = 6, device = "svg")
