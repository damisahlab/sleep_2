#### Configuration ----
library(tidyverse)
library(ggh4x)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep')

#### Parameters ----
spike_path <- 'Cache/S02_spikes.csv'
sws_corr_path <- 'Cache/S02_sws_corr.csv'
plot_path <- 'Results/S02_unit_raster.svg'

#### Parameters ----
spike_path <- 'Cache/Subject01/S01_spikes.csv'
sws_corr_path <- 'Cache/Subject01/S01_sws_corr.csv'
plot_path <- 'Results/S01_unit_raster.svg'

#### Spike Munging ----
## Load and Format Spikes
spikes <- 
  read.csv(spike_path) %>%
  mutate(time = round(seconds/60, digits = 1)) %>%
  select(unit_id, unit_laterality, unit_region, sign, time) %>%
  distinct() %>%
  group_by(unit_id) %>%
  mutate(sort_value = sum(n()))

spikes$side <- ifelse(spikes$unit_laterality == 'left', 'Left', 'Right')

## Merge with SW's to get distance correlation
sws_corr <-
  read.csv(sws_corr_path) %>%
  select(unit_id, dcor) %>%
  group_by(unit_id) %>%
  filter(dcor == max(dcor))

data <- merge(spikes, sws_corr, by = 'unit_id'); rm(sws_corr)

#### Plotting ----
ggplot(data, aes(x = time, y = reorder(unit_id, sort_value))) + 
  geom_raster(aes(fill = dcor)) + 
  #facet_wrap(side + unit_region ~ ., scales = 'free') + 
  facet_nested(side ~ ., 
               nest_line = element_line(linetype = 1),
               scale = 'free', 
               space = 'free_y') +
  scale_fill_continuous(#limit = c(0.2, 0.4),
    #oob = scales::squish,
    low = 'gray90', 
    high = 'red') +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Time (min)',
       y = 'Claustrum Unit (ID)',
       fill = 'Distance\nCorrelation',
       caption = 'Units sorted by spike count; color scale is cropped') +
  theme_minimal() + 
  theme(#legend.title = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 10),
    #axis.text.y = element_blank(),
    axis.ticks.x = element_line(size = 0.4),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = 'black', 
                                fill = NA, 
                                size = 0.5))

ggsave(file = plot_path, 
       width = 10, height = 6, dpi = 300)