#### Configuration ####
library(tidyverse)
library(gridExtra)

setwd('Z:/Layton/Sleep_051324')

#### Parameters ####
input_path <- 'Cache/waveform_stability.csv'
output_path <- 'Results/waveform_stability.svg'

micro_regions <- c('CLA', 'ACC', 'AMY')
micro_colors <- c('#5C92C0', '#80BA55', '#E69978')

#### Load Data ####
df <- read.csv(input_path)

#### Plotting ####
plot1 <- 
  ggplot(data = df, aes(x = fwhm_cv, y = max_amplitude_cv, color = unit_region)) + 
  geom_point(size = 2.5) + 
  labs(x = 'FWHM (CoV)', y = 'Amplitude (CoV)') + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

plot2 <-
  ggplot(data = df, aes(x = auc_cv, y = max_amplitude_cv, color = unit_region)) + 
  geom_point(size = 2.5) + 
  labs(x = 'AUC (CoV)', y = 'Amplitude (CoV)') + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

plot3 <- 
  ggplot(data = df, aes(x = euclidean_distance_cv, y = l_ratio, color = unit_region)) + 
  geom_point(size = 2.5) + 
  labs(x = 'Euclidean Distance from Mean (CoV)', y = 'L-Ratio') + 
  scale_color_manual(values = setNames(micro_colors, micro_regions)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', fill = NA, size = 1))

# Show plot
grid.arrange(plot1, plot2, plot3, ncol = 3)

# Save plot
grid_layout <- arrangeGrob(plot1, plot2, plot3, ncol = 3)
ggsave(output_path, grid_layout, width = 25, height = 7.5, units = "cm")
