#### Configuration ####
library(tidyverse)
library(ggsci)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
window_min <- -50
bin_width <- 5

input_path_tfce <- 'Cache/All_xcorr_tfce_b5_w10_n1000_nrem_allpairs.csv'

recordings <- c('Feb02', 'Jul11', 'Jul12', 'Jul13')
input_path_ogram <- c('Cache/S01_Feb02_ocorrelogram_b5_w10_n1000_nrem_allpairs.csv',
                      'Cache/S05_Jul11_ocorrelogram_b5_w10_n1000_nrem_allpairs.csv',
                      'Cache/S05_Jul12_ocorrelogram_b5_w10_n1000_nrem_allpairs.csv',
                      'Cache/S05_Jul13_ocorrelogram_b5_w10_n1000_nrem_allpairs.csv')

meta_path <- 'Data/cell_types.csv'
output_path <- 'Results/xcorr_ograms_allnights_allpairs_p05_s2000.svg'

#### Combine Original Correlograms ####
data_frames <- list()

for (i in seq_along(recordings)) {
  input_path <- paste(input_path_ogram[i], sep = "")
  data <- read.csv(input_path)
  data$unit_1 <- paste(data$unit_1, recordings[i], sep = "_")
  data$unit_2 <- paste(data$unit_2, recordings[i], sep = "_")
  data_frames[[i]] <- data 
}

data <- do.call(rbind, data_frames)

#### Munging ####

# Load data
tfce <- read.csv(input_path_tfce)
meta <- read.csv(meta_path)
meta <- meta[c("unit_id", "unit_region")]

# Merge with unit pair meta-data
data <- merge(data, tfce, by = c('unit_1', 'unit_2'))

# ... for unit_1
data <- merge(data, meta, by.x = "unit_1", by.y = "unit_id", all.x = TRUE)
names(data)[names(data) == "unit_region"] <- "unit_region_1"

# ... for unit_2
data <- merge(data, meta, by.x = "unit_2", by.y = "unit_id", all.x = TRUE)
names(data)[names(data) == "unit_region"] <- "unit_region_2"

# Create a spike count column
lag_columns <- grep("lag_\\d+", names(data), value = TRUE)
data$spikes <- rowSums(data[, lag_columns], na.rm = TRUE)

# Keep cross-correlograms with pre-FDR p-value < 0.05 & spikes > 5000
data <- data[data$perm_p_val < 0.05,]
data <- data[data$spikes > 2000,]

#### Cast to Long ####

# Cast the data to long format and extract bin numbers
data_long <- data %>%
  pivot_longer(cols = starts_with("lag_"), 
               names_to = "Bin", 
               names_prefix = "lag_",  # Remove 'lag_' prefix
               values_to = "Count") %>%
  mutate(Bin = as.integer(Bin))  # Convert bin numbers to integer

# Convert Bin numbers to lag in milliseconds
data_long$lag <- window_min + (data_long$Bin - 1) * bin_width

# Create a unique identifier for each unit pair and unit region pair
data_long$unit_pair <- paste(data_long$unit_1, data_long$unit_2, sep = " - ")
data_long$unit_region_pair <- paste(data_long$unit_region_1, data_long$unit_region_2, sep = " - ")

#### Plotting ####
p <- ggplot(data_long, aes(x = lag, y = Count)) +
  geom_bar(aes(fill = unit_region_pair), stat = "identity", alpha = 0.3) +
  geom_bar(aes(color = unit_region_pair), stat = "identity", fill = NA, alpha = 1) +
  scale_fill_jco() +
  scale_color_jco() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ unit_pair, scales = "free_y", ncol = 10) +
  labs(title = 'Cross-correlograms for unit pairs with min p_value and spike_counts parameters', 
       x = "Lag (ms)", y = "Spikes (n)") +
  theme_minimal() +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
        strip.text = element_text(size = 6))

# # Adjusted geom_text layer to correctly position the p-value text
# p <- p + geom_text(data = unique(data_long[c("unit_pair", "perm_p_val")]),
#                    aes(label = paste("p =", round(perm_p_val, 2)), x = 5, y = 5),
#                    stat = "identity", vjust = -0.5, size = 5, hjust = 1)

ggsave(output_path, plot = p, width = 30, height = 40, units = "in", limitsize = FALSE)
