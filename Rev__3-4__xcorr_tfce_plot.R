#### Configuration ####
library(tidyverse)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
window_min <- -50
bin_width <- 5

input_path <- 'Cache/All_xcorr_tfce_b5_w10_n1000_nrem_allpairs.csv'
meta_path <- 'Data/cell_types.csv'
output_path <- 'Results/xcorr_scatter_allnights_allpairs.svg'

alpha = 0.05

#### Munging ####

# Load data
data <- read.csv(input_path)

meta <- read.csv(meta_path)
meta <- meta[c("unit_id", "unit_region")]

# Merge for unit_1
data <- merge(data, meta, by.x = "unit_1", by.y = "unit_id", all.x = TRUE)
names(data)[names(data) == "unit_region"] <- "unit_region_1"

# Merge for unit_2
data <- merge(data, meta, by.x = "unit_2", by.y = "unit_id", all.x = TRUE)
names(data)[names(data) == "unit_region"] <- "unit_region_2"

# Convert bin lag to lag in milliseconds
data$max_tfce_lag_ms <- window_min + (data$max_tfce_lag - 1) * bin_width

# Find all unique region pairs
data$region_pair <- paste(data$unit_region_1, data$unit_region_2, sep = " - ")

# Log transform 
data$log10_max_tfce <- log10(data$max_tfce)

# Set significance level
data$signif <- ifelse(data$perm_p_val_fdr < alpha, 1, 0)
data$color_factor <- ifelse(data$signif == 1, "sig", "no_sig")

# Combine equivalent region pairs
data$region_pair <- sapply(data$region_pair, function(x) {
  pair <- strsplit(x, " - ")[[1]]
  sorted_pair <- sort(pair)
  paste(sorted_pair, collapse = " - ")
})

data$region_pair <- factor(data$region_pair, levels = c(
  "AMY - CLA", "ACC - CLA", "ACC - AMY", 
  "CLA - CLA", "AMY - AMY", "ACC - ACC"
))

# # Remove unwanted pair combinations
# unwanted_pairs <- c('ACC - AMY')
# 
# data <- data %>%
#   filter(!region_pair %in% unwanted_pairs)

#### Plotting ####

plot <- ggplot(data, aes(x = max_tfce_lag_ms, y = log10_max_tfce)) +
  geom_point(aes(color = color_factor, alpha = color_factor), size = 2.5) +
  scale_color_manual(values = c("sig" = "red", "no_sig" = "black")) +
  scale_alpha_manual(values = c("sig" = 1, "no_sig" = 0.15)) +
  scale_x_continuous(limits = c(-50, 50)) +
  #scale_y_continuous(limits = c(0, 7)) +
  facet_wrap(~region_pair) + # , scales = "free_x"
  labs(title = "Max TFCE for unit pairs (stratified by region pairs)", 
       x = "Lag (ms)", y = "Log10 Value (au)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75))

plot

ggsave(output_path, plot = plot, width = 8, height = 6, units = "in")
