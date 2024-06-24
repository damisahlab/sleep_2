#### Configuration ####
library(tidyverse)
#library(statsExpressions)
library(tune)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

#### Parameters ####
long_path <- c('Data/S01_Feb02_long_scatter.csv',
               'Data/S05_Jul11_long_scatter.csv',
               'Data/S05_Jul12_long_scatter.csv',
               'Data/S05_Jul13_long_scatter.csv')

short_path <- c('Data/S01_Feb02_short_scatter.csv',
                'Data/S05_Jul11_short_scatter.csv',
                'Data/S05_Jul12_short_scatter.csv',
                'Data/S05_Jul13_short_scatter.csv')

stage_out_path_1 <- 'Results/mini_barplot_stage.svg'
swa_out_path_1 <- 'Results/mini_barplot_swa.svg'
sws_out_path_1 <- 'Results/mini_barplot_sws.svg'

stage_out_path_2 <- 'Results/triangle_scatterplot_stage.svg'
swa_out_path_2 <- 'Results/triangle_scatterplot_swa.svg'
sws_out_path_2 <- 'Results/triangle_scatterplot_sws.svg'

micro_regions <- c('CLA', 'AMY', 'ACC')
micro_colors <- c('#E28DB8', '#7BA387', '#A67A77')

#### Munging ####
# Long scatterplot data (sleep stage)
long_data <- data.frame()

for(path in long_path){
  
  # Extract date from file path
  date <- gsub(".*/(.*)/.*\\.csv$", "\\1", path)
  
  temp_data <- 
    read.csv(path) %>%
    select(unit_id, unit_region, unit_laterality, fr, stage) %>%
    mutate(unit_id = paste(unit_id, date, sep="_"))
  
  long_data <- rbind(long_data, temp_data)
  
}

stage_data <-
  long_data %>%
  mutate(unit_region = as.factor(unit_region)) %>%
  mutate(unit_region = fct_relevel(unit_region, micro_regions)) %>%
  mutate(stage = ifelse(stage == 1, 'n23', 'not_n23')) %>%
  group_by(unit_id, unit_region, stage) %>%
  summarize(FR = mean(fr),
            logFR = log2(mean(fr))) %>%
  pivot_wider(id_cols = c(unit_id, unit_region),
              names_from = c(stage),
              values_from = c(logFR))

# Short scatterplot data (SWA, SWS)
short_data <- data.frame()

for(path in short_path){
  
  # Extract date from file path
  date <- gsub(".*/(.*)/.*\\.csv$", "\\1", path)
  
  temp_data <- 
    read.csv(path) %>%
    select(unit_id, unit_region, unit_laterality, fr, swa_zscore, sws) %>%
    mutate(unit_id = paste(unit_id, date, sep="_"))
  
  short_data <- rbind(short_data, temp_data)
  
}

sws_data <-
  short_data %>%
  mutate(unit_region = as.factor(unit_region)) %>%
  mutate(unit_region = fct_relevel(unit_region, micro_regions)) %>%
  mutate(SWS = case_when(sws >= mean(sws) ~ 'SWS',
                         sws == 0 ~ 'noSWS',
                         TRUE ~ 'midSWS')) %>% 
  group_by(unit_id, unit_region, SWS) %>%
  summarize(FR = mean(fr),
            logFR = log2(mean(fr))) %>%
  pivot_wider(id_cols = c(unit_id, unit_region),
              names_from = c(SWS),
              values_from = c(logFR))

swa_data <- 
  short_data %>%
  mutate(unit_region = as.factor(unit_region)) %>%
  mutate(unit_region = fct_relevel(unit_region, micro_regions)) %>%
  mutate(SWA = case_when(
    swa_zscore >= quantile(swa_zscore)[4] ~ 'SWA',
    swa_zscore <= quantile(swa_zscore)[2] ~ 'noSWA',
    TRUE ~ 'midSWA')) %>% 
  group_by(unit_id, unit_region, SWA) %>%
  summarize(FR = mean(fr),
            logFR = log2(mean(fr))) %>%
  pivot_wider(id_cols = c(unit_id, unit_region),
              names_from = c(SWA),
              values_from = c(logFR))

stage_data$unit_region <- factor(stage_data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
swa_data$unit_region <- factor(swa_data$unit_region, levels = c('CLA', 'AMY', 'ACC'))
sws_data$unit_region <- factor(sws_data$unit_region, levels = c('CLA', 'AMY', 'ACC'))

#### Statistics ####
# Sleep stage
chisq_test_result_stage <- stage_data %>%
  mutate(stats_group = ifelse(unit_region == "CLA", "CLA", "CTRL"),
         triangle = ifelse(n23 > not_n23, "lower", "upper")) %>%
  with(chisq.test(stats_group, triangle))

# SWA
chisq_test_result_swa <- swa_data %>%
  mutate(
    stats_group = ifelse(unit_region == "CLA", "CLA", "CTRL"),
    triangle = ifelse(SWA > noSWA, "lower", "upper")
  ) %>%
  with(chisq.test(stats_group, triangle))

# SWS
chisq_test_result_sws <- sws_data %>%
  mutate(
    stats_group = ifelse(unit_region == "CLA", "CLA", "CTRL"),
    triangle = ifelse(SWS > noSWS, "lower", "upper")
  ) %>%
  with(chisq.test(stats_group, triangle))

# Create a dataframe with raw p-values and labels
p_values_df <- data.frame(
  label = c("stage", "swa", "sws"),
  p_value = c(
    chisq_test_result_stage$p.value,
    chisq_test_result_swa$p.value,
    chisq_test_result_sws$p.value
  )
)

# Calculate FDR-corrected p-values
p_values_df$fdr_p_value <- p.adjust(p_values_df$p_value, method = "BH")

# Print Chi-Square results
print('Chi-Squared test | CLA vs. CTRL | Upper vs. Lower Triangle')
print(p_values_df)

#### Functions ####
# Function to compute the mean and 95% CI
compute_summary <- function(x) {
  n <- length(x)
  mu <- mean(x)
  stderr <- sd(x)/sqrt(n)
  ci_halfwidth <- stderr * qt(0.975, df = n-1)
  return(data.frame(mean = mu, lower = mu - ci_halfwidth, upper = mu + ci_halfwidth))
}

#### Mini Bar Plots ####
# Create a dataset that counts conditions for each unit_region
counts_data <- stage_data %>%
  group_by(unit_region) %>%
  summarise(
    n23_greater = sum(n23 > not_n23),
    n23_not_greater = sum(n23 <= not_n23),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(n23_greater, n23_not_greater),
               names_to = "condition",
               values_to = "count") %>%
  mutate(condition = ifelse(condition == "n23_greater", "n23 > not_n23", "n23 <= not_n23"))

# Create the mini bar plot with facets for each unit_region
mini_bar_stage <- ggplot(counts_data, aes(x = condition, y = count, fill = unit_region)) +
  geom_bar(alpha = 0.6, stat = "identity") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_text(aes(label = count, y = count/2), vjust = 0.5, color = "black") +
  scale_fill_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(~ unit_region) +  # Facet by unit_region
  theme_void() +
  theme(legend.position = "none")

mini_bar_stage

ggsave(file = stage_out_path_1, 
       width = 3, height = 1.5, 
       dpi = 300, limitsize = FALSE)

# Create a dataset that counts conditions for each unit_region
swa_counts <- swa_data %>%
  group_by(unit_region) %>%
  summarise(
    SWA_greater = sum(SWA > noSWA),
    SWA_not_greater = sum(SWA <= noSWA),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(SWA_greater, SWA_not_greater),
               names_to = "condition",
               values_to = "count") %>%
  mutate(condition = ifelse(condition == "SWA_greater", "SWA > noSWA", "SWA <= noSWA"))

# Create the mini bar plot with facets for each unit_region
mini_bar_swa <- ggplot(swa_counts, aes(x = condition, y = count, fill = unit_region)) +
  geom_bar(alpha = 0.6, stat = "identity") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_text(aes(label = count, y = count/2), vjust = 0.5, color = "black") +
  scale_fill_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(~ unit_region) +
  theme_void() +
  theme(legend.position = "none")

mini_bar_swa

ggsave(file = swa_out_path_1, 
       width = 3, height = 1.5, 
       dpi = 300, limitsize = FALSE)

# Create a dataset that counts conditions for each unit_region
sws_counts <- sws_data %>%
  group_by(unit_region) %>%
  summarise(
    SWS_greater = sum(SWS > noSWS),
    SWS_not_greater = sum(SWS <= noSWS),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(SWS_greater, SWS_not_greater),
               names_to = "condition",
               values_to = "count") %>%
  mutate(condition = ifelse(condition == "SWS_greater", "SWS > noSWS", "SWS <= noSWS"))

# Create the mini bar plot with facets for each unit_region
mini_bar_sws <- ggplot(sws_counts, aes(x = condition, y = count, fill = unit_region)) +
  geom_bar(alpha = 0.6, stat = "identity") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_text(aes(label = count, y = count/2), vjust = 0.5, color = "black") +
  scale_fill_manual(values = setNames(micro_colors, micro_regions)) +
  facet_wrap(~ unit_region) +
  theme_void() +
  theme(legend.position = "none")

mini_bar_sws

ggsave(file = sws_out_path_1, 
       width = 3, height = 1.5, 
       dpi = 300, limitsize = FALSE)

#### Triangle Scatterplots ####
# Summary
sum_data <- 
  stage_data %>%
  group_by(unit_region) %>%
  summarise(across(c(x = n23, y = not_n23), 
                   compute_summary, 
                   .names = "{.col}_")) %>%
  tidyr::unnest(cols = c(x_, y_), names_sep = "")

# Plot
stage_plot <- 
  ggplot() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed') + 
  geom_point(data = stage_data, 
             aes(x = n23, y = not_n23, 
                 group = unit_region,
                 color = unit_region),
             size = 3, alpha = 0.5) +
  geom_linerange(data = sum_data,
                 aes(y = y_mean,
                     xmin = x_lower,
                     xmax = x_upper),
                 color = 'black',
                 linewidth = 0.75) +
  geom_linerange(data = sum_data,
                 aes(x = x_mean,
                     ymin = y_lower,
                     ymax = y_upper),
                 color = 'black',
                 linewidth = 0.75) +
  scale_discrete_manual(values = micro_colors,
                        aesthetics = c('color')) + 
  facet_wrap(unit_region ~ .) + 
  coord_obs_pred() + 
  labs(x = 'log2 FR in NREM', 
       y = 'log2 FR not in NREM') +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.4),
        strip.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', 
                                    fill = NA, 
                                    size = 1))

stage_plot

ggsave(file = stage_out_path_2, 
       width = 9, height = 3, 
       dpi = 300, limitsize = FALSE)

# Summary
sum_data <- 
  swa_data %>%
  group_by(unit_region) %>%
  summarise(across(c(x = SWA, y = noSWA), 
                   compute_summary, 
                   .names = "{.col}_")) %>%
  tidyr::unnest(cols = c(x_, y_), names_sep = "")

# Plot
swa_plot <- 
  ggplot() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed') + 
  geom_point(data = swa_data, 
             aes(x = SWA, y = noSWA, 
                 group = unit_region,
                 color = unit_region),
             size = 3, alpha = 0.5) +
  geom_linerange(data = sum_data,
                 aes(y = y_mean,
                     xmin = x_lower,
                     xmax = x_upper),
                 color = 'black',
                 linewidth = 0.75) +
  geom_linerange(data = sum_data,
                 aes(x = x_mean,
                     ymin = y_lower,
                     ymax = y_upper),
                 color = 'black',
                 linewidth = 0.75) +
  scale_discrete_manual(values = micro_colors,
                        aesthetics = c('color')) + 
  facet_wrap(unit_region ~ .) + 
  coord_obs_pred() + 
  labs(x = 'log2 FR high SWA', 
       y = 'log2 FR low SWA') +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.4),
        strip.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', 
                                    fill = NA, 
                                    size = 1))

swa_plot

ggsave(file = swa_out_path_2, 
       width = 9, height = 3, 
       dpi = 300, limitsize = FALSE)

# Summary
sum_data <- 
  sws_data %>%
  group_by(unit_region) %>%
  summarise(across(c(x = SWS, y = noSWS), 
                   compute_summary, 
                   .names = "{.col}_")) %>%
  tidyr::unnest(cols = c(x_, y_), names_sep = "")

# Plot
sws_plot <- 
  ggplot() +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 'dashed') + 
  geom_point(data = sws_data, 
             aes(x = SWS, y = noSWS, 
                 group = unit_region,
                 color = unit_region),
             size = 3, alpha = 0.5) +
  geom_linerange(data = sum_data,
                 aes(y = y_mean,
                     xmin = x_lower,
                     xmax = x_upper),
                 color = 'black',
                 linewidth = 0.75) +
  geom_linerange(data = sum_data,
                 aes(x = x_mean,
                     ymin = y_lower,
                     ymax = y_upper),
                 color = 'black',
                 linewidth = 0.75) +
  scale_discrete_manual(values = micro_colors,
                        aesthetics = c('color')) + 
  facet_wrap(unit_region ~ .) + 
  coord_obs_pred() + 
  labs(x = 'log2 FR high SW presence', 
       y = 'log2 FR no SW presence') +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.4),
        strip.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black', 
                                    fill = NA, 
                                    size = 1))

sws_plot

ggsave(file = sws_out_path_2, 
       width = 9, height = 3, 
       dpi = 300, limitsize = FALSE)