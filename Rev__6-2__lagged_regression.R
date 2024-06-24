#### Configuration ####
library(tidyverse)
library(tidymodels)
library(easystats)
library(xgboost)

setwd('G:/My Drive/Residency/Research/Lab - Damisah/Project - Sleep/Revisions')

set.seed(42)

#### Parameters ####
lagged <- FALSE
min_units <- 4
nrem_only <- FALSE

dependent <- 'power'
independent_prefix <- 'tphate_'

input_data <- c('Cache/S01_Feb02_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul11_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul12_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul13_model_stack_10s_tphate3.csv')

output_plot_1 <- 'Results/regression_swa_performance_plot_individual___nolag_svm.svg'
output_plot_2 <- 'Results/regression_swa_performance_plot_grouped___nolag_svm.svg'
output_plot_3 <- 'Results/regression_swa_calibration_plot____nolag_svm.svg'

unit_regions <- c('cla', 'acc', 'amy')
unit_colors <- c('#E28DB8', '#A67A77', '#7BA387')

#### Munging ####
# Load and row-bind data
data <- bind_rows(lapply(input_data, read_csv))
data$subset <- paste0(data$subset, '_', data$recording)

# Remove subsets with unit_num < min_units
data <- data[data$unit_num >= min_units,]

# Select lagged vs. standard regression
if (!lagged) {
  # Pattern to match columns starting with independent_prefix and ending with a number not equal to 1
  pattern <- paste0("^", independent_prefix, ".*[^1]$")
  
  # Identify columns to remove
  cols_to_remove <- grep(pattern, names(data), value = TRUE)
  
  # Remove identified columns from the dataframe
  data <- data[ , !(names(data) %in% cols_to_remove)]
}

# Select NREM epochs only vs. all epochs
if (nrem_only) {
  data <- data[data$NREM == 1,]
}

#### Functions ####
regression_fx <- function(training_data, testing_data, dependent, indepedent_prefix) {
  
  # Identify independent variable columns based on the prefix
  independent <- grep(paste0("^", independent_prefix), names(training_data), value = TRUE)
  
  # Construct the formula string
  formula_str <- paste(dependent, "~", paste(independent, collapse = " + "))
  print(formula_str)
  formula <- as.formula(formula_str)
  
  # Define the model specification
  
  # model_spec <-
  #   boost_tree() %>%
  #   set_engine("xgboost") %>%
  #   set_mode("regression")
  
  # model_spec <-
  #   linear_reg() %>%
  #   set_engine("glm") %>%
  #   set_mode("regression")
  
  # model_spec <- gen_additive_mod() %>%
  #   set_engine("mgcv") %>%
  #   set_mode("regression")

  model_spec <-
    svm_rbf() %>%
    set_engine("kernlab") %>%
    set_mode("regression")
  
  # Define the recipe with the dynamically constructed formula
  recipe <- recipe(formula, data = training_data)
  
  # Create a 10-fold cross-validation resampling plan for the training data
  cv_folds <- vfold_cv(training_data, v = 10, repeats = 1)
  
  # Create the workflow
  workflow <- workflow() %>%
    add_model(model_spec) %>%
    add_recipe(recipe)
  
  # Train the model with 10-fold cross-validation on the training data
  fit <- workflow %>%
    fit_resamples(resamples = cv_folds)
  
  # Optionally, re-train the model on the entire training dataset and evaluate on the testing dataset
  final_fit <- workflow %>%
    fit(data = training_data)
  
  # Use the final model to make predictions on the testing set
  predictions <- 
    predict(final_fit, new_data = testing_data) %>%
    bind_cols(testing_data) %>% 
    rename(observed = !!sym(dependent), predicted = .pred)
  
  # Calculate performance metrics using the metrics() function
  testing_metrics <- 
    predictions %>%
    metrics(truth = observed, estimate = predicted)
  
  # Return a list containing both the performance metrics and the selected columns from predictions
  results_list <- list(testing_metrics, select(predictions, observed, predicted))
  
  return(results_list)
}

#### Modeling ####

# Initialize an empty list to store performance metrics for each subset
performance_list <- list()

# Initialize an empty DataFrame to store summarized data for all subsets
calibration <- tibble()

# Loop through each unique subset
for (subset_selection in unique(data$subset)) {
  
  # Filter for the current subset
  subset_data <- filter(data, subset == subset_selection)
  
  # Split the data into training and testing sets randomly (70/30 split)
  data_split <- initial_split(subset_data, prop = 0.7)
  training_data <- training(data_split)
  testing_data <- testing(data_split)
  
  # Call regression_fx and extract the results
  results <- regression_fx(training_data, testing_data, dependent, independent_prefix)
  performance_metrics <- results[[1]]
  predictions <- results[[2]]
  
  # Store performance metrics with the subset label
  performance_list[[subset_selection]] <- 
    performance_metrics %>% 
    mutate(subset = subset_selection)
  
  # Create the 'decile' column based on 'observed'
  predictions <- predictions %>%
    mutate(decile = ntile(observed, 10))
  
  # Summarize by 'decile' to get the average of 'observed' and 'predicted'
  calibration_data <- predictions %>%
    group_by(decile) %>%
    summarise(observed = mean(observed), predicted = mean(predicted), .groups = 'drop') %>%
    mutate(subset = subset_selection)
  
  # Append the summarized data to the final DataFrame
  calibration <- bind_rows(calibration, calibration_data)
}

# Combine all performance metrics into a single DataFrame
performance <- 
  bind_rows(performance_list, .id = "Subset") %>%
  select(subset, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

#### Plotting ####

# Function to create the new columns based on 'subset'
create_new_columns <- function(df) {
  # Split the 'subset' column by underscores
  split_subsets <- strsplit(as.character(df$subset), "_")
  
  # Extract 'unit_region' (the second element) and 'unit_subset' (the first two elements joined with an underscore)
  df$unit_region <- sapply(split_subsets, function(x) x[2])
  df$response <- sapply(split_subsets, function(x) x[1])
  df$unit_subset <- sapply(split_subsets, function(x) paste(x[1], x[2], sep = "_"))
  
  return(df)
}

# Apply the function to both 'performance' and 'calibration' dataframes
performance <- create_new_columns(performance)
calibration <- create_new_columns(calibration)

# Create summary of model performance
performance_sum <- performance %>%
  group_by(unit_region, response) %>%
  summarise(
    rmse = mean(rmse, na.rm = TRUE),
    mae = mean(mae, na.rm = TRUE)
  ) %>%
  ungroup()  

# Model Performance Plot
p1 <- 
  ggplot(performance, aes(x = rmse, y = mae, shape = response, color = unit_region)) +
  geom_point(size = 4) +
  scale_color_manual(values = setNames(unit_colors, unit_regions)) + 
  labs(title = 'Model performance by region, response, and recording',
       x = "Root Mean Square Error (RMSE)",
       y = "Mean Absolute Error (MAE)") + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
p1
ggsave(output_plot_1, plot = p1, width = 6, height = 4)

# Model Performance Plot - Mean
p2 <- 
  ggplot(performance_sum, aes(x = rmse, y = mae, shape = response, color = unit_region)) +
  geom_point(aes(x = rmse, y = mae, color = unit_region), 
             size = 6) +
  scale_color_manual(values = setNames(unit_colors, unit_regions)) + 
  labs(title = 'Mean model performance by region and response',
       x = "Root Mean Square Error (RMSE)",
       y = "Mean Absolute Error (MAE)") + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
p2
ggsave(output_plot_2, plot = p2, width = 6, height = 4)

# Model Calibration Plot
p3 <- ggplot(calibration, aes(x = predicted, y = observed, color = unit_region)) +
  geom_point(alpha = 0.3) +  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Unity line
  coord_fixed() +  # Ensure that the scales of the x and y axes are identical
  scale_color_manual(values = setNames(unit_colors, unit_regions)) + 
  facet_wrap(~ unit_subset) + 
  labs(title = 'Calibration stratified by decile of observed values within recordings', 
       x = "Mean predicted value in decile", y = "Mean observed value in decile") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, size = 1))
p3
ggsave(output_plot_3, plot = p3, width = 10, height = 6)
