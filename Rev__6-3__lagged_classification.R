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

dependent_fct <- 'NREM_fct'
dependent_num <- 'NREM_num'
independent_prefix <- 'tphate_'

input_data <- c('Cache/S01_Feb02_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul11_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul12_model_stack_10s_tphate3.csv',
                'Cache/S05_Jul13_model_stack_10s_tphate3.csv')

output_plot_1 <- 'Results/classification_nrem_performance_plot_individual___nolag_svm.svg'
output_plot_2 <- 'Results/classification_nrem_performance_plot_grouped___nolag_svm.svg'

unit_regions <- c('cla', 'acc', 'amy')
unit_colors <- c('#E28DB8', '#A67A77', '#7BA387')

#### Munging ####
# Load and row-bind data
data <- bind_rows(lapply(input_data, read_csv))
data$subset <- paste0(data$subset, '_', data$recording)

# Recode NREM to factor
data <- data %>% 
  mutate(NREM_num = NREM,
         NREM_fct = factor(ifelse(NREM == 1, "NREM", "WREM"), 
                           levels = c("NREM", "WREM")))

# Remove subsets with unit_num < min_units
data <- data[data$unit_num >= min_units,]

# Select lagged vs. standard
if (!lagged) {
  # Pattern to match columns starting with independent_prefix and ending with a number not equal to 1
  pattern <- paste0("^", independent_prefix, ".*[^1]$")
  
  # Identify columns to remove
  cols_to_remove <- grep(pattern, names(data), value = TRUE)
  
  # Remove identified columns from the dataframe
  data <- data[ , !(names(data) %in% cols_to_remove)]
}

#### Functions ####
classification_fx <- function(training_data, testing_data, dependent_fct, dependent_num, indepedent_prefix) {
  
  # Identify independent variable columns based on the prefix
  independent <- grep(paste0("^", independent_prefix), names(training_data), value = TRUE)
  
  # Construct the formula string
  formula_str <- paste(dependent_fct, "~", paste(independent, collapse = " + "))
  print(formula_str)
  formula <- as.formula(formula_str)
  
  # Define the model specification
  # model_spec <-
  #   boost_tree() %>%
  #   set_engine("xgboost") %>%
  #   set_mode("classification")
  
  # model_spec <-
  #   logistic_reg() %>%
  #   set_engine("glm") %>%
  #   set_mode("classification")
  
  model_spec <-
    svm_rbf() %>%
    set_engine("kernlab") %>%
    set_mode("classification")  
  
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
    predict(final_fit, new_data = testing_data, type = 'prob') %>%
    bind_cols(testing_data) %>%
    rename(observed = dependent_fct, predicted = .pred_NREM)
  
  # Calculate AUC and AP
  roc <- predictions %>%
    roc_auc(truth = observed, predicted)
  
  ap <- predictions %>%
    pr_auc(truth = observed, predicted)
  
  testing_metrics <- tibble(roc_auc = roc$.estimate, ap_auc = ap$.estimate)
  
  # Return a list containing both the performance metrics and the selected columns from predictions
  results_list <- list(testing_metrics, select(predictions, observed, predicted))
  
  return(results_list)
}

#### Modeling ####

# Initialize an empty list to store performance metrics for each subset
performance_list <- list()

# Loop through each unique subset
for (subset_selection in unique(data$subset)) {
  
  # Filter for the current subset
  subset_data <- filter(data, subset == subset_selection)
  
  # Split the data into training and testing sets randomly (70/30 split)
  data_split <- initial_split(subset_data, prop = 0.7)
  training_data <- training(data_split)
  testing_data <- testing(data_split)
  
  # Call classification_fx and extract the results
  results <- classification_fx(training_data, testing_data, dependent_fct, dependent_num, independent_prefix)
  performance_metrics <- results[[1]]
  predictions <- results[[2]]
  
  # Store performance metrics with the subset label
  performance_list[[subset_selection]] <- 
    performance_metrics %>% 
    mutate(subset = subset_selection)
}

# Combine performance metrics from all subsets into a single dataframe
performance <- bind_rows(performance_list, .id = "subset")
performance

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

# Apply the function to 'performance'
performance <- create_new_columns(performance)

# Create summary of model performance
performance_sum <- performance %>%
  group_by(unit_region, response) %>%
  summarise(
    roc_auc = mean(roc_auc, na.rm = TRUE),
    ap_auc = mean(ap_auc, na.rm = TRUE)
  ) %>%
  ungroup()  

# Model Performance Plot
p1 <- 
  ggplot(performance, aes(x = roc_auc, y = ap_auc, shape = response, color = unit_region)) +
  geom_point(size = 4) +
  scale_color_manual(values = setNames(unit_colors, unit_regions)) + 
  labs(title = 'Model performance by region, response, and recording',
       x = "AUC of ROC Curve",
       y = "AUC of PR Curve") + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
p1
ggsave(output_plot_1, plot = p1, width = 6, height = 4)

# Model Performance Plot - Mean
p2 <- 
  ggplot(performance_sum, aes(x = roc_auc, y = ap_auc, shape = response, color = unit_region)) +
  geom_point(aes(x = roc_auc, y = ap_auc), size = 6) +
  scale_color_manual(values = setNames(unit_colors, unit_regions)) + 
  labs(title = 'Mean model performance by region and response',
       x = "AUC of ROC Curve",
       y = "AUC of PR Curve") + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
p2
ggsave(output_plot_2, plot = p2, width = 6, height = 4)