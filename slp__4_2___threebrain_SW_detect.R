library(tidyverse)
library(threeBrain)
library(raveio)

#### Configuration ----
setwd('Z:/Layton/Sleep_083023')

# List of input files and unique prefix to make electrode numbers unique
# input <- list(
#   list(number = 1, 
#        elec_path = 'Cache/Subject01/Feb02/S01_electrodes_tkrRAS.csv', 
#        selected_path = 'Cache/Subject01/Feb02/S01_event_rates.csv'))

input <- list(
  list(number = 1, 
       elec_path = 'Cache/Subject05/S05_electrodes_tkrRAS.csv', 
       selected_path = 'Cache/Subject05/Jul13/S05_event_rates.csv'))

#### Load and merge elec info and coordinates ----

# Create blank tables to append data
elec_info_columns <- c('Electrode', 
                       'Coord_x', 
                       'Coord_y', 
                       'Coord_z', 
                       'Label', 
                       'MNI305_x', 
                       'MNI305_y', 
                       'MNI305_z')
elec_info <- data.frame(matrix(vector(), nrow = 0, ncol = length(elec_info_columns)))
colnames(elec_info) <- elec_info_columns

elec_values_columns <- c('Electrode', 'region')
elec_values <- data.frame(matrix(vector(), nrow = 0, ncol = length(elec_values_columns)))
colnames(elec_values) <- elec_values_columns

# Load each file of electrode data
for (subject in input) {
  elec_info_temp <- read.csv(subject$elec_path)
  
  # Increase prefix to electrode numbers
  elec_info_temp$Electrode <- 1000 * subject$number + elec_info_temp$Electrode
  
  #### Format electrode value object ----
  elec_values_temp <-
    read.csv(subject$selected_path) %>%
    filter(sw >= 3) %>%
    select(Channel, sw)
  
  elec_values_temp <- merge(elec_values_temp, 
                            elec_info_temp[,c('Electrode', 'Label')],
                            by.x = 'Channel', 
                            by.y = 'Label')
  
  elec_values_temp <- elec_values_temp[, c('Electrode', 'sw')]
  colnames(elec_values_temp) <- elec_values_columns
  
  # Append processed data to blank table
  elec_info <- rbind(elec_info, elec_info_temp)
  elec_values <- rbind(elec_values, elec_values_temp)
}

#### Import brain template ----

# download_template_subject('fsaverage')
# download_template_subject('cvs_avg35_inMNI152')
# set_default_template('N27')
template_path <- file.path(default_template_directory(), 'N27')

# freesurfer_brain2() will eventually be deprecated by threeBrain()
brain <- freesurfer_brain2(
  fs_subject_folder = template_path,
  subject_name = 'N27', 
  surface_types = c('pial',
                    'white',
                    'smoothwm'),
  atlas_types = c('aparc.DKTatlas+aseg')
)

brain$set_electrodes(electrodes = elec_info)
brain$set_electrode_values(elec_values)

#### Plot ----
plot(brain)

