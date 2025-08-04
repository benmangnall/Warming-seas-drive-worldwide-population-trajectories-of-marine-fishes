#================================
# Calculate Grid Cell Abundances
#================================

# Load required packages
library(dplyr)
library(data.table)
library(tidyr)
library(DMwR2)
library(rfishbase)
library(broom)
library(rfishbase)

#-------------------------------------
# 1. Read in data and extract surveys
#-------------------------------------

# Load the FishGlob database
fish_data <- fread('Data.csv', sep = ',') 

# Extract unique surveys
unique_surveys <- fish_data %>%
  select(survey) %>%
  distinct()

#-------------------------------------
# 2. Calculate mean abundances for each species in each grid cell in each year 
#-------------------------------------

# Outer loop through each survey
for (survey_name in unique_surveys$survey) {
  
  # Filter data for the current survey
  survey_data <- fish_data %>%
    filter(survey == survey_name)
  
  #-------------------------------------
  # 2.1 Create grid cells for each survey
  #-------------------------------------
  
  # Extract unique haul coordinates for grid cell creation
  haul_coordinates <- survey_data %>%
    select(region, survey, haul_id, latitude, longitude) %>%
    distinct() %>%
    mutate(latitude_gridcell = floor(latitude * 2) / 2,
           longitude_gridcell = floor(longitude * 2) / 2) 
  
  # Extract unique grid cells
  unique_gridcells <- haul_coordinates %>%
    select(latitude_gridcell, longitude_gridcell) %>%
    distinct()
  
  #-------------------------------------
  # 2.2 Loop for finding mean abundances
  #-------------------------------------
  
  # Loop through each grid cell
  for (i in seq_len(nrow(unique_gridcells))) {
    
    # Extract current grid cell's latitude and longitude
    current_gridcell <- unique_gridcells[i, ]
    lat_grid <- current_gridcell$latitude_gridcell
    lon_grid <- current_gridcell$longitude_gridcell
    
    # Filter haul coordinates for the current grid cell
    cell_haul_ids <- haul_coordinates %>%
      filter(latitude_gridcell == lat_grid, longitude_gridcell == lon_grid) %>%
      pull(haul_id)
    
    # Filter fish_data for the hauls in the current grid cell
    gridcell_data <- survey_data %>%
      filter(haul_id %in% cell_haul_ids) %>%
      select(haul_id, year, accepted_name, abundance)
    
    # Identify all species ever found in the grid cell
    all_species_grid <- survey_data %>%
      filter(haul_id %in% cell_haul_ids) %>%
      select(accepted_name) %>%
      distinct()
    
    # Create all possible combinations of (year, species, haul_id) within the grid cell
    complete_grid <- survey_data %>%
      left_join(haul_coordinates, by = c("region", "survey", "haul_id", "latitude", "longitude")) %>%
      filter(latitude_gridcell == lat_grid & longitude_gridcell == lon_grid) %>%
      group_by(latitude_gridcell, longitude_gridcell, year) %>%
      summarise(haul_ids = list(unique(haul_id)), 
                species = list(unique(all_species_grid$accepted_name)), 
                .groups = "drop") %>%
      rowwise() %>%
      mutate(grid = list(expand.grid(haul_id = haul_ids, accepted_name = species, stringsAsFactors = FALSE))) %>%
      select(-haul_ids, -species) %>%
      unnest(cols = c(grid))
    
    # Join the expanded grid with the gridcell data to ensure all combinations
    abundances_with_zero <- complete_grid %>%
      full_join(gridcell_data, by = c('haul_id', 'year', 'accepted_name')) %>%
      mutate(abundance = ifelse(is.na(abundance), 0, abundance))
    
    # Calculate mean abundance within each grid cell, year, and species
    mean_abundance <- abundances_with_zero %>%
      select(-latitude_gridcell, -longitude_gridcell) %>%
      left_join(haul_coordinates, by = "haul_id") %>%
      group_by(region, survey, accepted_name, latitude_gridcell, longitude_gridcell, year) %>%
      summarise(mean_num = sum(abundance) / n_distinct(haul_id), .groups = "drop")
    
    # Write output to CSV, naming the file by survey, year, and grid cell
    output_filename <- paste0(survey_name, "_gridcell_", lat_grid, "_", lon_grid, ".csv")
    fwrite(mean_abundance, output_filename)
  }
}

#-------------------------------------
# 3. Join loop outputs together 
#-------------------------------------

# Get a list of all the CSV files (adjust the path if needed)
csv_files <- list.files(pattern = "_gridcell_.*\\.csv$")

# Read each CSV file and combine them into one data frame
all_abundance_mean <- csv_files %>%
  lapply(read.csv) %>%   # Read each CSV file into a list of data frames
  bind_rows()            # Combine all data frames into one

# Save the combined data frame to a single CSV file
write.csv(all_abundance_mean, "abundance_means.csv", row.names = F)

# Remove all intermediate files
file.remove(csv_files)

rm(list = ls())
gc()

#-------------------------------------
# 4. Join climate to abundance data 
#-------------------------------------

# Load in master data frame
fish_data <- fread('Data.csv', sep = ',') 

# Load in abundance data
abundance <- read.csv('abundance_means.csv', head = T, sep = ',')

# Load in climate data
climate <- read.csv('HADISST.csv', head = T, sep = ',') %>%
  rename(latitude_gridcell = lat, longitude_gridcell = lon)

# Join climate data to abundance data
abundance_climate <- abundance %>%
  left_join(climate, by = c('latitude_gridcell', 'longitude_gridcell', 'year'))

# select the unique grid cells and the required columns for imputation
columns_for_imputation <- abundance_climate %>%
  select(longitude_gridcell, latitude_gridcell, year, max_temperature, mean_temperature) %>%
  distinct()

# Perform imputation to fill in missing temperature values
imputed_data <- knnImputation(columns_for_imputation, k = 1) 

# Replace the old temperature values with the imputed ones
abundance_with_temperature <- abundance_climate %>%
  select(-max_temperature, -mean_temperature) %>%
  left_join(imputed_data, by = c('latitude_gridcell', 'longitude_gridcell', 'year'))

# Add a temperature change per year column
abundance_temperature <- abundance_with_temperature %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  mutate(temp_slope = lm(mean_temperature ~ year) %>% tidy() %>% filter(term == "year") %>% pull(estimate)) %>%
  filter(!is.na(temp_slope))

#-------------------------------------
# 5. Assign position in species range 
#-------------------------------------

#-------------------------------------
# 5.1 Northern hemisphere grid cells
#-------------------------------------

# Load the filtered GBIF data
gbif_data <- read.csv('gbif_data_northern.csv', head = T, sep = ',') %>%
  select(species, decimalLatitude) %>%
  rename(latitude = decimalLatitude, accepted_name = species)

# Get the northern hemisphere species ranges
north <- fish_data[fish_data$latitude > 0, ] %>%
  select(accepted_name, latitude)

# Remove rows in 'north' that are already in gbif
filtered_north <- north %>%
  anti_join(gbif_data, by = c("accepted_name", "latitude"))

# Add the abundance data latitudes on to gbif data
gbif_data <- bind_rows(gbif_data, filtered_north)

gbif_data <- gbif_data %>%
  group_by(accepted_name) %>%
  filter(n_distinct(latitude) > 100)

# Filter 1% extreme latitude values per species
filtered_gbif_data <- gbif_data %>%
  group_by(accepted_name) %>%
  filter(latitude > quantile(latitude, 0.001, na.rm = TRUE),
         latitude < quantile(latitude, 0.999, na.rm = TRUE),
         latitude > 0) %>%
  ungroup()

# Calculate min and max latitudes for each species
species_ranges <- filtered_gbif_data %>%
  group_by(accepted_name) %>%
  summarise(
    equator_latitude = min(latitude, na.rm = TRUE),
    poleward_latitude = max(latitude, na.rm = TRUE))

# Northern data
abundance_temperature_north <- abundance_temperature[abundance_temperature$latitude_gridcell >= 0, ]

# Join with other data 
full_data_north <- abundance_temperature_north %>%
  left_join(species_ranges, by = 'accepted_name') %>%
  mutate(position_in_range = (latitude_gridcell - equator_latitude) / (poleward_latitude - equator_latitude)) %>%
  select(-equator_latitude, -poleward_latitude) %>%
  filter(!is.na(position_in_range), position_in_range != 'Inf', position_in_range != '-Inf')

#-------------------------------------
# 5.2 Southern hemisphere grid cells
#-------------------------------------

# Load the filtered GBIF data
gbif_data <- read.csv('gbif_data_southern.csv', head = T, sep = ',') %>%
  select(species, decimalLatitude) %>%
  rename(latitude = decimalLatitude, accepted_name = species)

# Get the southern hemisphere species ranges
south <- fish_data[fish_data$latitude < 0, ] %>%
  select(accepted_name, latitude)

# Remove rows in 'south' that are already in gbif
filtered_south <- south %>%
  anti_join(gbif_data, by = c("accepted_name", "latitude"))

# Add the abundance data latitudes on to gbif data
gbif_data <- bind_rows(gbif_data, filtered_south)

gbif_data <- gbif_data %>%
  group_by(accepted_name) %>%
  filter(n_distinct(latitude) > 100)

# Filter 1% extreme latitude values per species
filtered_gbif_data <- gbif_data %>%
  group_by(accepted_name) %>%
  filter(latitude > quantile(latitude, 0.001, na.rm = TRUE),
         latitude < quantile(latitude, 0.999, na.rm = TRUE),
         latitude < 0) %>%
  ungroup()

# Calculate min and max latitudes for each species
species_ranges <- filtered_gbif_data %>%
  group_by(accepted_name) %>%
  summarise(
    poleward_latitude = min(latitude, na.rm = TRUE),
    equator_latitude = max(latitude, na.rm = TRUE)) 

# Southern data
abundance_temperature_south <- abundance_temperature[abundance_temperature$latitude_gridcell < 0, ]

# Join with other data 
full_data_south <- abundance_temperature_south %>%
  left_join(species_ranges, by = 'accepted_name') %>%
  mutate(position_in_range = (latitude_gridcell - equator_latitude) / (poleward_latitude - equator_latitude)) %>%
  select(-poleward_latitude, -equator_latitude) %>%
  filter(!is.na(position_in_range), position_in_range != 'Inf', position_in_range != '-Inf')

#-------------------------------------
# 5.3 Join both together and write to file
#-------------------------------------

# Bind both hemispheres together
all <- bind_rows(full_data_north, full_data_south)

# Write to .csv
write.csv(all, 'abundances.csv', row.names = F)

#-------------------------------------
# 6. Obtain species trait information
#-------------------------------------

# Load in the data
all <- fread('abundances.csv', sep = ',')

# Get a species list
species_list <- unique(all$accepted_name)

# Validate the species names
validated_names <- unique(species(species_list))  %>%
  select(Species)

# Fetch size data for validated species
adult_lh <- species(validated_names$Species, fields = c("Species", "Length", "DemersPelag", "DepthRangeDeep")) %>%
  rename(Habitat = DemersPelag)

# Get larval phase data
larval_lh <- larvae(validated_names$Species, fields = c("Species", "PlaceofDevelopment"))

# Get order information
order_data <- load_taxa(server=getOption("fishbase")) # get taxonomic data from fishbase
order_data <- order_data[,c(2,6)]

# Get Fecundity information
fecundity_lh <- fecundity(validated_names$Species) %>%
  select(Species, FecundityMin, FecundityMax, FecundityMean) %>%
  group_by(Species) %>%
  summarise(
    FecundityMin  = mean(FecundityMin,  na.rm = TRUE),
    FecundityMax  = mean(FecundityMax,  na.rm = TRUE),
    FecundityMean = mean(FecundityMean, na.rm = TRUE)) %>%
  mutate(Fecundity = rowMeans(select(., FecundityMin, FecundityMax, FecundityMean),na.rm = TRUE)) %>%
  select(Species, Fecundity)

# Get trophic level information
trophic_lh <- estimate(validated_names$Species, fields = c("Species", "Troph"))

# Get growth rate information
growth_lf <- estimate(validated_names$Species, fields = c("Species", "K"))

# Join to full species list
life_history <- validated_names %>%
  left_join(adult_lh, by = "Species") %>%
  left_join(larval_lh, by = "Species") %>%
  left_join(order_data, by = "Species") %>%
  left_join(fecundity_lh, by = "Species") %>%
  left_join(trophic_lh, by = "Species") %>%
  left_join(growth_lf, by = "Species") %>%
  rename(accepted_name = Species)

# Correct naming formats
life_history <- life_history %>%
  mutate(Order = ifelse(Order == "Carangaria/misc", "Carangiformes", Order)) %>%
  mutate(Order = ifelse(Order == "Eupercaria/misc", "Acanthuriformes", Order)) %>%
  mutate(Order = sub("/.*", "", Order)) %>%
  mutate(Habitat = case_when(
    Habitat %in% ('bathydemersal') ~ 'demersal',
    Habitat %in% c('bathypelagic', 'pelagic-neritic', 'pelagic-oceanic') ~ 'pelagic',  TRUE ~ Habitat)) %>%
  mutate(PlaceofDevelopment = case_when(
    PlaceofDevelopment %in% c('attached to parental body', 'in brood pouch', 'in close association with substrate', 
                              'in closed nest (e.g. burrow or tunnel)', 'in female (livebearers)', 'in mouth (mouthbrooders)', 
                              'in open nest') ~ 'Demersal', TRUE ~ PlaceofDevelopment))

# Write to .csv file
write.csv(life_history, 'life_history.csv', row.names = F)

rm(list = ls())
gc()

#-------------------------------------
# 7. Model abundance temperature response
#-------------------------------------

# Load the data
abundance_mean <- fread('abundances.csv', sep = ',') 

# Only include grid cells with at least 5 years of sampling
filtered <- abundance_mean %>%
  group_by(survey, latitude_gridcell, longitude_gridcell) %>%
  filter(n_distinct(year) >= 5)

# Standardizing the data (Z transformation)
standardize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Apply Z transformation to mean abundance
abundance_standardised <- filtered %>%
  group_by(latitude_gridcell, longitude_gridcell, accepted_name, survey) %>%
  mutate(mean_abundance = standardize(mean_num)) %>%
  na.omit()

# Add a time series duration column 
abund_time <- abundance_standardised %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  mutate(duration = max(year) - min(year) + 1)

# Calculate slopes for abundance ~ year and abundance ~ temperature
slope <- abund_time %>%
  group_by(region, survey, latitude_gridcell, longitude_gridcell, accepted_name, position_in_range, duration, temp_slope) %>%
  summarise(abundance_temp_slope = lm(mean_num ~ mean_temperature) %>% tidy() %>% filter(term == 'mean_temperature') %>% pull(estimate)) %>% 
  mutate(abundance_temp_binary = as.integer(abundance_temp_slope > 0)) 

# Calculate coefficients for each species in each grid cell 
coefficients <- abund_time %>%
  group_by(region, survey, latitude_gridcell, longitude_gridcell, accepted_name, position_in_range, duration, temp_slope) %>%
  summarise(abundance_mean_temp_correlation = cor(mean_abundance, mean_temperature, use = 'complete.obs'))

# Join both methods of calculating abundance change into one table
full_change <- coefficients %>%
  left_join(slope, by = c('region', 'accepted_name', 'latitude_gridcell', 'longitude_gridcell', 'position_in_range', 'survey', 'duration', 'temp_slope')) %>%
  filter(position_in_range > 0, position_in_range < 1) # Only keep valid positions

#-------------------------------------
# 8. Add trait information to abundance responses
#-------------------------------------

# Load in life history data
life_history <- fread('life_history.csv', sep = ',')

# Add life history information to data set
full_change <- full_change %>%
  left_join(life_history, by = "accepted_name") 

# Convert accepted_name to numeric 
full_change$accepted_name_numeric <- as.numeric(as.factor(full_change$accepted_name))

# Create a log Length column
full_change$logLength <- scale(log10(full_change$Length))

# Create a log depth column
full_change$logDepth <- scale(log10(full_change$DepthRangeDeep))

# Create a log fecundity column
full_change$logFecundity <- scale(log10(full_change$Fecundity))

# Scale trophic level
full_change$Troph <- scale(full_change$Troph)

# Create a log growth rate column
full_change <- full_change %>% filter(K != 0)
full_change$logK <- scale(log10(full_change$K)) 

# Remove now redundant columns
full_change <- full_change %>%
  select(-Length, -DepthRangeDeep, -Fecundity, -K)

# Save file 
write.csv(full_change, 'final_abundance_changes.csv', row.names = F)

# end of script

############ Summary data ############

fish_data <- fread('data.csv', sep = ',')

num_hauls <- n_distinct(fish_data$haul_id)

num_surveys <- n_distinct(fish_data$survey)

num_abundance <- sum(fish_data$abundance)

num_species <- n_distinct(All_regions$accepted_name)


#number of hauls 
haul_coordinates <- fish_data %>%
  select(region, survey, haul_id, latitude, longitude, year) %>%
  distinct() %>%
  mutate(latitude_gridcell = floor(latitude * 2) / 2,
         longitude_gridcell = floor(longitude * 2) / 2) 

duration <- haul_coordinates %>%
  distinct() %>%
  group_by(latitude_gridcell, longitude_gridcell) %>%
  mutate(durations = max(year) - min(year) +1)

haul_num <- duration %>%
  filter(durations >= 5) 

haul_number <- n_distinct(haul_num$haul_id)
