#============================================
# Tidy Climate Data And Assign To Grid Cells
#============================================

library(ncdf4)
library(R.utils)
library(dplyr)
library(lubridate)
library(broom)

#-------------------------------------
# 1. Load in climate data
#-------------------------------------

# Set the path to your downloaded file
file_path <- "HadISST_sst.nc.gz"

# Decompress the file
decompressed_file_path <- gunzip(file_path, remove = FALSE)

# Open the NetCDF file
nc_data <- nc_open('HadISST_sst.nc')

#-------------------------------------
# 2. Tidy data
#-------------------------------------

# Extract variables
lat <- ncvar_get(nc_data, "latitude")
lon <- ncvar_get(nc_data, "longitude")
time <- ncvar_get(nc_data, "time")
temperature <- ncvar_get(nc_data, "sst")

# Convert time to a date format
start_date <- as.Date("1870-01-01")
time_dates <- start_date + time

#-------------------------------------
# 3. Split into grid cells
#-------------------------------------

# Create a data frame with latitude, longitude, and time combinations
grid <- expand.grid(lon = lon, lat = lat, time = time_dates)

# Flatten the temperature array into a vector
temperature_vector <- as.vector(temperature)

# Combine everything into a data frame
data_frame <- data.frame(
  lon = grid$lon,
  lat = grid$lat,
  time = grid$time,
  temperature = temperature_vector) %>%
  mutate(year = year(time)) %>%
  filter(year >= 1963, temperature >= -1.8) # values below -1.8 degrees appear to be ice measurements

# Remove rows with NA values (if any)
data_frame <- na.omit(data_frame)

# Calculate the average yearly temperature for each 1x1 grid cell
yearly_temperature <- data_frame %>%
  group_by(lon, lat, year) %>%
  summarise(max_temperature = max(temperature), mean_temperature = mean(temperature)) %>%
  select(lat, lon, year, max_temperature, mean_temperature)

# Function to create 4 new 0.5x0.5 grid cells from a 1x1 grid cell
split_grid <- function(lon, lat, max_temperature, mean_temperature, year) {
  # Create four new grid cells based on the original cell's boundaries:
  data.frame(
    lon = c(lon - 0.5, lon - 0.5, lon, lon),
    lat = c(lat - 0.5, lat, lat - 0.5, lat),
    max_temperature = rep(max_temperature, 4),
    mean_temperature = rep(mean_temperature, 4),
    year = rep(year, 4)
  )
}

# Apply the splitting function to each row of the data frame
expanded_data_frame <- yearly_temperature %>%
  rowwise() %>%
  do(split_grid(.$lon, .$lat, .$max_temperature, .$mean_temperature, .$year)) %>%
  ungroup()

#-------------------------------------
# 4. Write to file
#-------------------------------------

# Write data to .csv file
write.csv(expanded_data_frame, 'HADISST.csv', row.names = F)
