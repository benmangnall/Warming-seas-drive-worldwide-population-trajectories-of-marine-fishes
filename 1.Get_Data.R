#-----------------------------------------------------
# Join All Data Sets Together To Prepare For Analysis
#-----------------------------------------------------

library(lubridate)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

#-------------------------------------
# 1. Collate IMOS and NIWA data sets
#-------------------------------------

# Load in IMOS survey
IMOS_reef <- read.csv('IMOS_Reef_Abundance.csv', header = T, sep = ',', skip = 71) %>%
  mutate(survey = 'IMOS-reef')

# Load in IMOS survey
IMOS_cryptobentic <- read.csv('IMOS_Cryptobenthic_Abundance.csv', header = T, sep = ',', skip = 71) %>%
  mutate(survey = 'IMOS-cryptobenthic')

# Bind IMOS rows together
IMOS <- bind_rows(IMOS_reef, IMOS_cryptobentic)

# Treat survey_id as a string
IMOS$survey_id <- as.character(IMOS$survey_id) 

# Correct naming format
IMOS <- IMOS %>%
  mutate(year = year(survey_date), species_name = tolower(species_name)) %>% # same grid cell gets listed as different regions 
  rename(haul_id = survey_id, accepted_name = species_name, abundance = total, region = country) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  mutate(region = if_else(region == "Morocco", "Spain", region))

# contains the multiple abundance records for the same species in the same dive, therefore needed to be summed together
IMOS_summed <- IMOS %>%
  group_by(region, survey, haul_id, latitude, longitude, year, accepted_name, depth) %>%
  summarise(abundance = sum(abundance))

# Load in New Zealand Fish and squid data
NIWA_NZ <- read.csv('IMOS_NZ_Fish.csv', header = T, sep = ',')

# Extract longitude and latitude from 'geom'
NIWA <- NIWA_NZ %>%
  mutate(occurrenceremarks = as.numeric(gsub("[^0-9.]", "", occurrenceremarks))) %>%
  mutate(
    geom = gsub("POINT \\(", "", geom),  # Remove 'POINT ('
    geom = gsub("\\)", "", geom)        # Remove ')'
  ) %>%
  separate(geom, into = c("longitude", "latitude"), sep = " ", convert = TRUE)

# Select only the year 
NIWA$time <- format(as.Date(NIWA$time), "%Y")

# Standardize naming formats
NIWA <- NIWA %>%
  select(-year) %>%
  rename(region = waterbody, survey = institutioncode, haul_id = fieldnumber, year = time, accepted_name = scientificname, abundance = occurrenceremarks, depth = maximumdepthinmeters) %>%
  mutate(accepted_name = tolower(accepted_name), region = case_when(region == 'Pacific Ocean' ~ 'New Zealand')) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth)

# Ensure 'year' in NIWA is numeric
NIWA$year <- as.numeric(NIWA$year)

# Contains patchy data pre 1979, best to remove 
NIWA <- NIWA %>%
  filter(abundance != '', year > 1978)

# Bind all Reef life survey data together
RLS <- bind_rows(IMOS_summed, NIWA)

#-------------------------------------
# 2. Process FishGLOB data sets
#-------------------------------------

# Load the FishGlob database
fishglob <- fread('FishGlob_.csv', sep = ',') 

# Standardize accepted name format, and use the most frequent abundance metric for each survey
fishglob <- fishglob %>%
  mutate(
    accepted_name = tolower(accepted_name),
    num_cpua = case_when(
      survey == "NEUS" ~ num,
      (source == "DFO" & !(survey %in% c("gsl-s", "gmex"))) | survey %in% c("AI", "GOA", "WCANN", "WCTRI") ~ wgt_cpua,
      TRUE ~ num_cpua
    )
  ) %>%
  rename(abundance = num_cpua, region = continent)

# Add a continent/region column based on the survey
fishglob <- fishglob %>%
  mutate(region = case_when(
    # Western US surveys
    survey %in% c('EBS', 'AI', 'GOA', 'DFO-HS', 'DFO-QCS', 'DFO-SOG',
                  'DFO-WCHG', 'DFO-WCVI', 'WCANN', 'WCTRI') ~ 'Western_US',
    # Eastern US surveys
    survey %in% c('GSL-N', 'GSL-S', 'SCS', 'NEUS', 'SEUS', 'GMEX') ~ 'Eastern_US',
    # European surveys
    survey %in% c('BITS', 'EVHOE', 'FR-CGFS', 'IE-IGFS', 'NIGFS', 'Nor-BTS',
                  'NS-IBTS', 'PT-IBTS', 'ROCKALL', 'SP-ARSA', 'SP-NORTH',
                  'SP-PORC', 'SWC-IBTS') ~ 'Europe',
    TRUE ~ NA_character_  # Assign NA for any unmatched survey
  )) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth)

#-------------------------------------
# 3. Process biotime data sets
#-------------------------------------

Pelagic_Fish_Observations <- fread("Pelagic_Fish_Observations_1968-1999.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(),
    survey = "Pelagic_Fish_Observations", region = "Antartica") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit()

St_Croix_USVI_Fish_Assessment_and_Monitoring_Data <- fread("St_Croix_USVI_Fish_Assessment_and_Monitoring_Data.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(),
    survey = "St_Croix_USVI_Fish_Assessment_and_Monitoring_Data", region = "Caribbean") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit()

Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean <- fread("Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,5] %>%                           
      as.numeric(), 
    survey = "Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean", region = "Europe") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # qeustionable, is shorter than 10 years, so might just get filtered out 


Previous_fisheries_REVIZEE_Program_285 <- fread("Previous_fisheries_REVIZEE_Program_285.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # 

Previous_fisheries_REVIZEE_Program_284 <- fread("Previous_fisheries_REVIZEE_Program_284.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # 

Previous_fisheries_REVIZEE_Program_135 <- fread("Previous_fisheries_REVIZEE_Program_135.csv") %>%
  select(-ABUNDANCE) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program_135", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = BIOMAS, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # 

Alcatrazes_monitoring_program_660 <- fread("Alcatrazes_monitoring_program_660.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,2] %>%                           
      as.numeric(), 
    survey = "Alcatrazes_monitoring_program_660", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # has no clear depth, think its second column

ECNASAP_East_Coast_North_America_Strategic_Assessment <- fread("ECNASAP_East_Coast_North_America_Strategic_Assessment.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,5] %>%                           
      as.numeric(), 
    survey = "ECNASAP_East_Coast_North_America_Strategic_Assessment", region = "Eastern_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # 

MARMAP_Chevron_Trap_Survey_1990_2009 <- fread("MARMAP_Chevron_Trap_Survey_1990_2009.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "MARMAP_Chevron_Trap_Survey_1990_2009", region = "Eastern_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() # trap surveys, not trawl or dive

Pacific_Shrimp_Trawl_Survey <- fread("Pacific_Shrimp_Trawl_Survey(OBIS_Canada).csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Pacific_Shrimp_Trawl_Survey(OBIS_Canada)", region = "Western_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() 

CRED_Fish_Surveys_in_the_Pacific_Ocean <- fread("CRED_Fish_Surveys_in_the_Pacific_Ocean.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "CRED_Fish_Surveys_in_the_Pacific_Ocean", region = "Pacific") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() 

CSIRO_Marine_Data_Warehouse <- fread("CSIRO_Marine_Data_Warehouse.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "CSIRO_Marine_Data_Warehouse", region = "Oceania") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR) %>%
  select(region, survey, haul_id, latitude, longitude, year, accepted_name, abundance, depth) %>%
  na.omit() 

#-------------------------------------
# 4.Join individual data sets together
#-------------------------------------

# Join the data sets together
data <- bind_rows(fishglob, RLS, Pelagic_Fish_Observations, St_Croix_USVI_Fish_Assessment_and_Monitoring_Data, 
                  Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean, Previous_fisheries_REVIZEE_Program_285,
                  Previous_fisheries_REVIZEE_Program_284, Previous_fisheries_REVIZEE_Program_135, 
                  Alcatrazes_monitoring_program_660, ECNASAP_East_Coast_North_America_Strategic_Assessment, 
                  MARMAP_Chevron_Trap_Survey_1990_2009, Pacific_Shrimp_Trawl_Survey, CRED_Fish_Surveys_in_the_Pacific_Ocean,
                  CSIRO_Marine_Data_Warehouse, CSIRO_Marine_Data_Warehouse)

# Correct species naming format
data$accepted_name <- gsub("(^[a-z])", "\\U\\1", data$accepted_name, perl = TRUE) 

# Remove non species (e.g. genus) level data and filter out hauls deeper than 200m, and 0 or lower
data <- data %>%
  filter(str_detect(accepted_name, "\\s"), !str_detect(accepted_name, "spp\\.")) %>%
  na.omit() %>%
  filter(depth <= 200 & depth > 0) %>%
  select(-depth)

# Write full data set to a .csv file
write.csv(data, 'Data.csv', row.names = F)



################################################################################################################

Pelagic_Fish_Observations <- fread("Pelagic_Fish_Observations_1968-1999.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(),
    survey = "Pelagic_Fish_Observations", region = "Antartica") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit()

St_Croix_USVI_Fish_Assessment_and_Monitoring_Data <- fread("St_Croix_USVI_Fish_Assessment_and_Monitoring_Data.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(),
    survey = "St_Croix_USVI_Fish_Assessment_and_Monitoring_Data", region = "Caribbean") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit()

Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean <- fread("Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,5] %>%                           
      as.numeric(), 
    survey = "Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean", region = "Europe") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # qeustionable, is shorter than 10 years, so might just get filtered out 


Previous_fisheries_REVIZEE_Program_285 <- fread("Previous_fisheries_REVIZEE_Program_285.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # 

Previous_fisheries_REVIZEE_Program_284 <- fread("Previous_fisheries_REVIZEE_Program_284.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # 

Previous_fisheries_REVIZEE_Program_135 <- fread("Previous_fisheries_REVIZEE_Program_135.csv") %>%
  select(-ABUNDANCE) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Previous_fisheries_REVIZEE_Program_135", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = BIOMAS, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # 

Alcatrazes_monitoring_program_660 <- fread("Alcatrazes_monitoring_program_660.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,2] %>%                           
      as.numeric(), 
    survey = "Alcatrazes_monitoring_program_660", region = "South_America") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # has no clear depth, think its second column

ECNASAP_East_Coast_North_America_Strategic_Assessment <- fread("ECNASAP_East_Coast_North_America_Strategic_Assessment.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,5] %>%                           
      as.numeric(), 
    survey = "ECNASAP_East_Coast_North_America_Strategic_Assessment", region = "Eastern_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # 

MARMAP_Chevron_Trap_Survey_1990_2009 <- fread("MARMAP_Chevron_Trap_Survey_1990_2009.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "MARMAP_Chevron_Trap_Survey_1990_2009", region = "Eastern_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() # trap surveys, not trawl or dive

Pacific_Shrimp_Trawl_Survey <- fread("Pacific_Shrimp_Trawl_Survey(OBIS_Canada).csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "Pacific_Shrimp_Trawl_Survey(OBIS_Canada)", region = "Western_US") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() 

CRED_Fish_Surveys_in_the_Pacific_Ocean <- fread("CRED_Fish_Surveys_in_the_Pacific_Ocean.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "CRED_Fish_Surveys_in_the_Pacific_Ocean", region = "Pacific") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() 

CSIRO_Marine_Data_Warehouse <- fread("CSIRO_Marine_Data_Warehouse.csv") %>%
  select(-BIOMAS) %>% 
  mutate(
    depth = SAMPLE_DESC %>%
      str_split("_", simplify = TRUE) %>% 
      .[,6] %>%                           
      as.numeric(), 
    survey = "CSIRO_Marine_Data_Warehouse", region = "Oceania") %>%
  rename(haul_id = SAMPLE_DESC, accepted_name = valid_name, abundance = ABUNDANCE, latitude = LATITUDE, longitude = LONGITUDE, year = YEAR, month = MONTH, day = DAY) %>%
  select(region, survey, haul_id, latitude, longitude, year, month, day, accepted_name, abundance, depth) %>%
  na.omit() 

library(data.table)
library(dplyr)
library(lubridate)

data <- bind_rows(
  Pelagic_Fish_Observations,
  St_Croix_USVI_Fish_Assessment_and_Monitoring_Data,
  Marine_Fish_Underwater_Surveys_in_The_Israeli_Mediterranean,
  Previous_fisheries_REVIZEE_Program_285,
  Previous_fisheries_REVIZEE_Program_284,
  Previous_fisheries_REVIZEE_Program_135,
  Alcatrazes_monitoring_program_660,
  ECNASAP_East_Coast_North_America_Strategic_Assessment,
  MARMAP_Chevron_Trap_Survey_1990_2009,
  Pacific_Shrimp_Trawl_Survey,
  CRED_Fish_Surveys_in_the_Pacific_Ocean,
  CSIRO_Marine_Data_Warehouse)

data3 <- data %>%
  # create a proper Date column first
  mutate(survey_date = make_date(year, month, day)) %>%
  
  # group by survey (and region if you like)
  group_by(region, survey) %>%
  
  # compute start/end as Dates (not the formatted strings)
  summarize(
    survey_start_date = min(survey_date, na.rm = TRUE),
    survey_end_date   = max(survey_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # now compute decimal-year duration
  mutate(
    duration_years = as.numeric(interval(survey_start_date, survey_end_date) / years(1))
  ) %>%
  select(-region)


table <- data %>%
  group_by(survey) %>%
  summarise(min_month = min(month), max_month = max(month), mean_month = mean(month), mean_lat = mean(latitude), min_depth = min(depth), max_depth = max(depth), mean_depth = mean(depth)) 

table_done <- table %>%
  left_join(data3, by = ('survey'))


write.csv(table_done, 'table_done.csv', row.names = F)
