#===============================================================================
# Range Data Extraction from GBIF
#===============================================================================

#-------------------------------------
# 1. Load required packages
#-------------------------------------

library(dplyr)      
library(data.table)  
library(rgbif)       
library(stringr)     

#-------------------------------------
# 2. GBIF credentials
#-------------------------------------

options(
  gbif_user  = "ben.mangnall",
  gbif_email = "vr21539@bristol.ac.uk",
  gbif_pwd   = "123456"
)

#-------------------------------------
# 3. Read data set
#-------------------------------------

input_csv <- 'Data.csv'
master_df <- fread(input_csv)

#-------------------------------------
# 4. Split by hemisphere and extract species lists
#-------------------------------------

north_df          <- filter(master_df, latitude >  0)
south_df          <- filter(master_df, latitude <  0)
species_northern  <- unique(north_df$accepted_name)
species_southern  <- unique(south_df$accepted_name)

#-------------------------------------
# 5. Function to retrieve and save GBIF occurrences
#-------------------------------------

fetch_gbif_by_hemisphere <- function(species_vec, hemi = c("northern", "southern"), output_file) {
  hemi <- match.arg(hemi)
  
  # 5.1 Resolve species names to GBIF taxon keys
  keys <- name_backbone_checklist(species_vec) %>%
    filter(matchType == "EXACT") %>%
    pull(usageKey)
  
  # 5.2 Request download from GBIF 
  download_id <- occ_download(
    pred_in("taxonKey", keys),
    pred("hasCoordinate", TRUE),
    format = "SIMPLE_CSV"
  )
  
  message(sprintf("Submitted GBIF download (%s) for %s hemisphere", download_id, hemi))
  
  # 5.3 Wait for download to complete 
  occ_download_wait(download_id)
  
  # 5.4 Import and filter
  temp_file <- occ_download_get(download_id, overwrite = TRUE)
  gbif_df   <- fread(temp_file, sep = "\t", fill = TRUE, quote = "") %>%
    select(species, decimalLatitude) %>%
    filter(
      species %in% species_vec,
      if (hemi == "northern") decimalLatitude >  0 else decimalLatitude < 0)
  
  # 5.5 Save output
  fwrite(gbif_df, file = output_file)
  message(sprintf("Wrote %s records to %s", nrow(gbif_df), output_file))
  
  # 5.6 Cleanup temporary files
  file.remove(temp_file)
  return(invisible(gbif_df))
}

#-------------------------------------
# 6. Execute for each hemisphere
#-------------------------------------

fetch_gbif_by_hemisphere(
  species_vec = species_northern,
  hemi        = "northern",
  output_file = "gbif_data_northern.csv"
)

fetch_gbif_by_hemisphere(
  species_vec = species_southern,
  hemi        = "southern",
  output_file = "gbif_data_southern.csv"
)

# End of script