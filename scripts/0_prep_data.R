# This script wrangles and prepares the data used in this analysis
# Author: Nick McManus

## Read in packages
library(tidyverse)            # always
library(here)                 # easier file paths
library(terra)                # GIS functions
library(sf)                   # vector functions
library(rnaturalearth)        # get countries oulintes
library(rnaturalearthdata)    # " "
library(httr)                 # download GLAD via GEE url
library(furrr)                # faster parallel processing


# ROI --------------------------------------------------------------------

# First generate shapefiles defining region of interest. 
# These are used to crop and mask other data.
# **NOTE: Originally looked at both continent and sub-Saharan Africa (SSA). 
# **Moving forward, likely only considering SSA. 

## Get continent from Natural Earth dataset
africa_sf <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(region_un == "Africa") %>%
  ## filter out small and/or island territories
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  filter(area_km2 > 5000,
         admin != "Madagascar",
         type != "Dependency")

## get only sub saharan countries
ssa_sf <- africa_sf %>% 
  filter(region_wb == "Sub-Saharan Africa") 


## Save for later use
dir_out <- here("data/ROI")
if (!dir.exists(dir_out)) dir.create(dir_out)

write_sf(africa_sf, file.path(dir_out, "africa_ne.shp"))
write_sf(ssa_sf, file.path(dir_out, "ssa_ne.shp"))

# NatureBase --------------------------------------------------------------

# This section loops through the desired NCS layers and crops/masks to ROI.
# Direct data download from: https://app.naturebase.org/data (all pathways)

## Define which NCS pathways to consider
layers <- c(
  "grs_agc",  #avoided grassland conversion
  "grs_asc",  #avoided shrubland conversions
  "grs_grr",  #grassland restoration
  "grs_scg",  #increased soil carbon in grazing lands
  "grs_sfm",  #savannah fire management
  "wet_ipm"   #improved peatland mgmt
)

## Assign directories
dir_in <- here("data/naturebase/all_pathways")
dir_out <- here("data/naturebase/africa")
if (!dir.exists(dir_out)) dir.create(dir_out)


## For each raster, crop & mask to ROI and export
for (lyr in layers) {
  ## get file name
  fname <- paste0(lyr, "_tco2eha.tif")
  
  ## read in raster
  r <- rast(file.path(dir_in, fname))
  
  # ## match vector to raster crs
  # af_v <- africa_sf %>% 
  #   vect() %>% 
  #   project(., crs(r))
  
  ssa_v <- ssa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  ## crop & mask raster
  # r_af <- crop(r, af_v, mask = TRUE)
  r_ssa <- crop(r, ssa_v, mask = TRUE)
  
  ## export
  # writeRaster(r_af, 
  #             file.path(dir_out, 
  #                       paste0(tools::file_path_sans_ext(fname), "_afr.tif")),
  #             overwrite = TRUE)
  
  writeRaster(r_ssa, 
              file.path(dir_out, 
                        paste0(tools::file_path_sans_ext(fname), "_ssa.tif")), 
              overwrite = TRUE)
  
}


# GLAD --------------------------------------------------------------------

# This section downloads all the tiles of LULC data touching our ROI,
# then merges them together and masks to ROI.
# Direct download available via: https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/download.html

## All possible lats/longs for SSA tiles
latitudes <- seq(-30, 30, by = 10)
longitudes <- seq(-20, 40, by = 10)  

## Generate the combinations of latitude and longitude
coords <- expand.grid(lat = latitudes, lon = longitudes)

## Format the lons/lats to match file naming convention
coords_NS <- coords %>%
  mutate(
    lon_formatted = case_when(lon < 0 ~ sprintf("%03dW", abs(lon)), # format neg lons with W
                              .default = sprintf("%03dE", lon)),    # format positive lons with E
    lat_formatted = case_when(lat < 0 ~ sprintf("%02dS", abs(lat)), 
                              .default = sprintf("%02dN", lat))
    )


## Create list of all possible file names for SSA
ssa_coords <- paste0(
  "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020/",
  coords_NS$lat_formatted,
  "_",
  coords_NS$lon_formatted,
  ".tif"
)

## Get list of all the 2020 tiles worldwide
url_2020 <- "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020.txt"
all_2020_tiles <- readLines(url_2020)

## Only return tiles in SSA
ssa_tiles <- intersect(all_2020_tiles, ssa_coords)

## Download all raster tiles
purrr::map(ssa_tiles, function(url) {
  GET(url, 
      write_disk(file.path(here("data/GLAD/tiles/raw"), basename(url)), 
                 overwrite = TRUE))
})


