# This script wrangles multiple data used in this analysis

## Read in packages
library(tidyverse)            # always
library(here)
library(terra)                # GIS functions
library(sf)                   # vector functions
library(rnaturalearth)        # get countries oulintes
library(rnaturalearthdata)    # ' '
library(httr)                 # download GLAD via GEE url

# ROI --------------------------------------------------------------------
# First generate shapefile defining region of interest to crop
# and mask all other data

## get continent from Natural Earth dataset
africa_sf <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(region_un == "Africa") %>%
  ## filter out small/island territories
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  filter(area_km2 > 5000,
         admin != "Madagascar",
         type != "Dependency")

## also get only sub saharan countries
ssa_sf <- africa_sf %>% 
  filter(region_wb == "Sub-Saharan Africa")


# NatureBase --------------------------------------------------------------
# This section loops through the desired NCS layers and crops/masks to
# Sub-Sahara or continental Africa, then exports
# Direct download from: https://app.naturebase.org/data (all pathways)


## Only considering these NCS pathways for carbon (for now)
layers <- c(
  "grs_agc",  #avoided grassland conversion
  "grs_asc",  #avoided shrubland conversions
  "grs_grr",  #grassland restoration
  "grs_scg",  #increased soil carbon in grazing lands
  "grs_sfm",  #savannah fire management
  "wet_ipm"   #improved peatland mgmt... maybe remove?
)

## assign directories
dir_in <- here("data/naturebase/all_pathways")
dir_out <- here("data/naturebase/africa")
if (!dir.exists(dir_out)) dir.create(dir_out)

## loop through each
for (lyr in layers) {
  fname <- paste0(lyr, "_tco2eha.tif")
  
  ## read in raster
  r <- rast(file.path(dir_in, fname))
  
  ## match vector to raster crs
  v <- africa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  v_ssa <- ssa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  ## crop & mask raster
  r_crop <- crop(r, v, mask = TRUE)
  r_ssa_crop <- crop(r, v_ssa, mask = TRUE)
  
  ## export
  writeRaster(r_crop, 
              file.path(dir_out, paste0(fname, "_tco2eha_afr.tif")),
              overwrite = TRUE)
  
  writeRaster(r_ssa_crop, 
              file.path(dir_out, paste0(fname, "_tco2eha_ssa.tif")), 
              overwrite = TRUE)
  
}


# GLAD --------------------------------------------------------------------
# This section downloads all the tiles of LULC data touching our ROI
# then merges them together and masks to ROI.
# Define the latitude and longitude ranges
# Direct download available via: https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/download.html

## All possible lats/longs for continent tiles
latitudes <- seq(-30, 40, by = 10)
longitudes <- seq(-30, 50, by = 10)  

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


## Create list of all possible file names for African continetn
africa_coords <- paste0(
  "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020/",
  coords_NS$lat_formatted,
  "_",
  coords_NS$lon_formatted,
  ".tif"
)

## Get list of all the 2020 tiles worldwide
url_2020 <- "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020.txt"
all_2020_tiles <- readLines(url_2020)

## Only return tiles in African continent
africa_tiles <- intersect(all_2020_tiles, africa_coords)

## Download all raster tiles
map(africa_tiles, function(url) {
  GET(url, 
      write_disk(file.path(here("data/GLAD/tiles"), basename(url)), 
                 overwrite = TRUE))
})

## Merge together and crop/mask to Africa and SSA
tiles <- list.files(here("data/GLAD/tiles/"),
                    pattern = ".tif",
                    full.names = TRUE)
rasters <- lapply(tiles, rast)

merged_r <- do.call(merge, rasters)

africa_v <- vect(africa_sf) %>% 
  project(., merged_r)
ssa_v <- vect(ssa_sf) %>% 
  project(., merged_r)

glad_crop <- crop(merged_r, africa_v, mask = TRUE)
writeRaster(glad_crop, here("data/GLAD/glad_africa_2020.tif"))

glad_ssa_crop <- crop(glad_crop, ssa_v, mask = TRUE)
writeRaster(glad_ssa_crop, here("data/GLAD/glad_ssa_2020.tif"))





