# This script wrangles multiple data used in this analysis

## Read in packages
library(tidyverse)            # always
library(here)
library(terra)                # GIS functions
library(sf)                   # vector functions
library(rnaturalearth)        # get countries oulintes
library(rnaturalearthdata)    # ' '

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















