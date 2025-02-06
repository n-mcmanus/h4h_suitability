# This script wrangles and prepares the data used in this analysis
# Author: Nick McManus

## Read in packages
library(tidyverse)            # always
library(here)                 # easier file paths
library(terra)                # GIS functions
library(sf)                   # vector functions
library(rnaturalearth)        # get countries oulintes
library(rnaturalearthdata)    # ' '
library(httr)                 # download GLAD via GEE url
library(furrr)                # faster parallel processing

# ROI --------------------------------------------------------------------
# First generate shapefiles defining region of interest. 
# Used to crop and mask other data.
# NOTE: Originally looked at both continent and sub-Saharan Africa (SSA). 
# Moving forward, likely only considering SSA. 

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


# NatureBase --------------------------------------------------------------
# This section loops through the desired NCS layers and crops/masks to ROI
# Direct data download from: https://app.naturebase.org/data (all pathways)


## Only considering the following NCS pathways for carbon (for now...)
layers <- c(
  "grs_agc",  #avoided grassland conversion
  "grs_asc",  #avoided shrubland conversions
  "grs_grr",  #grassland restoration
  "grs_scg",  #increased soil carbon in grazing lands
  "grs_sfm",  #savannah fire management
  "wet_ipm"   #improved peatland mgmt
)

## assign directories
dir_in <- here("data/naturebase/all_pathways")
dir_out <- here("data/naturebase/africa")
if (!dir.exists(dir_out)) dir.create(dir_out)


## Loop through each raster
for (lyr in layers) {
  ## get file name
  fname <- paste0(lyr, "_tco2eha.tif")
  
  ## read in raster
  r <- rast(file.path(dir_in, fname))
  
  ## match vector to raster crs
  af_v <- africa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  ssa_v <- ssa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  ## crop & mask raster
  r_af <- crop(r, af_v, mask = TRUE)
  r_ssa <- crop(r, ssa_v, mask = TRUE)
  
  ## export
  writeRaster(r_af, 
              file.path(dir_out, 
                        paste0(tools::file_path_sans_ext(fname), "_afr.tif")),
              overwrite = TRUE)
  
  writeRaster(r_ssa, 
              file.path(dir_out, 
                        paste0(tools::file_path_sans_ext(fname), "_ssa.tif")), 
              overwrite = TRUE)
  
}


# GLAD --------------------------------------------------------------------
# This section downloads all the tiles of LULC data touching our ROI (SSA),
# then merges them together and masks to ROI.
# Direct download available via: https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/download.html

## Download rasters ----------------------
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
map(ssa_tiles, function(url) {
  GET(url, 
      write_disk(file.path(here("data/GLAD/tiles/raw"), basename(url)), 
                 overwrite = TRUE))
})



## Reclassify & Resample ----------------------
## For each GLAD tile, reclassify to only keep land covers associated with H4H.
## then change resolution (to match naturebase ~90m) and save that version as well

## Get names and paths for each raster
tiles <- list.files(here("data/GLAD/tiles/raw"),
                    pattern = ".tif",
                    full.names = TRUE)

## match CRS of ROI and GLAD rasters
glad_crs <- crs(rast(tiles[1]))
ssa_v <- vect(ssa_sf) %>% 
  project(., y = glad_crs)

## checking if CRS needs to change. Luckily no!
naturebase_r <- rast(here("data/naturebase/africa/grs_agc_tco2eha_ssa.tif"))
# crs(temp_r) == crs(naturebase_r)
# [1] TRUE

## Create output directories for reclassified rasters
dir_30m <- file.path(here("data/GLAD/tiles/reclassified/30m"))
if (!dir.exists(dir_30m)) dir.create(dir_30m)
dir_1km <- file.path(here("data/GLAD/tiles/reclassified/1km"))
if (!dir.exists(dir_1km)) dir.create(dir_1km)


## binary reclass matrix
reclass_m <- matrix(c(
  -Inf, 1, 0,     # remove true desert
  2, 18, 1,       # keep semi-arid
  19, 24, 1,      # keep dense short veg
  25, 27, 1,      # keep trees under 5m
  28, 48, 0,      # remove taller trees
  100, 101, 0,    # remove salt pans
  102, 124, 1,    # keep short veg wetlands???
  125, 127, 1,    # keep <5m trees wetlands???
  128, 148, 0,    # remove tall tree wetlands
  200, 207, 0,    # remove open surface water
  241, 241, 0,    # remove snow/ice
  244, 244, 1,    # keep cropland for now?
  250, Inf, 0     # remove built-up and ocean
), ncol = 3, byrow = TRUE)

## function for reclassifying
reclass_30m_glad <- function(tile) {
    ## read in raster
    r <- rast(tile)
    ## get file names
    fname <- paste0(basename(tools::file_path_sans_ext(tile)),
                    "_rcl_binary_30m.tif")
    
    ## reclassify
    r_rcl <- classify(r, reclass_m, include.lowest = TRUE, right = NA)
    ## save
    writeRaster(r_rcl, 
                file.path(dir_30m, fname),
                overwrite = TRUE) 
}

## run fxn for all tiles
map(tiles, reclass_30m_glad, .progress = TRUE)



## Get list of reclassified 30m tiles
tiles30m_list <- list.files(dir_30m,
                            pattern = ".tif",
                            full.names = TRUE)

## Aggregate and resample to match NatureBase (~1km)
## NOTE: Do we want to go with max or modal function? 
resample_1km_glad <- function(tile, fun) {
  ## read in
  r <- rast(tile)
  
  ## get file name
  fname <- basename(tile) %>% 
    # gsub("30m", "1km", x=.)
    gsub("30m", "1km_modal", x=.)
  
  ## match resolution of NatureBase
  factor <- ceiling(res(naturebase_r)[1] / res(r)[1])
  r_agg <- terra::aggregate(r,
                            fact = factor,
                            # fun = "max")
                            fun = "modal") ## if any part is 1, make whole thing 1
  r_resample <- resample(r_agg, naturebase_r, method = "near")
  
  ## save
  writeRaster(r_resample, 
              file.path(dir_1km, fname),
              overwrite = TRUE)
}

plan(multisession, workers = 6) 
map(tiles30m_list, resample_1km_glad, .progress = TRUE)
plan(sequential)

## Merge and mask to Africa ----------------------
## get all tiles
# fnames <- list.files(dir_1km, pattern = "*.tif$", full.names = TRUE)
fnames <- list.files(dir_1km, pattern = "*modal.tif", full.names = TRUE)

## stack them & merge
# stack <- lapply(fnames, rast)
# merged <- do.call(terra::merge, stack)

### Not sure why, terra::merge results in continuous layer despite 
### all inputs being binary. terra::vrt combines them better as binary
### Can dig deeper into difference here later...
### NOTE: may have been due to saving the agg rather than resampled raster in function...
### can re-run and test if the resampled version works better or if it just makes the whole thing take longer
merged_r <- vrt(fnames) 

## Crop/mask and export
# crs(ssa_sf) == crs(merged_r)
# [1] TRUE
ssa_v <- vect(ssa_sf)

merged_ssa_r <- crop(merged_r, ssa_v, mask = TRUE)

writeRaster(merged_ssa_r, 
            # here("data/GLAD/glad_2020_ssa_maxbinary_1km.tif"),
            here("data/GLAD/glad_2020_ssa_modalbinary_1km.tif"),
            overwrite = TRUE)
