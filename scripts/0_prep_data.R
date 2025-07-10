## Script: 0_prep_data.R
##
## Purpose: This script wrangles and prepares the data used for other scripts
## in this analysis. This includes downloading, reclassifying, resampling, and more.
##
## Last updated: 02 June 2025
##
## Author: Nick McManus
## Email: nmcmanus@conservation.org
##
## Notes: Subject to updates as more datasets are incorporated and suitability
## analysis is further developed.

## Set up ------------------------------------------------------------------
## Install package 'pacman' if needed
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  terra,              # GIS functions
  sf,                 # vector functions
  rnaturalearth,      # get country outlines
  rnaturalearthdata,  # " "
  httr,               # Helps download GLAD data via GEE URL
  furrr)              # loads both future and purrr packages


## ROI --------------------------------------------------------------------
## First generate shape files defining region of interest. 
## These are used to crop and mask other data.
##
## NOTES: Currently only evaluating sub-Saharan Africa (SSA) using data from
## Natural Earth. May change ROI with more specific polygon for SSA in future.

## Get continent from Natural Earth dataset
africa_sf <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  ## keep only countries in Africa
  filter(region_un == "Africa") %>%
  ## filter out small and/or island territories for this analysis
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%
  filter(area_km2 > 5000,
         admin != "Madagascar", #remove for now. maybe include later
         type != "Dependency")

## Get only SSA countries
ssa_sf <- africa_sf %>% 
  filter(region_wb == "Sub-Saharan Africa") 

## Save for later use
dir_out <- here("data/ROI")
if (!dir.exists(dir_out)) dir.create(dir_out) #creates directory locally if it doesn't exist

write_sf(ssa_sf, file.path(dir_out, "ssa_ne.shp"))


## NatureBase --------------------------------------------------------------
## This section loops through the desired NCS layers and crops/masks to ROI.
## Direct data download from: https://app.naturebase.org/data (all pathways)

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
dir_in <- here("data/naturebase/all_pathways") #where data is saved locally
dir_out <- here("data/naturebase/africa")      #desired output directory
if (!dir.exists(dir_out)) dir.create(dir_out)


## Loop through each global raster, crop & mask to ROI, then export
for (lyr in layers) {
  print(paste("Working on:", lyr))
  
  ## get file name
  fname <- paste0(lyr, "_tco2eha.tif")
  
  ## read in raster
  r <- rast(file.path(dir_in, fname))
  
  ## match vector to raster crs
  ssa_v <- ssa_sf %>% 
    vect() %>% 
    project(., crs(r))
  
  ## crop & mask raster
  r_ssa <- crop(r, ssa_v, mask = TRUE)
  
  ## export
  writeRaster(r_ssa, 
              file.path(dir_out, 
                        paste0(tools::file_path_sans_ext(fname), "_ssa.tif")),  #amend original file name
              overwrite = TRUE)
}


## GLAD --------------------------------------------------------------------
## This section downloads all the tiles of LULC data intersecting our ROI,
## then merges them together and masks to ROI.
## Direct download available via: https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/download.html

## Create output directory
dir_out <- here("data/GLAD/tiles/raw")
if (!dir.exists(dir_out)) dir.create(dir_out)

## All possible lats/longs for SSA tiles
latitudes <- seq(-30, 30, by = 10)
longitudes <- seq(-20, 40, by = 10)  

## Generate the combinations of latitude and longitude
coords <- expand.grid(lat = latitudes, lon = longitudes)

## Format the lon/lat to match file naming convention
coords_formatted <- coords %>%
  mutate(
    lon = case_when(lon < 0 ~ sprintf("%03dW", abs(lon)), # format neg lons with W
                    .default = sprintf("%03dE", lon)),    # format positive lons with E
    lat = case_when(lat < 0 ~ sprintf("%02dS", abs(lat)),
                    .default = sprintf("%02dN", lat))
  )

## Create list of all possible file names for SSA
ssa_coords <- paste0(
  "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020/",
  coords_formatted$lat,
  "_",
  coords_formatted$lon,
  ".tif"
)

## Get list of all the 2020 GLAD tiles worldwide
url_2020 <- "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2020.txt"
all_2020_tiles <- readLines(url_2020)

## Only return tiles in SSA
ssa_tiles <- intersect(all_2020_tiles, ssa_coords)

## Download all raster tiles
purrr::map(ssa_tiles, function(url) {
  GET(url, write_disk(file.path(dir_out, basename(url)), overwrite = TRUE))
})


## Impervious Surfaces --------------------------------------------------------
## Process tiles from global impervious surface map.These will be used to
## remove roads and other surfaces not classified as urban in GLAD LULC.
## Data can be directly downloaded here: https://zenodo.org/records/3505079

## Set filepath for locally saved impervious surface data
fpath <- here("data/Global Impervious Surfaces products")

## Set output directory
dir_out <- file.path(fpath, "reclassified")
if (!dir.exists(dir_out)) dir.create(dir_out)

## Similar process as GLAD, but with slightly different formatting
latitudes <- seq(-30, 30, by = 10)
longitudes <- seq(-20, 40, by = 10)  

## Generate the combinations of latitude and longitude
coords <- expand.grid(lat = latitudes, lon = longitudes)

## Create df to match file naming conventions
coords_df <- coords %>%
  mutate(
    ## get format used for imp surface data
    lon_lat = paste0(case_when(lon < 0 ~ sprintf("W%d", abs(lon)), # format neg lons with W
                               .default = sprintf("E%d", lon)), 
                     case_when(lat < 0 ~ sprintf("S%d", abs(lat)), 
                               .default = sprintf("N%d", lat))),
    ## list of potential file names
    fpath = file.path(fpath, "raw", paste0("ImperviousMap_", lon_lat, ".tif")),
    ## add formatting for GLAD data
    lat_lon_glad = paste0(case_when(lat < 0 ~ sprintf("%02dS", abs(lat)),
                                    .default = sprintf("%02dN", lat)),
                          "_",
                          case_when(lon < 0 ~ sprintf("%03dW", abs(lon)), # format neg lons with W
                                    .default = sprintf("%03dE", lon)))
  ) %>% 
  dplyr::select(c(fpath, lat_lon_glad))


## List all files in dataset
imp_files <- data.frame("fpath" = list.files(
  file.path(fpath, "raw"),
  pattern = "*.tif$",
  full.names = TRUE
))

## Only keep those within SSA
imp_tiles_df <- semi_join(coords_df, imp_files, by = "fpath")

## Binary reclassification matrix
reclass_m <- matrix(c(
  1, 0,   #1 becomes 0
  2, 1    #2 becomes 1
), ncol = 2, byrow = TRUE)

## Reclassify all rasters
plan(multisession, workers = 8) #use parallel processing if possible
purrr::pmap(imp_tiles_df, .progress = TRUE, function(fpath, lat_lon_glad) {
  ## read in raster and reclassify 1s to 0s
  r <- rast(fpath)
  r_rcl <- classify(r, reclass_m)
  
  ## generate filename and save
  fname <- paste0("ImperviousMap_", lat_lon_glad, "_rcl.tif")
  writeRaster(r_rcl, 
              file.path(dir_out, fname), 
              datatype = "INT1U",  #reduce file size
              overwrite = TRUE)
})

## return to sequential processing
plan(sequential) 

