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
library(furrr)                # loads future and purrr


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


# Impervious Surfaces --------------------------------------------------------
# Process tiles from impervious surface map.
# Can download here: https://zenodo.org/records/3505079
# NOTE: May be incorporated into other script, but for now reclassifying here to
# make next steps faster and comparison easier


## set filepath for impervious surface data
imp_fpath <- here("data/Global Impervious Surfaces products")

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
    fpath = file.path(imp_fpath, "raw", paste0("ImperviousMap_", lon_lat, ".tif")),
    ## add formatting for GLAD data
    lat_lon_glad = paste0(case_when(lat < 0 ~ sprintf("%02dS", abs(lat)),
                                    .default = sprintf("%02dN", lat)),
                          "_",
                          case_when(lon < 0 ~ sprintf("%03dW", abs(lon)), # format neg lons with W
                                    .default = sprintf("%03dE", lon)))
  ) %>% 
  dplyr::select(c(fpath, lat_lon_glad))


## list all files in dataset
imp_files <- data.frame("fpath" = list.files(
  file.path(imp_fpath, "raw"),
  pattern = "*.tif$",
  full.names = TRUE
))

## only keep those within SSA
imp_tiles_df <- semi_join(coords_df, imp_files, by = "fpath")


## set output directory
out_dir <- file.path(imp_fpath, "reclassified")
if (!dir.exists(out_dir)) dir.create(out_dir)

## binary reclassification matrix
reclass_m <- matrix(c(
  1, 0,   #1 becomes 0
  2, 1    #2 becomes 1
), ncol = 2, byrow = TRUE)


## Reclassify all rasters
plan(multisession, workers = 8)
purrr::pmap(imp_tiles_df, .progress = TRUE, function(fpath, lat_lon_glad) {
  ## read in raster and reclassify 1s to 0s
  r <- rast(fpath)
  r_rcl <- classify(r, reclass_m)
  
  ## generate filename and save
  fname <- paste0("ImperviousMap_", lat_lon_glad, "_rcl.tif")
  writeRaster(r_rcl, file.path(out_dir, fname), datatype = "INT1U")
})
plan(sequential)
