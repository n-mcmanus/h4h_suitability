# This script generates a raster of baseline H4H activity availability
# based on criteria such as LULC, slope, etc.

## Read in packages
library(tidyverse)            # always
library(here)                 # easier file paths
library(terra)                # GIS functions
library(sf)                   # vector functions
library(purrr)                

## Read in ROI
ssa_sf <- read_sf(here("data/ROI/ssa_ne.shp"))

# --------------- 1. Isolate LULC classes suitable for H4H ---------------------

# This section reclassifies GLAD LULC data into a binary H4H availability layer,
# aggregates and resamples to match the carbon data, then finally merges
# all tiles and masks them to the ROI


## 1a. Reclassify --------------------------
## For each GLAD tile, reclassify to only keep land covers associated with H4H.

## Get names and paths for each raster
tile_fnames <- list.files(here("data/GLAD/tiles/raw"),
                    pattern = "*.tif$",
                    full.names = TRUE)

## match CRS of ROI and GLAD rasters
glad_crs <- crs(rast(tile_fnames[1]))
ssa_v <- vect(ssa_sf) %>% 
  project(., y = glad_crs)

## checking if CRS needs to change. Luckily no!
naturebase_r <- rast(here("data/naturebase/africa/grs_agc_tco2eha_ssa.tif"))
# crs(glad_crs) == crs(naturebase_r)
# [1] TRUE

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

## Create output folder
dir_30m <- here("data/GLAD/tiles/reclassified/30m")
if (!dir.exists(dir_30m)) dir.create(dir_30m)

## Function for reclassifying
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
map(.x = tile_fnames, .f = reclass_30m_glad, .progress = TRUE)


## 1b. Aggregate & Resample --------------------------
## For reclassified tiles, match naturebase resolution (~900m) and save

## Output directory for 1km tiles
dir_1km <- here("data/GLAD/tiles/reclassified/1km")
if (!dir.exists(dir_1km)) dir.create(dir_1km)

## Get list of reclassified 30m tiles
tiles_rcl_list <- list.files(dir_30m, pattern = ".tif", full.names = TRUE)

## Aggregate and resample to match NatureBase (~1km)
resample_1km_glad <- function(tile, fun) {
  ## read in
  r <- rast(tile)
  
  ## match resolution & extent of NatureBase
  ## first aggregate to get close, then resample to match exactly
  factor <- ceiling(res(naturebase_r)[1] / res(r)[1])
  r_agg <- terra::aggregate(r,
                            fact = factor,
                            fun = fun) 
  r_resample <- resample(r_agg, naturebase_r, method = "near")
  
  ## create file name
  fname <- basename(tile) %>% 
    gsub("30m", paste0("1km"), x=.)
  
  ## save
  writeRaster(r_resample, 
              file.path(dir_1km, fname),
              overwrite = TRUE)
}

## run fxn 
map(.x = tiles_rcl_list, 
    .f = resample_1km_glad, 
    fun = "modal", # aggregation method
    .progress=TRUE)



## 1c. Merge and mask to SSA ------------------------

## Get all tiles
fnames <- list.files(dir_1km, 
                     pattern = "*.tif$", 
                     full.names = TRUE)

## stack them & merge
tiles_stack <- lapply(fnames, rast)
tiles_merged <- do.call(terra::merge, tiles_stack)

## Crop/mask and export
# crs(ssa_sf) == crs(merged_r)
# [1] TRUE
ssa_v <- vect(ssa_sf)

tiles_merged_ssa <- crop(tiles_merged, ssa_v, mask = TRUE)

writeRaster(tiles_merged_ssa, 
            here("data/GLAD/glad_2020_ssa_binary_1km.tif"),
            overwrite = TRUE)




# ---------------- 2. Relative suitability from carbon ------------------------

# Here, we sum and normalize all the carbon data, then mask to the baseline
# H4H availability (output of previous section)

## Read in binary GLAD LULC (generated in previous section)
glad_r <- rast(here("data/GLAD/glad_2020_ssa_binary_1km.tif"))

## Get list of all NatureBase carbon layers
## For now, focusing just on SSA 
carbon_fnames <- list.files(here("data/naturebase/africa/"),
                            pattern = "ssa.tif$",
                            full.names = TRUE)


## stack all the carbon rasters
carbon_list <- rast(carbon_fnames)

## add them together
carbon_sum <- app(carbon_list, sum, na.rm = TRUE)

## normalize the carbon raster
min <- minmax(carbon_sum)[1, ]
max <- minmax(carbon_sum)[2, ]

carbon_norm <- (carbon_sum - min) / (max - min)

## mask to LULC
carbon_masked <- mask(carbon_norm, glad_r, maskvalues = 0, updatevalue = NA)

## export
writeRaster(carbon_masked, 
            here("data/h4h_suitability_carbon_1km.tif"),
            overwrite = TRUE)


# ### IF NORMALIZING FIRST, THEN SUMMING:  -----------
# ## create list of each pathway raster
# carbon_list <- lapply(carbon_fnames, rast)
# 
# ## normalize all rasters in list
# normalize_rast <- function(r) {
#   min <- global(r, "min", na.rm = TRUE)[1, ]
#   max <- global(r, "max", na.rm = TRUE)[1, ]
#   (r - min) / (max - min)
# }
# 
# carbon_norm <- map(carbon_list, normalize_rast, .progress = TRUE)
# 
# ## convert from list to stack and add together
# carbon_stack <- rast(carbon_norm)
# 
# carbon_sum <- app(carbon_stack, sum, na.rm = TRUE)
# 
# ## mask to LULC
# carbon_masked <- mask(carbon_sum, glad_r, maskvalues = 0, updatevalue = NA)
# 
# ## export
# writeRaster(carbon_masked, 
#             here("data/h4h_suitability_carbon_max_1km_normfirst.tif"),
#             overwrite = TRUE)





