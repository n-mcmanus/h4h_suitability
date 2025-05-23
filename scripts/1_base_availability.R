# This script generates a raster of baseline H4H activity availability
# based on criteria such as LULC, slope, etc.

## Read in packages
library(tidyverse)            # always
library(here)                 # easier file paths
library(terra)                # GIS functions
library(sf)                   # vector functions
library(furrr)                

## Read in ROI
ssa_v <- vect(here("data/ROI/ssa_ne.shp"))

# --------------- 1. Isolate LULC classes suitable for H4H ---------------------

# This section reclassifies GLAD LULC data into a binary H4H availability layer,
# aggregates and resamples to match the carbon data, masks out any additional
# impervious surfaces, and then finally merges all tiles and masks them to the ROI


## 1a. Reclassify GLAD & mask --------------------------
## For each GLAD tile, reclassify to only keep land covers associated with H4H.
## Then mask out impervious surfaces

## Get file names and paths for each raster
tiles_df <- data.frame(
  ## raw GLAD LULC files
  "glad" = list.files(
    here("data/GLAD/tiles/raw"),
    pattern = "*.tif$",
    full.names = TRUE
  ),
  ## impervious surface data
  "imp_surface" = list.files(
    here("data/Global Impervious Surfaces products/reclassified"),
    pattern = "*.tif$",
    full.names = TRUE
  )
)

## match CRS of ROI and GLAD rasters
glad_crs <- crs(rast(tiles_df$glad[1]))
ssa_v <- project(ssa_v, y = glad_crs)

## checking if CRS needs to change. Luckily no!
naturebase_r <- rast(here("data/naturebase/africa/grs_agc_tco2eha_ssa.tif"))
imp_r <- rast(tiles_df$imp_surface[1])
# glad_crs == crs(naturebase_r)
# [1] TRUE
# glad_crs == crs(imp_r)
# [1] TRUE

## binary reclassification matrix
reclass_m <- matrix(c(
  -Inf, 1, 0,     # remove true desert
  2, 18, 1,       # keep semi-arid
  19, 24, 1,      # keep dense short veg
  25, 27, 1,      # keep trees under 5m
  28, 48, 0,      # remove taller trees
  100, 101, 0,    # remove salt pans
  102, 124, 1,    # keep short veg wetlands
  125, 127, 1,    # keep <5m trees wetlands
  102, 124, 0,    # remove short veg wetlands
  125, 127, 0,    # remove <5m trees wetlands
  128, 148, 0,    # remove tall tree wetlands
  200, 207, 0,    # remove open surface water
  241, 241, 0,    # remove snow/ice
  244, 244, 1,    # keep cropland for now?
  250, Inf, 0     # remove built-up and ocean
), ncol = 3, byrow = TRUE)

## Create output folder
dir_30m <- here("data/GLAD/tiles/reclassified/30m")
if (!dir.exists(dir_30m)) dir.create(dir_30m)


## Testing if for loop is faster
for (i in 1:nrow(tiles_df)) {
  glad <- tiles_df$glad[i]
  imp_surface <- tiles_df$imp_surface[i]
  
  print(paste0("Working on tile: ", basename(glad), " (", i, " out of ", nrow(tiles_df), ")"))
  
  ## assign file name
  fname <- paste0(basename(tools::file_path_sans_ext(glad)),
                  "_rcl_binary_30m.tif")
  
  ## Only generate if tile not already processed
  ## NOTE: this function can take several hours so may get interrupted
  if (file.exists(file.path(dir_30m, fname))) {
    next
  } else {
    ## read in LULC raster
    glad_r <- rast(glad)
    
    ## reclassify
    glad_rcl <- classify(glad_r, reclass_m, include.lowest = TRUE, right = NA)
    
    ## read in and resample impervious surface tile to match GLAD
    imp_r <- rast(imp_surface) %>% 
      resample(., y = glad_rcl, method = "near")
    
    ## now mask out impervious surfaces
    glad_mask <- mask(glad_rcl, imp_r, maskvalue = 1, updatevalue = 0)
    
    ## save
    writeRaster(glad_mask, 
                file.path(dir_30m, fname),
                datatype = "INT1U",
                overwrite = TRUE) 
  } 
}

# ## Function for reclassifying
# reclass_30m_glad <- function(glad, imp_surface) {
#   ## read in LULC raster
#   glad_r <- rast(glad)
#  
#   ## assign file name
#   fname <- paste0(basename(tools::file_path_sans_ext(glad)),
#                   "_rcl_binary_30m.tif")
#   
#   ## only generate if tile not already processed
#   ## (this function can take several hours and may get interrupted)
#   if (!file.exists(file.path(dir_30m, fname))) {
#     ## reclassify
#     glad_rcl <- classify(glad_r, reclass_m, include.lowest = TRUE, right = NA)
#     
#     ## read in and resample impervious surface tile to match GLAD
#     imp_r <- rast(imp_surface) %>% 
#       resample(., y = glad_rcl, method = "near")
#     
#     ## now mask out impervious surfaces
#     glad_mask <- mask(glad_rcl, imp_r, maskvalue = 1, updatevalue = 0)
#     
#     ## save
#     writeRaster(glad_mask, 
#                 file.path(dir_30m, fname),
#                 datatype = "INT1U",
#                 overwrite = TRUE) 
#   } 
# }

## run fxn for all tiles
# plan(multisession, workers = 8)
# pmap(.x = tiles_df, .f = reclass_30m_glad, .progress = TRUE)
# plan(sequential)



## 1b. Aggregate --------------------------
## For reclassified tiles, match naturebase resolution (~900m) and save

## Output directory for 1km tiles
dir_1km <- here("data/GLAD/tiles/reclassified/1km")
if (!dir.exists(dir_1km)) dir.create(dir_1km)

## Get list of reclassified 30m tiles (previous step)
tiles_rcl_list <- list.files(dir_30m, pattern = "*.tif$", full.names = TRUE)

## Fxn to aggregate and resample to match NatureBase (~900m)
resample_1km_glad <- function(tile, fun) {
  ## read in raster
  r <- rast(tile)
  
  ## match resolution & extent of NatureBase
  ## first aggregate to get close, then resample to match exactly
  factor <- ceiling(res(naturebase_r)[1] / res(r)[1])
  r_agg <- terra::aggregate(r, fact = factor, fun = fun) 
  r_resample <- resample(r_agg, naturebase_r, method = "near") #categorical data
  
  ## assign file name
  fname <- basename(tile) %>% 
    gsub("30m", paste0("1km"), x=.)
  
  ## save
  writeRaster(r_resample, 
              file.path(dir_1km, fname),
              datatype = "INT1U",
              overwrite = TRUE)
}

## run fxn 
map(.x = tiles_rcl_list, 
    .f = resample_1km_glad, 
    fun = "modal", # aggregation method (categorical)
    .progress = TRUE)


## 1c. Merge and mask to SSA ------------------------

## Get all tiles
fnames <- list.files(dir_1km, pattern = "*.tif$", full.names = TRUE)

## stack them & merge
tiles_stack <- lapply(fnames, rast)
tiles_merged <- do.call(terra::merge, tiles_stack)

## Crop/mask and export
tiles_merged_ssa <- crop(tiles_merged, ssa_v, mask = TRUE)

writeRaster(tiles_merged_ssa, 
            here("data/GLAD/glad_2020_ssa_binary_1km.tif"),
            overwrite = TRUE)


# ---------------- 2. Remove steep terrains ------------------------
# Here, we mask out any areas with a slope over 30% 
# Data comes from USGS HMDA global dataset (Africa subset available)
# Can download here: https://www.sciencebase.gov/catalog/item/591f6d17e4b0ac16dbdde1cd

## Pathway to data
slope_fpath <- here("data/USGS_HMDA")

## List files
slope_fnames = list.files(file.path(slope_fpath, "raw"),
                          pattern = "*.tif$", 
                          full.names = T)

## Merge all rasters
slope_stack <- lapply(slope_fnames, rast)
slope_merged <- do.call(terra::merge, slope_stack)

## Save 
writeRaster(slope_merged, 
            file.path(slope_fpath, "merged.tif"),
            overwrite = TRUE)

## Read in merged 1km GLAD data
lulc_r <- rast(here("data/GLAD/glad_2020_ssa_binary_1km.tif"))

## Aggregate/resample slope data to match resolution
factor <- ceiling(res(lulc_r)[1] / res(slope_merged)[1])

slope_agg <- aggregate(slope_merged, 
                       fact = factor, 
                       fun = "mean") ##maybe Max? depends

slope_resamp <- resample(slope_agg, lulc_r, "bilinear")

## Reclassify to binary
## Anything over 30% slope is not suitable (NA)
slope_rcl <- classify(slope_resamp,
                      matrix(c(30, Inf, NA), ncol = 3, byrow = TRUE))

slope_ssa <- crop(slope_rcl, ssa_v, mask = TRUE)

writeRaster(slope_ssa,
            here("data/USGS_HMDA/slope_30_1km.tif"),
            datatype = "INT1U",
            overwrite = TRUE)

## Now use slope as mask for availability
avail_mask <- mask(lulc_r, slope_ssa, updatevalue = 0) %>% 
  mask(ssa_v)

writeRaster(avail_mask, 
            here("data/outputs/glad_2020_slope30_ssa_binary_1km.tif"),
            overwrite = TRUE)




# ---------------- 3. Relative suitability from carbon ------------------------




## order of operations
layers <- c(
  "grs_agc",  #avoided grassland conversion
  "grs_asc",  #avoided shrubland conversions
  "wet_ipm",   #improved peatland mgmt
  "grs_scg",  #increased soil carbon in grazing lands
  "grs_sfm",  #savannah fire management
  "grs_grr"  #grassland restoration
)

fpath <- here("data/naturebase/africa")
r_ssa[r_ssa == 0] <- NA ## ADD THIS IN BELOW SOMEWHERE

for (i in seq_along(layers)) {
  layer <- layers[i]
  
  r <- rast(file.path(fpath, paste0(layer, "_tco2eha_ssa_na.tif")))
  
  if (i == 1) {
    r_total <- r
  } else {
    r_mask <- mask(r, r_total, inverse = TRUE)
    r_total <- sum(r_total, r_mask, na.rm = TRUE)
  }
}


## Read in binary GLAD LULC (generated in previous section)
glad_r <- rast(here("data/outputs/glad_no_imp_slope30_2020_ssa_binary_1km.tif"))
glad_r[glad_r == 0] <- NA

carbon_masked <- mask(r_total, glad_r)

glad_r[glad_r == 1] <- 0
test <- cover(carbon_masked, glad_r)

writeRaster(test,
            here("data/h4h_availability_slope_imp_carbon_1km.tif"),
            overwrite = TRUE)
















# Here, we sum and normalize all the carbon data, then mask to the baseline
# H4H availability (output of previous section)

## Read in binary GLAD LULC (generated in previous section)
glad_r1 <- rast(here("data/GLAD/glad_2020_ssa_binary_1km.tif"))

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
            here("data/h4h_suitability_carbon_nowet_1km.tif"),
            overwrite = TRUE)


## difference
nowet <- rast(here("data/h4h_suitability_carbon_nowet_1km.tif"))
nowet <- ifel(!is.na(nowet), 1, 0)

wet <- rast(here("data/outputs/h4h_suitability_carbon_modal_1km_sumfirst.tif"))
wet <- ifel(!is.na(wet), 1, 0)
diff <- wet - nowet %>% mask(., ssa_v)
writeRaster(diff, here("data/outputs/wetland_areas_1km.tif"))

## viz wetland on top 
wet_diff <- wet + diff
writeRaster(wet_diff, here("data/outputs/wetland_baseline_1km.tif"))

### were the areas cut out with high suitability?
diff <- rast(here("data/outputs/wetland_areas_1km.tif"))
diff[diff == 0] <- NA
wet = rast(here("data/outputs/h4h_suitability_carbon_modal_1km_sumfirst.tif"))
wet_vals <- wet %>% mask(diff)
writeRaster(wet_vals, here("data/outputs/wetland_area_suitability_1km.tif"))

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





