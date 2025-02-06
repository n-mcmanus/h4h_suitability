# This script generates a raster of baseline H4H activity availability
# based on LULC

## Read in packages
library(tidyverse)            # always
library(here)                 # easier file paths
library(terra)                # GIS functions
library(sf)                   # vector functions
library(purrr)



# Can go about this two different ways:
# 1) Normalize each carbon layer, THEN add them together
# 2) Add all carbon together, THEN normalize
# Second method may better capture relative difference between areas....




## Read in binary GLAD LULC (generated in 0_prep_data.R)
glad_r <- rast(here("data/GLAD/glad_2020_ssa_modalbinary_1km.tif"))

## Get list of all NatureBase carbon layers
## For now, focusing just on SSA 
carbon_fnames <- list.files(here("data/naturebase/africa/"),
                            pattern = "ssa.tif",
                            full.names = TRUE)

carbon_list <- lapply(carbon_fnames, rast)

### NOTE!: need to play around with when to resample, seems to affect wetlands the most!
### when reading in, or at end. Difference? 

# carbon_list <- map(carbon_list, function(r) {
#   r_resamp <- resample(r, glad_r, method = "bilinear")
# })
# 
# carbon_list <- rast(carbon_list)
# min <- minmax(carbon_list)[1, ]
# max <- minmax(carbon_list)[2, ]
# 
# carbon_norm <- carbon_list - min / (max - min)
# 
# 
# c <- lapply(carbon_list, rast)
# Function to normalize a raster


normalize_rast <- function(r) {
  min <- global(r, "min", na.rm = TRUE)[1, ]
  max <- global(r, "max", na.rm = TRUE)[1, ]
  (r - min) / (max - min)
}

# Normalize each raster
carbon_norm <- map(carbon_list, normalize_rast, .progress = TRUE)

## convert from list to stack
carbon_stack <- rast(carbon_norm)

## add together
carbon_sum <- app(carbon_stack, sum, na.rm = TRUE)

## Mask carbon to suitable areas
# crs(carbon_sum) == crs(glad_r)
# [1] TRUE
# res(carbon_sum) == res(glad_r)
# [1] FALSE FALSE

carbon_masked <- resample(carbon_sum, glad_r, method = "bilinear") %>% 
  mask(., glad_r, maskvalues = 0, updatevalue = NA)

## export
writeRaster(carbon_masked, 
            here("data/h4h_suitability_carbon_modal_1km.tif"),
            overwrite = TRUE)




### 2) add together THEN normalize


## Read in binary GLAD LULC (generated in 0_prep_data.R)
glad_r <- rast(here("data/GLAD/glad_2020_ssa_modalbinary_1km.tif"))

## Get list of all NatureBase carbon layers
## For now, focusing just on SSA 
carbon_fnames <- list.files(here("data/naturebase/africa/"),
                            pattern = "ssa.tif",
                            full.names = TRUE)

carbon_list <- rast(carbon_fnames)

carbon_sum <- app(carbon_list, sum, na.rm = TRUE)

min <- minmax(carbon_sum)[1, ]
max <- minmax(carbon_sum)[2, ]

carbon_norm <- (carbon_sum - min) / (max - min)


## Mask carbon to suitable areas
# crs(carbon_sum) == crs(glad_r)
# [1] TRUE
# res(carbon_sum) == res(glad_r)
# [1] FALSE FALSE

carbon_masked <- resample(carbon_norm, glad_r, method = "bilinear") %>% 
  mask(., glad_r, maskvalues = 0, updatevalue = NA)

## export
writeRaster(carbon_masked, 
            here("data/h4h_suitability_carbon_modal_1km_sumfirst.tif"),
            overwrite = TRUE)


