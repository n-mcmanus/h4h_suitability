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
glad_r <- rast(here("data/GLAD/glad_2020_ssa_modal_binary_1km.tif"))

## Get list of all NatureBase carbon layers
## For now, focusing just on SSA 
carbon_fnames <- list.files(here("data/naturebase/africa/"),
                            pattern = "ssa.tif",
                            full.names = TRUE)



### 1) IF NORMALIZING FIRST, THEN SUMMING:  -----------
## create list of each pathway raster
carbon_list <- lapply(carbon_fnames, rast)

## normalize all rasters in list
normalize_rast <- function(r) {
  min <- global(r, "min", na.rm = TRUE)[1, ]
  max <- global(r, "max", na.rm = TRUE)[1, ]
  (r - min) / (max - min)
}

carbon_norm <- map(carbon_list, normalize_rast, .progress = TRUE)

## convert from list to stack and add together
carbon_stack <- rast(carbon_norm)

carbon_sum <- app(carbon_stack, sum, na.rm = TRUE)

## mask to LULC
carbon_masked <- mask(carbon_sum, glad_r, maskvalues = 0, updatevalue = NA)

## export
writeRaster(carbon_masked, 
            here("data/h4h_suitability_carbon_modal_1km_normfirst.tif"),
            overwrite = TRUE)



### 2) IF SUMMING FIRST, THEN NORMALIZING:  ------------
## stack all the carbon rasters and add up
carbon_list <- rast(carbon_fnames)
carbon_sum <- app(carbon_list, sum, na.rm = TRUE)

## normalize the carbon raster
min <- minmax(carbon_sum)[1, ]
max <- minmax(carbon_sum)[2, ]

carbon_norm <- (carbon_sum - min) / (max - min)

## mask to LULC
carbon_masked <- mask(carbon_norm, glad_r, maskvalues = 0, updatevalue = NA)

## export
writeRaster(carbon_masked, 
            here("data/h4h_suitability_carbon_modal_1km_sumfirst.tif"),
            overwrite = TRUE)



