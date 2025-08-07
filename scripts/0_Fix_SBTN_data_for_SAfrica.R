##### Purpose: Script to adjust the SBTN natural land layer in South Africa
## SBTN uses the SA NLC as a data input in South Africa but misinterpreted some of the built-up classes
## e.g. residential area with trees (as per SA NLC) were classed as non-natural tree cover rather than built-up in SBTN
## To address this, a reclassification of the SBTN layer in S Africa using the input 2020 SA NLC was done

## Author: Luke Wilson
## Email: lwilson@conservation.org

## Last updated 7 July 2025

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
  furrr,              # loads both future and purrr packages
  RCGLS,              # download Copernicus LC and fractional cover data
  RCurl,
  dplyr)              

#### Step 1: Load and merge SBTN data for SSA ####
## Directory for the SBTN datasets (change as needed)
sbtn_dir <- ("C:/Users/lwilson/OneDrive - Conservation International Foundation/GIS_data/Land_Cover/SBTN_Natural_Lands_Data")

## Quick fxn to merge and crop/mask SBTN tiles
sbtn_prep <- function(files) {
  ## first read in all rasters as a stack
  stack <- lapply(files, rast)
  
  ## then merge
  merge <- do.call(terra::merge, stack)
  
  ## now crop & mask to ROI
  ssa_v <- ssa_sf %>% 
    vect() %>% 
    project(., crs(merge))
  
  cropped <- crop(merge, ssa_v, mask = TRUE)
  
  return(cropped)
}

## Get list of SBTN tiles
sbtn_files <- list.files(
  path = file.path(sbtn_dir, "Tiles"),
  pattern = "*.tif$",
  full.names = TRUE
)

## Run fxn
sbtn_merge <- sbtn_prep(sbtn_files)
plot(sbtn_merge) # check output

## Export
writeRaster(sbtn_merge, 
            file.path(sbtn_dir, "SBTN_SSA.tif"),
            overwrite = TRUE)


#### Step 2: Load, merge and resample SA NLC data to match SBTN ####
## Directory for SA NLC tiles (change as needed)
sanlc_dir <- ("C:/Users/lwilson/OneDrive - Conservation International Foundation/GIS_data/Land_Cover/SA_NLC")

## Merge SA NLC tiles
sanlc_prep <- function(files) {
  ## first read in all rasters as a stack
  stack <- lapply(files, rast)
  
  ## then merge
  merge <- do.call(terra::merge, stack)
  
  return(merge)
}

## Get list of SA NLC tiles
sanlc_files <- list.files(
  path = file.path(sanlc_dir, "Tiles"),
  pattern = "*.tif$",
  full.names = TRUE
)


## Run fxn
sanlc_merge <- sanlc_prep(sanlc_files)
plot(sanlc_merge) # check output

## Export
writeRaster(sanlc_merge, 
            file.path(sanlc_dir, "SA_NLC_2020_merged.tif"),
            overwrite = TRUE)


## Resample SA NLC to match SBTN
sbtn_crop <- terra::crop(sbtn_merge, sanlc_merge) # crop SBTN layer to SA NLC to reduce processing extent
plot(sbtn_crop) # check crop
sanlc_resampled <- terra::resample(sanlc_merge, sbtn_crop, method = "mode")



#### Step 3: Update the SBTN layer using the SA NLC in South Africa ####

## Mask SBTN to the resampled SA NLC extent (maybe not necessary but want to ensure it is masked exactly 
## to the resampled SA NLC to avoid errors later)
sbtn_SA_crop <- terra::mask(sbtn_crop, sanlc_resampled)
plot(sbtn_SA_crop)


## Correct built-up areas in SBTN based on SA NLC
built_mask <- sanlc_resampled >= 47 & sanlc_resampled <= 72 # these are various built LC classes in SA NLC 
sbtn_SA_edit <- sbtn_SA_crop # creating a copy to keep the original layer 
sbtn_SA_edit[built_mask] <- 13  # changes areas within the masked area to the value for built-up in SBTN


## Mosaic edited SBTN values for SA back on to the SSA raster
# This assigns NA values to the edited raster outside of the edited region 
sbtn_edit_aligned <- resample(sbtn_SA_edit, sbtn_merge, method="near")
# This gives NA areas the value from the original SBTN data
updated_sbtn <- cover(sbtn_edit_aligned, sbtn_merge)
plot(updated_sbtn) # check new layer

## Export edited SBTN layer
writeRaster(updated_sbtn, 
            file.path(sbtn_dir, "SBTN_SSA_edited.tif"),
            overwrite = TRUE)
