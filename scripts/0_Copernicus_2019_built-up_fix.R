## Script 0_Copernicus_2019_built-up_fix.R

### Purpose: fixing false built up LC in 2019 Copernicus ###
#  Erroneous swathes of built-up LC occur over Tanzania and Mozambique in the 
#  discrete Copernicus classification for 2019.
#  These exist in 2018 and 2017 layers too, but not 2016. So this script updates
#  the 2019 LC values in these swathes with the 2016 LC values.   

## Set up ------------------------------------------------------------------
## Install package 'pacman' if needed
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  terra,              # GIS functions
  sf,                 # vector functions
  rnaturalearth,      # get country outlines
  rnaturalearthdata,  # " "
  httr,               # Helps download GLAD data via GEE URL
  furrr)              # loads both future and purrr packages


## Set directories (change as needed)
setwd("C:/GIS_data/H4H")
dir_out <- ("C:/GIS_data/H4H")

## Load the discrete land cover rasters for 2019 and 2016
lc_2019 <- rast("C:/Users/lwilson/OneDrive - Conservation International Foundation/GIS_data/Land_Cover/Coperncius/Discrete_Classification/Copernicus_LC_SSA_2019.tif")
lc_2016 <- rast("C:/Users/lwilson/OneDrive - Conservation International Foundation/GIS_data/Land_Cover/Coperncius/Discrete_Classification/Copernicus_Discrete_2016_fix.tif")

## Load the shapefile representing the bounding boxes for the regions with erroneous data to be fixed
## These polygons were created manually in ArcGIS around the obvious swathes of false built-up pixels
fix_region <- vect("C:/Users/lwilson/Conservation International Foundation/SPARCLER-CC - Documents/SPARCLE 2.0/H4H/Data/Land_Cover/Copernicus/Discrete_Classification/Region_to_fix/built_up_fix_region_2019.shp")

## Crop and mask the 2019 and 2016 LC raster to the fix region extent
lc_2019_crop <- crop(lc_2019, fix_region)
lc_2019_mask <- mask(lc_2019_crop, fix_region)
plot(lc_2019_mask)

lc_2016_crop <- crop(lc_2016, fix_region)
lc_2016_mask <- mask(lc_2016_crop, fix_region)


## Reclassify built-up areas in the fix region (value 50) in 2019
## to assume the value from the 2016 LC.
lc_2019_mask_reclass <- ifel(lc_2019_mask == 50, lc_2016_mask, lc_2019_mask)
plot(lc_2019_mask_reclass)

## Align the reclassified layer with the full 2019 LC extent
## This assings NA values to areas outside of the fixed region
reclass_aligned <- resample(lc_2019_mask_reclass, lc_2019, method="near")


## mosaic the fixed values from the update region back on to the full 2019 layer  
## NA areas from the reclassified layer are assigned the value from the full 2019 layer
## while keeping the reclassified values in the fixed zone
updated_raster <- cover(reclass_aligned, lc_2019)
plot(updated_raster)  # to check reclassification has worked

terra::writeRaster(updated_raster, filename=file.path(dir_out, "Copernicus_fix/copernicus_discrete_LC_SSA_2019_fixed_r.tif"), overwrite = TRUE)
