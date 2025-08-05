## Script: 1_base_availability.R
##
## Purpose: This script generates a binary raster of where H4H activity may 
## or may not occur based on biophysical characteristics.
##
## Last updated: 01 August 2025
##
## Authors: Nick McManus (nmcmanus@conservation.org)
##          Luke Wilson (lwilson@conservation.org)
##
## Notes: Still reviewing methodology for this layer. 

## Set up -------------------------------------------------------------------
## Install package 'pacman' if needed
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  terra,              # GIS functions
  sf,                 # vector functions
  furrr)              # loads both future and purrr packages


## Set working directories (in SharePoint)
## **NOTE:** change this first one to match local SharePoint shortcut
h4h_shortcut <- "C:/Users/nmcmanus/OneDrive - Conservation International Foundation/Documents/Projects/SPARCLE"
data <- file.path(h4h_shortcut, "SPARCLER-CC - H4H/Data") #this should automatically work if properly synced
out_dir <- file.path(data, "outputs")

## Read in ROI
ssa_v <- vect(file.path(data, "ROI/ssa_ne.shp"))


## 1. Isolate LULC classes suitable for H4H ------------------------------------
## We will use several different layers of Copernicus data for this

## Directory for the Copernicus datasets
cop_dir <- file.path(data, "Land_Cover/Copernicus_nick")

### 1A. Reclassify Discrete LULC data -----------------------
## Read in
discrete_r <- rast(file.path(cop_dir, "discrete/copernicus_discrete_LC_SSA_2019_fixed.tif"))

## binary reclassification matrix
rcl_m <- matrix(c(
  0, 0,     # remove "unknown" 
  20, 1,    # keep woody shrubs
  30, 1,    # keep herbaceous veg
  40, 1,    # keep cultivated & managed veg/ag
  50, 0,    # remove urban & built up
  60, 1,    # keep bare/sparse veg (<10% veg); more finely masked later
  70, 0,    # remove snow/ice
  80, 0,    # remove perm water
  90, 1,    # keep herbaceous wetland
  100, 0,   # remove moss/lichen
  111, 0,   # remove closed forest (canopy > 70%), evergreen needle leaf 
  112, 0,   # remove closed forest, evergreen broad leaf 
  113, 1,   # keep closed forest, deciduous need leaf
  114, 1,   # keep closed forest, deciduous broad leaf
  115, 1,   # keep closed forest, mixed
  116, 1,   # keep closed forest, not matching other definitions
  121, 1,   # keep open forest (canopy < 70%), evergreen needle leaf
  122, 1,   # keep open forest, evergreen broad leaf 
  123, 1,   # keep open forest, deciduous needle leaf
  124, 1,   # keep open forest, deciduous broad leaf
  125, 1,   # keep open forest, mixed
  126, 1,   # keep open forest, not matching other def
  200, 0    # remove oceans
), ncol = 2, byrow = TRUE)

## Reclassify discrete lulc to binary availability
discrete_rcl_r <- classify(discrete_r, rcl_m, include.lowest = TRUE)

## Save intermediate raster
writeRaster(discrete_rcl_r, 
            file.path(cop_dir, "discrete/copernicus_discrete_LC_SSA_2019_binary.tif"),
            datatype = "INT1U", #binary layer
            overwrite = TRUE)


### 1B. Reclassify Fractional Coverages -------------------
## Now evaluate Copernicus fractional tree and bare ground cover to get more precise thresholds
tree_r <- rast(file.path(cop_dir, 'fractional_coverage/trees', 'copernicus_tree_fraction_SSA_2019.tif'))
bare_r <- rast(file.path(cop_dir, 'fractional_coverage/bare', 'copernicus_bare_fraction_SSA_2019.tif'))

## Reclassify to only consider areas with <80% tree canopy
rcl_tree_m <- matrix(c(
  0, 79, 1,
  80, 100, 0
), ncol = 3, byrow = TRUE)

tree_rcl_r <- classify(tree_r, rcl_tree_m, right = NA)

## Save intermediate raster
writeRaster(tree_rcl_r,
            file.path(cop_dir, 'fractional_coverage/trees', 'copernicus_tree_binary_SSA_2019.tif'),
            datatype = "INT1U",
            overwrite = TRUE)


## Reclassify to isolate areas with 100% bare ground
bare_rcl_m <- matrix(c(
  0, 99, 1,
  100, Inf, 0
), ncol = 3, byrow = TRUE)

bare_rcl_r <- classify(bare_r, bare_rcl_m, right = NA)

## save intermediate
writeRaster(bare_rcl_r,
            file.path(cop_dir, "fractional_coverage/bare", "copernicus_bare_binary_SSA_2019.tif"),
            datatype = "INT1U",
            overwrite = TRUE)


## Want to only mask out areas that are both 100% bare in Copernicus
## and not part of rangeland extent in Ramona
range_ext_r <- rast(file.path(data, 'livestock_grazing/Ramona_Data/rangeland_max_extent_LT_NDVI_smoothed_3_AFRICA.tiff'))

## crop and mask to ROI
ssa_reproj <- project(ssa_v, range_ext_r)
range_ext_crop <- crop(range_ext_r, ssa_reproj, mask = TRUE)

## Reproject/disaggregate/resample ramona data to match
range_ext_reproj <- project(range_ext_crop, crs(bare_rcl_r), method = "near")

factor <- floor(res(range_ext_reproj)[1] / res(bare_rcl_r)[1])

range_ext_resamp <- disagg(range_ext_reproj, factor, "near") %>% 
  resample(., bare_rcl_r, "near")

## Isolate areas totally bare AND outside grazing extent
bare_mask <- as.numeric((range_ext_resamp == 0) & (bare_rcl_r == 0))

lulc_avail <- 
  ## Mask out bare areas
  mask(discrete_rcl_r, bare_mask, maskvalues = 1, updatevalue = 0) %>% 
  ## Mask out high tree canopy
  mask(., tree_rcl_r, maskvalues = 0, updatevalue = 0)

## Save intermediate
writeRaster(lulc_avail,
            file.path(cop_dir, "copernicus_SSA_availability.tif"),
            datatype = "INT1U",
            overwrite = TRUE)

## Clear up memory a bit
gc();

## 2. Remove steep terrains ------------------------------------------------
## Here, we mask out any areas with a slope over 30% 
## Data comes from USGS HMDA global dataset (Africa subset available)
## Can download here: https://www.sciencebase.gov/catalog/item/591f6d17e4b0ac16dbdde1cd

## Pathway to data
slope_dir <- file.path(data, "USGS_HMDA_slope")

## List files
slope_fnames = list.files(file.path(slope_dir, "raw"),
                          pattern = "*.tif$", 
                          full.names = TRUE)

## Merge all rasters
slope_stack <- lapply(slope_fnames, rast)
slope_merged <- do.call(terra::merge, slope_stack)

## Save merged at native resolution
writeRaster(slope_merged, 
            file.path(slope_dir, "af_slope_merged_90m.tif"),
            overwrite = TRUE)


## Match CRS and resolution of LULC data
# lulc_avail <- rast(file.path(cop_dir, "copernicus_SSA_availability.tif")) #read in if needed

slope_rcl <- 
  ## First reduce area by cropping to only SSA
  crop(slope_merged, ssa_v, mask = TRUE) %>%
  ## Then match CRS and extent of Copernicus layer
  ## (very similar resolutions (~95m vs 100m), so no need to aggregate or resample separately)
  project(., lulc_avail, method = "bilinear") %>% 
  ## Finally, reclassify to remove areas with slope > 30%
  classify(., matrix(c(30, Inf, NA), ncol = 3, byrow = TRUE))


## Save binary slope output
writeRaster(slope_rcl, 
            file.path(slope_dir, "slope_ssa_100m_binary.tif"),
            datatype = "INT1U",
            overwrite = TRUE)

ssa_v <- project(ssa_v, crs(lulc_avail))

## Finally, use binary slope as mask for availability
avail_mask <- mask(lulc_avail, slope_rcl, updatevalue = 0) %>% 
  mask(ssa_v)

## Save output
writeRaster(avail_mask, 
            file.path(out_dir, "availability_ssa_slope_100m.tif"),
            datatype = "INT1U",
            overwrite = TRUE)



## 3. Remove additional built-up/impervious areas --------------------------




## 4. Carbon benefits in available areas ----------------------------------
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


## Find carbon benefits for availability layer
carbon_benefits <- function(res, carbon_r) {
  ## read in raster
  avail_r <- rast(here(sprintf("data/outputs/availability_ssa_%s.tif", res)))
  
  ## remove values of 0
  avail_r[avail_r == 0] <- NA
  
  ## mask carbon data to available areas
  carbon_masked <- mask(carbon_r, avail_r)
  
  writeRaster(carbon_masked,
              here(sprintf("data/outputs/availability_ssa_carbon_%s.tif", res)),
              overwrite = TRUE)
}

## Run fxn
carbon_benefits(res = "1km", carbon_r = r_total)








## Read in binary availability layer (generated in previous section)
avail_r <- rast(here("data/outputs/availability_ssa_"))
glad_r <- rast(here("data/outputs/glad_no_imp_slope30_2020_ssa_binary_1km.tif"))
glad_r[glad_r == 0] <- NA

carbon_masked <- mask(r_total, glad_r)

glad_r[glad_r == 1] <- 0
test <- cover(carbon_masked, glad_r)

writeRaster(test,
            here("data/h4h_availability_slope_imp_carbon_1km.tif"),
            overwrite = TRUE)














