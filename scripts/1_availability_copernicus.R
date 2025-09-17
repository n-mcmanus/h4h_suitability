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
cop_dir <- file.path(data, "Land_Cover/Copernicus")

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


### 1B. Reclassify Fractional Tree Coverage -------------------
## Now evaluate Copernicus fractional tree and cover to get more precise thresholds
tree_r <- rast(file.path(cop_dir, 'fractional_coverage/trees', 'copernicus_tree_fraction_SSA_2019.tif'))

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


### 1C. Isolate non-herding bare areas ---------------------
### Finally, remove areas that are considered "bare" in Copernicus & 
### not within Ramona rangeland extent

## Read in Ramona data
range_ext_r <- rast(file.path(data, 'livestock_grazing/Ramona_Data/rangeland_max_extent_LT_NDVI_smoothed_3_AFRICA.tiff'))

## crop and mask to ROI
ssa_reproj <- project(ssa_v, range_ext_r)
range_ext_crop <- crop(range_ext_r, ssa_reproj, mask = TRUE)

## Reproject/disaggregate/resample ramona data to match
range_ext_reproj <- project(range_ext_crop, crs(discrete_r), method = "near")

factor <- floor(res(range_ext_reproj)[1] / res(discrete_r)[1])

range_ext_resamp <- disagg(range_ext_reproj, factor, "near") %>% 
  resample(., discrete_r, "near")

## Isolate areas totally bare AND outside grazing extent
bare_mask <- as.numeric((range_ext_resamp == 0) & (discrete_r == 60))

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
            file.path(slope_dir, "availability_ssa_slope_100m.tif"),
            datatype = "INT1U",
            overwrite = TRUE)



## 3. Remove additional built-up & impervious areas --------------------------

## Read in most recent LULC avail if needed
avail_r <- rast(file.path(slope_dir, "availability_ssa_slope_100m.tif"))

## Read in ROI if needed
ssa_v <- vect(file.path(data, "ROI/ssa_ne.shp")) %>% 
  project(., avail_r)


## We'll remove mining areas from the Maus et al., 2022 polygons
## Read in and crop to only ROI
mines_v <- vect(file.path(data, "Land_Cover/Maus_2022_mining", "global_mining_polygons_v2.gpkg")) %>% 
  project(., ssa_v) %>% 
  crop(., ssa_v)

### We also want to remove built-up areas from SBTN data layer. These data capture more
### urban and built-up environments than Copernicus LULC alone.
sbtn_r <- rast(file.path(data, "Land_Cover/SBTN/SBTN_SSA_edited.tif")) %>% 
  ## match projection and then resolution of avail_r 
  project(., y = crs(avail_r)) %>% 
  aggregate(., fact = 3, fun = "modal") %>% #SBTN at ~30m; avail lyr is 100m
  resample(., avail_r, method = "mode")


## Now mask availability layer with mines & SBTN
avail_no_built <-
  mask(avail_r, mines_v, inverse = TRUE, touches = TRUE, updatevalue = 0) %>%
  mask(., sbtn_r, maskvalues = 13, updatevalue = 0)


## Export availability layer
writeRaster(avail_no_built,
            file.path(out_dir, "availability_ssa_100m.tif"),
            datatype = "INT1U",
            overwrite = TRUE)




## 5. Carbon benefits in available areas ----------------------------------

### 5A. Naturebase ------------------------------
### According to priority order of Naturebase's logic rules, need to mask and
### add relevant layers in proper order

layers <- c(
  "grs_agc",  #avoided grassland conversion
  "grs_asc",  #avoided shrubland conversions
  "wet_ipm",   #improved peatland mgmt
  "grs_scg",  #increased soil carbon in grazing lands
  "grs_sfm",  #savannah fire management
  "grs_grr"  #grassland restoration
)

## directory of prepped naturebase layers
nb_dir <- file.path(data, "carbon/naturebase/africa")

## loop through each layer
for (i in seq_along(layers)) {
  layer <- layers[i]
  
  ## read in layer and make 0 values NA for easier masking
  r <- rast(file.path(nb_dir, paste0(layer, "_tco2eha_ssa.tif")))
  r[r==0] <- NA
  
  if (i == 1) {
    carbon_total <- r
  } else {
    r_mask <- mask(r, carbon_total, inverse = TRUE)
    carbon_total <- sum(carbon_total, r_mask, na.rm = TRUE)
  }
}


## Read in availability layer (if needed)
avail_r <- rast(file.path(out_dir, "availability_ssa_100m.tif"))
avail_r[avail_r == 0] <- NA # remove values of 0

## mask carbon data to available areas
carbon_masked <- mask(carbon_total, avail_r, maskvalue = 0, updatevalue = NA)

## save output
writeRaster(carbon_masked, 
            file.path(out_dir, "availability_carbon_naturebase_100m.tif", ), 
            overwrite = TRUE)



### 5B. CSA's (Heidi's) Carbon Data ----------------------------
### Also mask availability w/carbon data generated from Heidi Hawkins

## Read in availability layer (if needed)
avail_r <- rast(file.path(out_dir, "availability_ssa_100m.tif"))

ssa_v <- project(ssa_v, avail_r)

## Directory for data
afc_dir <- file.path(data, "carbon/Heidi_AfRange")

## Only want the 3 rasters of mean difference of carbon
       ## **NOTE**: Double-check this is right!
afc_fnames <- list.files(afc_dir,
                         pattern = "MEAN_DIFF.*\\.tif$",
                         full.names = TRUE)

## Add up ABG, BGB, and SOC rasters
afc_r <- rast(afc_fnames)
afc_sum_r <- app(afc_r, sum, na.rm = TRUE) %>% 
  project(., y = crs(avail_r))


## Approach 1: Keep larger resolution for carbon mask ---

## aggregate availability layer to match resolution of carbon data
factor <- floor(res(afc_sum_r)[1] / res(avail_r)[1])
avail_agg_r <- aggregate(avail_r, factor, fun = "modal")

## resample carbon data to match availability, then mask
avail_agg_r[avail_agg_r == 0] <- NA
afc_mask <- resample(afc_sum_r, avail_agg_r, method = "near") %>% 
  mask(., avail_agg_r)

plot(afc_mask)
plot(ssa_v, add = TRUE)

## save output
writeRaster(afc_mask, 
            file.path(out_dir, "carbon_meandiff_avail_lowres.tif"),
            overwrite = TRUE)


## Approach 2: Get smaller resolution outputs for carbon mask ---
factor <- ceiling(res(afc_sum_r)[1] / res(avail_r)[1])
afc_disagg <- disagg(afc_sum_r, factor, method = "near") ## using near to retain "original resolution" as best as possibe

## resample carbon data to match availability, then mask 
afc_resamp <- resample(afc_disagg, avail_r, method = "bilinear")
afc_mask <- mask(afc_resamp, ssa_v) %>% 
  mask(., avail_r, maskvalues = 0)

plot(afc_mask)
plot(ssa_v, add = TRUE)

## save output
writeRaster(afc_mask,
            file.path(out_dir, "carbon_meandiff_avail_100m.tif"),
            overwrite = TRUE)




