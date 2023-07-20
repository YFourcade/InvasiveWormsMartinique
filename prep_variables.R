# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# PREPARATION OF VARIABLES
#
# ================================= #

# Load libraries ====
library(tidyverse)
library(terra)
library(landscapemetrics)
setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

# Load data ====

## environmental variables ====
# elevation
elevation <- rast("../../Data_GIS/MNT/MNT_Martinique.grd")

# could cover
nebulosity <- rast("../../Data_GIS/Nebulosite/Cloud_cover_Martinique_UTM.asc")
crs(nebulosity) <- "epsg:5490"
nebulosity <- project(nebulosity, elevation)
nebulosity[nebulosity == 0] <- NA

# mask elevation
elevation <- mask(elevation, nebulosity)

# land cover
OS <- vect("../../Data_GIS/OCCUPATION_SOL/OCCUPATION_SOL.shp")
OS <- project(OS, elevation)

# proportion agriculture
OS.agri <- OS
OS.agri[grepl("US1.1", OS$CODE_US),"CODE_US"] <- 1
OS.agri[OS.agri$CODE_US != 1, "CODE_US"] <- 0
Agri <- rasterize(OS.agri, elevation, field = "CODE_US")
Agri.prop <- focal(Agri, 41, na.rm = T, fillvalue = 0, expand = T) # moving windows of resolution = 20 grid cells (ca. 500 m)
Agri.prop <- mask(Agri.prop, elevation)

# proportion forest
OS.forest <- OS
OS.forest[grepl("CS2.1", OS.forest$CODE_CS),"CODE_CS"] <- 1
OS.forest[OS.forest$CODE_CS != 1, "CODE_CS"] <- 0
Forest <- rasterize(OS.forest, elevation, field = "CODE_CS")
Forest.prop <- focal(Forest, 41, na.rm = T, fillvalue = 0, expand = T) 
Forest.prop <- mask(Forest.prop, elevation)

# water course
water <- vect("../../Data_GIS/Hydrographie/COURS_D_EAU.shp")
water <- project(water, elevation)
dist.water <- distance(aggregate(elevation, 10), water)
dist.water <- mask(resample(dist.water, elevation), elevation)

# soil types
soil <- vect("../../Data_GIS/Soil/Soil_Martinique.shp")
soil$LéGENDE__E <- gsub("SOL BRUN", "SOLS BRUN", soil$LéGENDE__E)
soil <- project(soil, elevation)
soil <- rasterize(soil, elevation, field = "LéGENDE__E")

# merge and write
env <- c(elevation, nebulosity, soil, Forest.prop, Agri.prop, dist.water)
names(env) <- c("Elevation", "Nebulosity", "Soil", "Forest", "Agriculture", "Water")
plot(env)

terra::writeRaster(env, ".././env.tiff", overwrite = T)

