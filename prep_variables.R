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
library(fasterRaster)
faster(grassDir = "C:/Program Files/GRASS GIS 8.4")

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

# distance to agriculture
OS.agri <- OS
OS.agri[grepl("US1.1", OS$CODE_US),"CODE_US"] <- 1
OS.agri[OS.agri$CODE_US != 1, "CODE_US"] <- 0

Agri <- rasterize(OS.agri, elevation, field = "CODE_US")
NAflag(Agri) <- 0

Agri_fst <- fast(as.numeric(Agri))
dist.agri <- distance(Agri_fst)
dist.agri <- rast(dist.agri)
dist.agri <- mask(dist.agri, elevation)

NonAgri <- rasterize(OS.agri, elevation, field = "CODE_US")
NAflag(NonAgri) <- 1

NonAgri_fst <- fast(as.numeric(NonAgri))
dist.NonAgri <- distance(NonAgri_fst)
dist.NonAgri <- rast(dist.NonAgri)
dist.NonAgri <- mask(dist.NonAgri, elevation)

dist.agri_intExt <- dist.agri - dist.NonAgri
plot(dist.agri_intExt)
hist(dist.agri_intExt)

# Agri.prop <- focal(Agri, 41, na.rm = T, fillvalue = 0, expand = T) # moving windows of resolution = 20 grid cells (ca. 500 m)
# Agri.prop <- mask(Agri.prop, elevation)

# distance to forest
OS.forest <- OS
OS.forest[grepl("CS2.1", OS.forest$CODE_CS),"CODE_CS"] <- 1
OS.forest[OS.forest$CODE_CS != 1, "CODE_CS"] <- 0

Forest <- rasterize(OS.forest, elevation, field = "CODE_CS")
NAflag(Forest) <- 0

Forest_fst <- fast(as.numeric(Forest))
dist.forest <- distance(Forest_fst)
dist.forest <- rast(dist.forest)
dist.forest <- mask(dist.forest, elevation)

NonForest <- rasterize(OS.forest, elevation, field = "CODE_CS")
NAflag(NonForest) <- 1

NonForest_fst <- fast(as.numeric(NonForest))
dist.NonForest <- distance(NonForest_fst)
dist.NonForest <- rast(dist.NonForest)
dist.NonForest <- mask(dist.NonForest, elevation)

dist.forest_intExt <- dist.forest - dist.NonForest
plot(dist.forest_intExt)
hist(dist.forest_intExt)

# Forest.prop <- focal(Forest, 41, na.rm = T, fillvalue = 0, expand = T)
# Forest.prop <- mask(Forest.prop, elevation)

# water course
# water <- vect("../../Data_GIS/Hydrographie/COURS_D_EAU.shp")
# water <- project(water, elevation)
# dist.water <- distance(elevation, water)
# dist.water <- mask(dist.water, elevation)

# soil types
soil <- vect("../../Data_GIS/Soil/Soil_Martinique.shp")
soil$LéGENDE__E <- gsub("SOL BRUN", "SOLS BRUN", soil$LéGENDE__E)
soil <- project(soil, elevation)
soil <- rasterize(soil, elevation, field = "LéGENDE__E")

# merge and write
env <- c(elevation, nebulosity, soil, dist.forest_intExt, dist.agri_intExt)
names(env) <- c("Elevation", "Nebulosity", "Soil", "Forest", "Agriculture")

plot(env)
terra::writeRaster(env, ".././env.tiff", overwrite = T)

