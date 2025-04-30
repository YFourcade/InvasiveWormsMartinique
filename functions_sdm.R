# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# FUNCTIONS TO BE LOADED
# ================================= #
range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}


sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- terra::extract(sdm, occs)[,2] %>% na.omit
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- 0
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}
