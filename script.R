# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
# ================================= #

# Load libraries ====
library(tidyverse)
library(readxl)
library(terra)
library(rnaturalearth)
library(ecospat)
library(FactoMineR)
library(ENMeval)
library(ENMTools)
library(enmSdmX)
library(sf)
library(retry)

# Load data ====
## occurrence records ====
occ <- read_csv("../species_records.csv") %>% filter(Y<15) %>% 
  distinct()

### Select species ====
occ %>% group_by(Species) %>% 
  summarise(n = n())

occ <- occ %>% left_join(
  occ %>% group_by(Species) %>% 
    summarise(n = n())
) %>% filter(n > 10) %>% dplyr::select(-n)

## environmental variables ====
elevation <- rast("../../Data_GIS/MNT/MNT_Martinique.grd")

nebulosity <- rast("../../Data_GIS/Nebulosite/Cloud_cover_Martinique_UTM.asc")
crs(nebulosity) <- "epsg:5490"
nebulosity <- project(nebulosity, elevation)
nebulosity[nebulosity == 0] <- NA

elevation <- mask(elevation, nebulosity)

soil <- vect("../../Data_GIS/Soil/Soil_Martinique.shp")
soil$LéGENDE__E <- gsub("SOL BRUN", "SOLS BRUN", soil$LéGENDE__E)
soil <- project(soil, elevation)
soil <- rasterize(soil, elevation, field = "LéGENDE__E")

OS <- vect("../../Data_GIS/OCCUPATION_SOL/OCCUPATION_SOL.shp")
OS[grepl("CS2.1", OS$CODE_CS),"CODE_CS"] <- 1
OS[OS$CODE_CS != 1, "CODE_CS"] <- 0
OS <- project(OS, elevation)
Forest <- rasterize(OS, elevation, field = "CODE_CS")

### Merge and plot ====
env <- c(elevation, nebulosity, soil, Forest)
plot(env)

# Map data ====
martinique <- ne_countries(scale = 10, type = "map_units", country = 'france', returnclass = "sf") %>%
  filter(name_fr == "Martinique") %>% .$geometry

ggplot() + 
  geom_sf(data = martinique) +
  geom_point(data = occ, aes(x = X, y = Y, color = Species)) +
  facet_wrap(~ Species) +
  theme_void() +
  theme(legend.position = 'none')

# ENM / SDM ====
## set parameters ====
# env <- aggregate(env, 20, fun = "modal") #for testing only
nb.back.pts <- 1000 # for final models, use more background points

## launch loop over all species ====
results_all <- c()
pred_all <- c()
pred.bin_all <- c()
for(i in setdiff(unique(occ$Species), results_all$Species)){
  print(paste("Modelling:", i))
  
  # select occurrence data
  occ.tmp <- occ %>% filter(Species == i) %>% dplyr::select(X, Y)
  occ.tmp <- st_as_sf(occ.tmp, coords = c("X", "Y"), crs = getCRS('WGS84'))
  
  # extract environmental values at occurrences and background data
  occEnv <- terra::extract(env, occ.tmp)[,-1] %>% na.omit()
  occ.tmp <- occ.tmp[as.numeric(rownames(occEnv)),]
  
  back <- spatSample(vect(martinique), size = nb.back.pts)
  backEnv <- terra::extract(env, back)[,-1] %>% na.omit()
  back <- back[as.numeric(rownames(backEnv))]
  
  presBg <- data.frame(
    presBg = c(
      rep(1, nrow(occEnv)),
      rep(0, nrow(backEnv))
    )
  )
  
  dat.tmp <- rbind(occEnv, backEnv)
  dat.tmp <- cbind(presBg, dat.tmp)
  
  # dat.tmp <- dat.tmp %>% 
  #   fastDummies::dummy_cols(
  #     select_columns = "LéGENDE__E",
  #     remove_selected_columns = T,
  #     omit_colname_prefix = T,
  #     remove_most_frequent_dummy  = T
  #   )
  # names(dat.tmp) <- gsub(" ", "_", names(dat.tmp))
  # 
  # names(dat.tmp) <- gsub('[^[:alnum:] ]','',names(dat.tmp))
  # names(dat.tmp) <- substr(names(dat.tmp), nchar(names(dat.tmp)) - 6 + 1, nchar(names(dat.tmp)))
  # names(dat.tmp) <- substr(names(dat.tmp), )
  
  # define cross-validation folds
  k = 4
  # folds <- blockCV::cv_spatial(
  #   st_as_sf(
  #     bind_rows(
  #       st_coordinates(occ.tmp) %>% as_tibble, 
  #       geom(back)[,c("x", "y")] %>% as_tibble %>% rename(X = x, Y = y)
  #     ) %>% mutate(presBg), 
  #     coords = c("X", "Y"), 
  #     crs = getCRS('WGS84')
  #   ), 
  #   column = "presBg",
  #   k = k,
  #   extend = .5,
  #   hexagon = F
  # )
  
  folds <- get.block(st_coordinates(occ.tmp), geom(back)[,c("x","y")])
  folds.all <- c(folds$occs.grp, folds$bg.grp)
  
  # check data & folds 
  plot(martinique, border = 'gray', main = paste(k, 'geo-folds'))
  plot(back, pch = 3, col = folds$bg.grp + 1, add = TRUE)
  plot(st_geometry(occ.tmp), pch = 20 + folds$occs.grp, bg = folds$occs.grp + 1, add = TRUE)
  
  # MaxNet
  print("MaxEnt")
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      mx.cv <- trainByCrossValid(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        folds = folds.all,
        trainFx = trainMaxNet,
        classes = "default",
        regMult = c(.5,5,.5),
        forceLinear = T,
        testClasses = T,
        cores = 8,
        verbose = 1
      )
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  res.mx <- bind_rows(mx.cv$tuning) %>% 
    group_by(regMult, classes) %>% 
    summarise(cbiTest = mean(cbiTest, na.rm = T), 
              aucTest = mean(cbiTest, na.rm = T),
              tssTest = mean(tssTest, na.rm = T))
  sel.mx <- res.mx %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)
  
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      mx <- trainMaxNet(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        classes = sel.mx$classes,
        regMult = sel.mx$regMult
      )      
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  
  mxMap <- predictEnmSdm(mx, env)
  
  # Random Forest
  print("Random Forest")
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      m.rf.cv <- trainByCrossValid(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        folds = folds.all,
        trainFx = trainRF,
        numTrees = c(250, 500, 750, 1000),
        cores = 8,
        verbose = 1
      )      
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  res.rf <- bind_rows(m.rf.cv$tuning) %>% 
    group_by(numTrees) %>% 
    summarise(cbiTest = mean(cbiTest, na.rm = T), 
              aucTest = mean(cbiTest, na.rm = T),
              tssTest = mean(tssTest, na.rm = T))
  sel.rf <- res.rf %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)
  
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      m.rf <- trainRF(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        numTrees = res.rf$numTrees,
        cores = 8
      )      
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  rfMap <- predictEnmSdm(m.rf, env)
  
  # # Boosted regression trees
  # print("BRT")
  # boolFalse<-F
  # while(boolFalse==F){
  #   tryCatch({
  #     mbrt.cv <- trainByCrossValid(
  #       data = dat.tmp,
  #       resp = "presBg",
  #       preds = names(dat.tmp)[2:ncol(dat.tmp)],
  #       folds = folds.all,
  #       trainFx = trainBRT,
  #       learningRate = c(.0001, .001, .01, .1),
  #       treeComplexity = c(1,3,5,7,9,11),
  #       bagFraction = c(5:7)/10,
  #       cores = 8,
  #       verbose = 1
  #     )      
  #     boolFalse = T
  #   },error=function(e){
  #   },finally={})
  # }
  # 
  # res.mbrt <- bind_rows(mbrt.cv$tuning) %>% 
  #   group_by(learningRate, treeComplexity, bagFraction) %>% 
  #   summarise(cbiTest = mean(cbiTest, na.rm = T), 
  #             aucTest = mean(cbiTest, na.rm = T),
  #             tssTest = mean(tssTest, na.rm = T))
  # sel.mbrt <- res.mbrt %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)
  # 
  # boolFalse<-F
  # while(boolFalse==F){
  #   tryCatch({
  #     mbrt <- trainBRT(
  #       data = dat.tmp,
  #       resp = "presBg",
  #       preds = names(dat.tmp)[2:ncol(dat.tmp)],
  #       learningRate = sel.mbrt$learningRate,
  #       treeComplexity = sel.mbrt$treeComplexity,
  #       bagFraction = sel.mbrt$bagFraction
  #     )      
  #     boolFalse = T
  #   },error=function(e){
  #   },finally={})
  # }
  # 
  # brtMap <- predictEnmSdm(mbrt, env)
  
  # GLM
  print("GLM")
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      mglm.cv <- trainByCrossValid(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        folds = folds.all,
        trainFx = trainGLM,
        quadratic = F,
        interaction = T,
        interceptOnly = F,
        cores = 8,
        verbose = 1,
        presPerTermInitial = 1,
        presPerTermFinal = 1,
        maxTerms = 5
      )
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  res.mglm <- bind_rows(mglm.cv$tuning) %>% 
    group_by(k) %>% 
    summarise(cbiTest = unique(cbiTest[AICc == min(AICc)]), 
              aucTest = unique(aucTest[AICc == min(AICc)]),
              tssTest = unique(tssTest[AICc == min(AICc)]))
  sel.mglm <- res.mglm %>% ungroup %>% 
    summarise(cbiTest = mean(cbiTest, na.rm = T), 
              aucTest = mean(aucTest, na.rm = T),
              tssTest = mean(tssTest, na.rm = T))
  
  boolFalse<-F
  while(boolFalse==F){
    tryCatch({
      mglm <- trainGLM(
        data = dat.tmp,
        resp = "presBg",
        preds = names(dat.tmp)[2:ncol(dat.tmp)],
        quadratic = F,
        interaction = T,
        interceptOnly = F,
        presPerTermInitial = 1,
        presPerTermFinal = 1,
        maxTerms = 5
      )      
      boolFalse = T
    },error=function(e){
    },finally={})
  }
  
  glmMap <- predictEnmSdm(mglm, env)
  
  # ensemble predictions
  # preds <- c(mxMap,brtMap,rfMap,glmMap)
  preds <- rast(lapply(
    list(mxMap,rfMap,glmMap),
    climateStability::rescale0to1
  ))
  weights <- c(sel.mx$cbiTest, sel.rf$cbiTest, sel.mglm$cbiTest)
  pred.ens <- weighted.mean(preds, weights)
  pred.sd <- stdev(preds)
  
  # Find and apply 10% presence threshold
  n10 <- ceiling(length(terra::extract(pred.ens, occ.tmp)[,2]) * 0.1)
  pred.ens.bin <- pred.ens > sort(terra::extract(pred.ens, occ.tmp)[,2])[n10]
  
  # merge evaluations / settings
  results_tmp <- bind_cols(
    Species = i,
    bind_rows(
      cbind.data.frame(
        algorithm = "Maxnet",
        sel.mx %>% dplyr::select(c("cbiTest", "aucTest", "tssTest")),
        sel.mx %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% unite(setting_values),
        settings = sel.mx %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% names %>% paste(collapse = ' ')
      ),
      cbind.data.frame(
        algorithm = "RF",
        sel.rf %>% dplyr::select(c("cbiTest", "aucTest", "tssTest")),
        sel.rf %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% unite(setting_values),
        settings = sel.rf %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% names %>% paste(collapse = ' ')
      ),
      # cbind.data.frame(
      #   algorithm = "BRT",
      #   sel.mbrt %>% dplyr::select(c("cbiTest", "aucTest", "tssTest")),
      #   sel.mbrt %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% unite(setting_values),
      #   settings = sel.mbrt %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest")) %>% names %>% paste(collapse = ' ')
      # ),
      cbind.data.frame(
        algorithm = "GLM",
        sel.mglm %>% dplyr::select(c("cbiTest", "aucTest", "tssTest")),
        setting_values = NA,
        settings = NA
      )
    )
  )
  
  # merge
  results_all <- rbind.data.frame(results_all, results_tmp)
  pred_all <- c(pred_all, pred.ens)
  pred.bin_all <- c(pred.bin_all, pred.ens.bin)
  
}

# Restructure results ====
results_all <- results_all %>% left_join(
  occ %>% dplyr::select(1,4) %>% 
    group_by(Species) %>% 
    summarise(Introduced = unique(Introduced))
)

pred_all <- rast(pred_all)
names(pred_all) <- unique(occ$Species)

pred.bin_all <- rast(pred.bin_all)
names(pred.bin_all) <- unique(occ$Species)

## write rasters on disk ====
writeRaster(pred_all, "../SDM_predictions.tif", overwrite = T)
writeRaster(pred.bin_all, "../SDM_predictions.binary.tif", overwrite = T)

## reload rasters ====
pred_all <- rast("../SDM_predictions.tif")
pred.bin_all <- rast("../SDM_predictions.binary.tif")

# Check model evaluation ====
# ggplot(results_all, aes(y = cbiTest, x = n.occ)) +
#   geom_point() + theme_bw()
# ggplot(results_all, aes(y = auc.val.avg, x = n.occ)) +
#   geom_point() + theme_bw()
# ggplot(results_all, aes(y = cbi.val.avg, x = auc.val.avg)) +
#   geom_point() + theme_bw()

results_all %>% 
  group_by(Introduced) %>% 
  summarise(cbi = mean(cbiTest, na.rm = T),
            cbi.se = plotrix::std.error(cbiTest, na.rm = T)) %>% 
  ggplot(aes(y = cbi, x = Introduced, fill = Introduced,
             ymin = cbi - cbi.se, ymax = cbi + cbi.se)) +
  geom_bar(stat = 'identity') + 
  geom_errorbar(width = .3) +
  scale_y_continuous("Continuous Boyce index", limits = c(0, .8), expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = "../eval.pdf", height = 4, width = 5)

# Plot model predictions ====
## continuous ====
as.data.frame(pred_all, xy = T) %>% as_tibble %>% 
  pivot_longer(cols = -c(1:2), names_to = "Species", values_to = "Suitability") %>% 
  ggplot(aes(x = x, y = y, fill = Suitability)) +
  geom_raster() +
  facet_wrap(~ Species) +
  scale_fill_continuous("Suitability", type = "viridis") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf()

ggsave(filename = "../Suitability_maps.pdf", height = 7, width = 6.5)
ggsave(filename = "Suitability_maps.png", height = 7, width = 6.5)

## binary ====
as.data.frame(pred.bin_all, xy = T) %>% as_tibble %>% 
  pivot_longer(cols = -c(1:2), names_to = "Species", values_to = "Suitability") %>% 
  mutate(Suitability = as.factor(Suitability)) %>% 
  ggplot(aes(x = x, y = y, fill = Suitability)) +
  geom_raster() +
  facet_wrap(~ Species) +
  scale_fill_manual(
    "Suitability", 
    values = c("lightgrey", "orange3"), 
    labels = c("Unsuitable", "Suitable")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf()

# Overlap between introduced and native species ====
## all species ====
as.data.frame(pred.bin_all, xy = T) %>% as_tibble %>% 
  pivot_longer(cols = -c(1:2), names_to = "Species", values_to = "Suitability") %>% 
  group_by(x, y) %>% 
  summarise(Overlap = sum(Suitability)) %>% 
  mutate(Overlap = factor(Overlap, levels = 0:max(Overlap))) %>%
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  geom_raster() +
  scale_fill_manual(
    "No. species",
    values = c("lightgrey", "orange", "orange2", "orange3", "red4", "red3", "red")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf()

## introduced and native ====
as.data.frame(pred.bin_all, xy = T) %>% as_tibble %>% 
  pivot_longer(cols = -c(1:2), names_to = "Species", values_to = "Suitability") %>% 
  left_join(
    occ %>% dplyr::select(1,4) %>% 
      group_by(Species) %>% 
      summarise(Introduced = unique(Introduced))
  ) %>% group_by(x, y, Introduced) %>% 
  summarise(Overlap = sum(Suitability)) %>% 
  mutate(Overlap = factor(Overlap, levels = 0:max(Overlap)),
         Introduced = ifelse(Introduced == "Yes", "Introduced", "Native")) %>%
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  facet_grid(~ Introduced) +
  geom_raster() +
  scale_fill_manual(
    "No. species",
    values = c("lightgrey", "orange", "orange2", "orange3", "red4", "red3", "red")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf()


## introduced vs. native ====
overlap.plot <- c()
for(i in unique(as.data.frame(occ)[occ$Introduced == "Yes","Species"])){
  for(j in unique(as.data.frame(occ)[occ$Introduced == "No","Species"])){
    rast.temp <- subset(pred.bin_all, c(i,j))
    overlap.plot.tmp <- as.data.frame(rast.temp, xy = T) %>% as_tibble %>% 
      rename(Introduced = 3, Native = 4) %>% 
      mutate(Overlap = ifelse(Native == 1 & Introduced == 1, "Overlap",
                              ifelse(Native == 1 & Introduced == 0, "Native only",
                                     ifelse(Native == 0 & Introduced == 1, "Introduced only", 
                                            "No species")))) %>% 
      mutate(Introduced = i, Native = j)
    overlap.plot <- bind_rows(overlap.plot, overlap.plot.tmp)
  }
}

### plot a matrix of overlap area ====
overlap.plot %>% group_by(Introduced, Native) %>% 
  summarise(n.overlap = length(Overlap[Overlap == "Overlap"])) %>% 
  ggplot(aes(x = Introduced, y = Native, 
             fill = n.overlap, 
             label = n.overlap)) +
  geom_tile() +
  geom_text(size = 3, col = "white") +
  scale_x_discrete("") +
  scale_y_discrete("") +
  scale_fill_gradient(low = "lightgrey", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.pos = "none") +
  coord_fixed()

### plot overal areas ====
overlap.plot %>% mutate(Overlap = factor(Overlap, levels = c("No species", "Native only", "Introduced only", "Overlap"))) %>% 
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  facet_grid(Native ~ Introduced) +
  geom_raster() +
  scale_fill_manual(
    "Overlap",
    values = c("lightgrey", "green4", "yellow3", "red3")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf()

ggsave(filename = "../overlap_maps.pdf", width = 11, height = 6)
ggsave(filename = "overlap_maps.png", width = 11, height = 6)


### synthetic invasion risk ====
overlap.plot %>% group_by(x, y) %>% 
  summarise(n.overlap = length(Overlap[Overlap == "Overlap"])) %>% 
  ggplot(aes(x = x, y = y, fill = n.overlap)) +
  geom_raster() +
  scale_fill_viridis_c(
    "No. pairs of\noverlapping species",
    option = "magma"
  ) +
  theme_void() +
  coord_sf()

ggsave(filename = "../synthetic_maps.pdf", width = 6, height = 5)
ggsave(filename = "synthetic_maps.png", width = 6, height = 5)

# Metrics of geographical and environmental overlap ====
## extract environmental values ====
### at occurrences ====
occ_env <- bind_cols(
  Species = occ[,c(1,4)], 
  terra::extract(env, occ[,3:2])
)

### at background points ====
# env <- env[[c(1,2,4)]] # do not keep
back.env <- spatSample(env, na.rm = T, size = nb.back.pts)

## overlap in environmental space (D_env) ====
famd.cal <- FAMD(
  bind_rows(back.env %>% na.omit, 
            occ_env %>% dplyr::select(colnames(back.env)) %>% na.omit),
  ind.sup = (nrow(back.env %>% na.omit)+1):
    (nrow(back.env %>% na.omit)+nrow(occ_env %>% dplyr::select(colnames(back.env)) %>% na.omit)),
  ncp = 2,
  graph = F
)

scores.back <- famd.cal$ind$coord[1:nrow(back.env %>% na.omit),]
scores.sp <- famd.cal$ind.sup$coord %>% as_tibble %>% 
  mutate(Species = occ_env %>% na.omit %>% pull(Species))

ggplot() + 
  geom_point(data = as_tibble(scores.back), aes(x = Dim.1, y = Dim.2), color = "lightgrey") +
  geom_point(data = scores.sp, aes(x = Dim.1, y = Dim.2, color = Species)) +
  theme_bw()

overlap_env <- c()
for(i in occ %>% filter(Introduced == "Yes") %>% pull(Species) %>% unique){
  for(j in occ %>% filter(Introduced == "No") %>% pull(Species) %>% unique){
    if(i != j){
      scores.sp1 <- scores.sp %>% filter(Species == i)
      scores.sp2 <- scores.sp %>% filter(Species == j)
      
      sp1_ecospat <- ecospat.grid.clim.dyn(
        glob = scores.back,
        glob1 = scores.back,
        sp = scores.sp1[,1:2],
        R = 100
      )
      
      sp2_ecospat <- ecospat.grid.clim.dyn(
        glob = scores.back,
        glob1 = scores.back,
        sp = scores.sp2[,1:2],
        R = 100
      )
      
      D <- ecospat.niche.overlap(sp1_ecospat, sp2_ecospat, cor = T)$D
      
      overlap_env <- rbind.data.frame(
        overlap_env,
        cbind.data.frame(
          Introduced = i , 
          Native = j, 
          D_env = D
        ) 
      )
    }
  }
}


p_d_env <- ggplot(data = overlap_env, aes(y=Introduced, x=Native, fill=D_env, label = round(D_env,2))) +
  geom_tile() +
  geom_text(size = 2.5, col = "white") +
  scale_y_discrete("Introduced species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[env])), low = "green3", high = "dark red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_fixed()


## overlap in geographical space (D_geo) ====

overlap_geo <- c()
for(i in occ %>% filter(Introduced == "Yes") %>% pull(Species) %>% unique){
  for(j in occ %>% filter(Introduced == "No") %>% pull(Species) %>% unique){
    if(i != j){
      sdm1 <- pred_all[[names(pred_all) == i]]
      sdm2 <- pred_all[[names(pred_all) == j]]
      
      D <- raster.overlap(sdm1, sdm2)$D
      
      overlap_geo <- rbind.data.frame(
        overlap_geo,
        cbind.data.frame(
          Introduced = i , 
          Native = j, 
          D_geo = D
        ) 
      )
      
    }
  }
}

p_d_geo <- ggplot(data = overlap_geo, aes(y=Introduced, x=Native, fill=D_geo, label = round(D_geo,2))) +
  geom_tile() +
  geom_text(size = 2.5, col = "white") +
  scale_y_discrete("Introduced species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[geo])), low = "blue", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_fixed()

## plot ====
cowplot::plot_grid(p_d_env, p_d_geo, nrow = 1, labels = "AUTO")

ggsave(filename = "../overlap_indices.pdf", width = 8, height = 4.5)

# merge
overlap <- left_join(overlap_env, overlap_geo)
overlap_unique <- overlap %>% 
  mutate(pair = purrr::map2_chr(Introduced, Native, ~paste(sort(c(.x, .y)), collapse = " - "))) %>%
  group_by(pair) %>%
  summarise(Introduced = first(Introduced),
            Native = first(Native),
            D_env = unique(D_env),
            D_geo = unique(D_geo)) %>% 
  dplyr::select(-pair)

plot(overlap_unique %>% dplyr::select(D_env, D_geo), pch = 16)
cor.test(overlap_unique %>% dplyr::select(D_env) %>% pull(D_env),
         overlap_unique %>% dplyr::select(D_geo) %>% pull(D_geo),
         method = "spearman")
