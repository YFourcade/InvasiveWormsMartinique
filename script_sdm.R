# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# MODELS
# ================================= #

# Load libraries ====
library(tidyverse)
library(terra)
library(rnaturalearth)
library(sf)
library(biomod2)
library(readxl)
# library(blockCV)

# Load data ====
## occurrence records ====
# read files
occ <- list.files(".././SP_NICHE_OVERLAP", full.names = T)[-1] %>% 
  lapply(.,read_excel) %>% bind_rows() %>% select(nom_espece, coord_X, coord_Y) %>% 
  rename(Species = 1, X = coord_X, Y = coord_Y) %>% 
  mutate(Species = stringr::word(Species, 1, 2)) 

# %>% 
#   mutate(X = ifelse(X < -50, -X, X))
  
# read species characteristics
sp_char <- read_excel(".././SP_NICHE_OVERLAP/SELECTION_SP.xlsx")
names(sp_char) <- sp_char[1,]
sp_char_a <- sp_char[2:9,] %>% mutate(Origin = "Exotic") %>% 
  rename(Species = 3) %>% select(-DATA)
sp_char_b <- sp_char[12:24,] %>% mutate(Origin = "Exotic") %>% 
  rename(Species = 3) %>% select(-DATA)
sp_char_tot <- bind_rows(sp_char_a, sp_char_b) %>% rename(Arboreal = 1, Soil = 2) %>% 
  mutate(Arboreal = ifelse(is.na(Arboreal), 0 ,1),
         Soil = ifelse(is.na(Soil), 0 , 1))

# merge
occ %>% left_join(sp_char_tot)

# occ <- read_csv() %>% filter(Y<15, Y>14, X < -60.8) %>% 
#   distinct()

### Select species ====
occ %>% group_by(Species) %>% 
  summarise(n = n())

# occ <- occ %>% left_join(
#   occ %>% group_by(Species) %>% 
#     summarise(n = n())
# ) %>% filter(n > 10) %>% dplyr::select(-n)

occ %>% group_by(Habitat, Origin) %>% 
  summarise(n.sp = length(unique(Species)))

## environmental variables ====
env <- rast(".././env.tiff")

legend_soil <- structure(
  list(
    ID = 0:8, 
    LéGENDE__E = c("FERRISOLS", 
                   "SOLS A ALLOPHANE ( ANDOSOLS )", 
                   "SOLS A ALLUVIONS", "SOLS BRUN - ROUILLE A HALLOYSITE", "SOLS FERSIALLITIQUES", 
                   "SOLS PEU EVOLUES SUR CENDRES", "VERTISOLS", "Zone urbaine", 
                   NA)),
  class = "data.frame", 
  row.names = c(NA, -9L)
)
env <- categories(env, layer = 3, value = legend_soil)
names(env) <- c("Elevation", "Nebulosity", "Soil", "Forest", "Agriculture", "Water")
env <- env[[-3]]

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
nb.back.pts <- 5000 # for final models, use more background points
back <- spatSample(vect(martinique), size = nb.back.pts)
backEnv <- extract(env, back)[,-1] %>% na.omit()
back <- back[as.numeric(rownames(backEnv))]

## launch loop over all species ====
library(doFuture)
registerDoFuture()

nc = 6 # no. of parallel stuff

results_all <- vector('list', length(unique(occ$Species)))
pred_all <- c()
pred.bin_all <- c()
n = 0
for(i in setdiff(unique(occ$Species), results_all$Species)){
  n = n + 1
  print(paste("Modelling:", i))
  
  
  ### 1. select occurrence data ====
  occ.tmp <- occ %>% filter(Species == i) %>% dplyr::select(X, Y)
  occ.tmp <- vect(occ.tmp, geom = c("X", "Y"), crs = enmSdmX::getCRS('WGS84'))
  
  ### 2. extract environmental values at occurrences and background data ====
  occEnv <- extract(env, occ.tmp)[,-1] %>% na.omit()
  occ.tmp <- occ.tmp[as.numeric(rownames(occEnv)),]
  
  presBg <- data.frame(
    presBg = c(
      rep(1, nrow(occEnv)),
      rep(0, nrow(backEnv))
    )
  )
  dat.tmp <- rbind(occEnv, backEnv)
  
  ###  3.  Formatting  Data ====
  myBiomodData  <-  BIOMOD_FormatingData(
    resp.name  =  i,
    resp.var  =  presBg,
    expl.var  =  dat.tmp,
    resp.xy = rbind.data.frame(geom(occ.tmp)[,3:4], geom(back)[,3:4]),
    dir.name = "D:/Models_Martinique/"
  )
  
  ### 4. Prepare stratification into training / test datasets ====
  cv.r <- bm_CrossValidation(bm.format = myBiomodData,
                             strategy = "block")
  
  ### 5.  Doing  Modelling ====
  # set options
  my.opts <- bm_ModelingOptions(data.type = 'binary', 
                     bm.format = myBiomodData,
                     strategy = 'bigboss')
  my.opts@options$MAXNET.binary.maxnet.maxnet@args.values$`_allData_allRun`$addsamplestobackground <- F
  my.opts@options$MAXNET.binary.maxnet.maxnet@args.values$`_allData_allRun`$regmult <- 2

  plan("multisession", workers = nc)
  
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    modeling.id = 'AllModels',
    models = c('GLM', 'RF', 'GBM', 'XGBOOST', 'MAXNET'),
    metric.eval = c('TSS', 'ROC', 'BOYCE'),
    scale.models = F,
    CV.strategy = "user.defined",
    CV.user.table = cv.r,
    CV.do.full.models = F,
    OPT.strategy = 'tuned',
    OPT.user.base = "bigboss",
    # OPT.user = my.opts,
    nb.cpu = nc
  )
  
  ### 6.  Doing  Ensemble  Modelling ====
  myBiomodEM  <-  BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen  =  'all',
    em.by  =  'all',
    metric.select  =  'BOYCE',
    metric.eval  =  c('TSS','ROC', 'BOYCE'),
    em.algo = c('EMcv', 'EMwmean'), 
    var.import = 10,
    nb.cpu = nc
  )
  
  ### 7. evaluation ====
  eval.temp <- get_evaluations(myBiomodModelOut)
  eval.temp <- eval.temp %>% group_by(algo, metric.eval) %>% 
    mutate(Species = i, .before = 1, n.occ = nrow(occEnv))
  
  ### 8. variable importance ====
  varImp <- get_variables_importance(myBiomodEM)
  varImp.summary.temp <- varImp %>% group_by(Variable = expl.var) %>% 
    mutate(Species = i, .before = 1)
  
  ### 9. response curves ====
  respCurves <- bm_PlotResponseCurves(myBiomodModelOut, do.plot = F)$tab
  respCurves.summary.temp <- respCurves %>% separate(pred.name, into = c("1","2","3","Algorithm"), sep = "_") %>% 
    select(expl.name, expl.val, Algorithm, pred.val) %>% 
    group_by(expl.name, expl.val, Algorithm) %>% 
    summarise(pred.val.mean = mean(pred.val),
              pred.val.se = plotrix::std.error(pred.val))
  
  ### 10.  Spatial Projection ====
  # single models
  myBiomodProjection  <-  BIOMOD_Projection(
    bm.mod  =  myBiomodModelOut,
    proj.name = "projection",
    new.env  =  env,
    models.chosen = 'all',
    build.clamping.mask  =  F,
    nb.cpu = nc
  )
  
  plan(sequential)
  # ensemble model
  myBiomodProjectionEns <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,
    bm.proj = myBiomodProjection
  )
  
  ### 11. Find and apply 10% presence threshold ====
  pred.ens <- unwrap(myBiomodProjectionEns@proj.out@val)
  n10 <- ceiling(length(terra::extract(pred.ens, occ.tmp)[,2]) * 0.1)
  pred.ens.bin <- pred.ens[[2]] > sort(terra::extract(pred.ens[[2]], occ.tmp)[,2])[n10]
  
  ### 12. merge lots of stuff ====
  results_tmp <- list(eval.temp, varImp.summary.temp, respCurves.summary.temp)
  
  # merge
  results_all[[n]] <- results_tmp
  pred_all <- c(pred_all, pred.ens)
  pred.bin_all <- c(pred.bin_all, pred.ens.bin)
  
}

# Restructure and write numerical results ====
Var_importance <- bind_rows(
  lapply(
    1:length(results_all),
    function(x)pluck(results_all, x, 2)
  )
)
write_csv(Var_importance, ".././Var_importance.csv")

SDM_evaluations <- bind_rows(
  lapply(
    1:length(results_all),
    function(x)pluck(results_all, x, 1)
  )
)
write_csv(SDM_evaluations, ".././SDM_evaluations.csv")

Var_response <- bind_rows(
  lapply(
    1:length(results_all),
    function(x)pluck(results_all, x, 3) %>% mutate(Species = unique(occ$Species)[[x]], .before = 1)
  )
)
write_csv(Var_response, ".././Var_response.csv")


# write rasters on disk ====
pred_all <- rast(pred_all) 
pred_all.cv <- rast(pred_all) 

pred_all.mean <- subset(pred_all, grepl("mean", names(pred_all)))
names(pred_all.mean) <- unique(occ$Species)

writeRaster(pred_all, "../SDM_predictions.tif", overwrite = T)

pred_all.cv <- subset(pred_all, grepl("cv", names(pred_all.cv)))
names(pred_all.cv) <- unique(occ$Species)

pred.bin_all <- rast(pred.bin_all)
names(pred.bin_all) <- unique(occ$Species)
writeRaster(pred.bin_all, "../SDM_predictions.binary.tif", overwrite = T)


