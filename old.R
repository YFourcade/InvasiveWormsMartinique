
###############
##### OLD #####
###############

ggplot() +
  geom_spatraster(data = env_all[[1]], na.rm = T, maxcell = 1E5, alpha = .7) +
  # geom_sf(data = martinique, fill = NA) +
  geom_point(data = occ %>%
               mutate(Origin = ifelse(Origin == "Native", "Native species", "Exotic species")) %>%
               mutate(Origin = factor(Origin, levels = c("Native species", "Exotic species"))),
             aes(x = X, y = Y, color = Habitat), size = .5) +
  ggh4x::facet_nested_wrap(Origin ~ Species, strip = ggh4x::strip_split(position = c("left", "top")),
                           nest_line = element_line(color = 'black')) +
  scale_fill_distiller("Elevation", na.value = "white", palette = "Greys", guide = "none") +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  theme(strip.text.x = element_text(face = "italic", margin=ggplot2::margin(b=5), size = unit(7, 'pt')),
        strip.text.y = element_text(face = "bold", margin=ggplot2::margin(r= 5), angle = 90),
        legend.position = "bottom")

ggsave(file = ".././Fig1.pdf", width = 5.8, height = 7)

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





############
# TIDY SDM 
############

### 1. select occurrence data ====
occ.tmp <- st_as_sf(occ %>% filter(Species == i), coords = c("X", "Y"))
st_crs(occ.tmp) <- 4326

### 2. sample background data ====
occ.bg.tmp <- occ.tmp %>% sample_pseudoabs(
  n = nb.back.pts,
  raster = env,
  method = c("dist_min", km2m(.5))
)
# ggplot() +
#   geom_spatraster(data = env[[1]], maxcell = 5000) +
#   geom_sf(data = occ.bg.tmp, aes(col = class))

### 3. extract environmental values at occurrences and background data ====
occ.bg.env.tmp <- bind_cols(occ.bg.tmp, extract(env, occ.bg.tmp, ID = FALSE)) %>% na.omit

### 4. prepare SDM ====
sdm_rec <- recipe(occ.bg.env.tmp, formula = class ~ .)

sdm_workflow <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = sdm_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf(tune = "all"),
      # boosted tree model (gbm) specs with tuning
      gbm = sdm_spec_boost_tree(tune = "all"),
      # maxent specs with tuning
      maxent = sdm_spec_maxent(tune = "all")
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

### 4. Prepare stratification into training / test datasets ====
cv.table <- vfold_cv(occ.bg.env.tmp, v = 4)
# autoplot(cv.table)
# check_splits_balance(cv.table, "class")

### 5. tuning ====
sdm_workflow <-
  sdm_workflow %>%
  workflow_map("tune_grid",
               resamples = cv.table,
               grid = 20,
               metrics = sdm_metric_set(), verbose = TRUE
  )

# autoplot(sdm_workflow)
# 
# explainer_SDM <- explain_tidysdm(sdm_ensemble)

### 5. ensemble model ====
sdm_ensemble <- simple_ensemble() %>%
  add_member(sdm_workflow, metric = "boyce_cont")

# sdm_ensemble %>% collect_metrics()


### 6. project ensemble SDM ====
pred.ens <- predict_raster(sdm_ensemble, aggregate(env,10), fun = "weighted_mean")
ggplot() +
  geom_spatraster(data = pred.ens, aes(fill = weighted_mean)) +
  scale_fill_terrain_c() +
  # plot presences used in the model
  geom_sf(data = occ.bg.env.tmp %>% filter(class == "presence"), size = 1)

### 11. Find and apply 10% presence threshold ====
n10 <- ceiling(length(extract(pred.ens, occ.tmp)[,2]) * 0.1)
pred.ens.bin <- pred.ens > sort(extract(pred.ens, occ.tmp)[,2])[n10]










##########


#### MaxEnt ====
print(paste("Modelling:", i, "; MaxEnt"))

mx.cv <- ENMevaluate(
  occs = occ.tmp,
  envs = env,
  bg = back %>% rename(X = x, Y = y),
  algorithm = "maxnet",
  partitions = "user",
  tune.args = list(
    fc = c("L","LQ","LQH","H","LQHP"), 
    rm = seq(.5,5,.5)
  ),
  user.grp = cv.r,
  categoricals = categ,
  parallel = T,
  numCores = nc
)


res.mx <- eval.results(mx.cv) %>% 
  rename(cbiTest = cbi.val.avg, 
         aucTest = auc.val.avg)
sel.mx <- res.mx %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)

mods.maxnet <- eval.results(mx.cv)
mx <- enm.maxnet@predict(mods.maxnet$fc.H_rm.0.5, env)



mx.cv <- trainByCrossValid(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  folds = cv.r,
  trainFx = trainMaxNet,
  classes = "default",
  regMult = seq(.5,5,.5),
  forceLinear = T,
  testClasses = T,
  cores = nc,
  verbose = T
)

res.mx <- bind_rows(mx.cv$tuning) %>% 
  group_by(regMult, classes) %>% 
  summarise(cbiTest = mean(cbiTest, na.rm = T), 
            aucTest = mean(aucTest, na.rm = T),
            tssTest = mean(tssTest, na.rm = T),
            msssThold = mean(msssThold, na.rm = T))
sel.mx <- res.mx %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)

mx <- trainMaxNet(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  classes = sel.mx$classes,
  regMult = sel.mx$regMult
)

mxMap <- predictEnmSdm(mx, env)
mxPred <- data_pred %>% mutate(predicted = predictEnmSdm(mx, as.data.frame(data_pred)))

#### RF ====
print(paste("Modelling:", i, "; Random forest"))

www <- table(dat.tmp$presBg) / nrow(dat.tmp)
www <- c(rep(www[1]/nrow(occEnv), nrow(occEnv)),
         rep(www[2]/nrow(backEnv), nrow(backEnv)))

# calculating the class weights and sample size
prNum <- as.numeric(table(dat.tmp$presBg)["1"]) # number of presences
bgNum <- as.numeric(table(dat.tmp$presBg)["0"]) # number of backgrounds
cwt <- c("1" = 1, "0" = prNum / bgNum)
cwt <- c(rep(cwt["1"], prNum),  rep(cwt["0"], bgNum))
samsize <- c("0" = prNum, "1" = prNum)


rf.cv <- trainByCrossValid(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  folds = cv.r,
  trainFx = trainRF,
  numTrees = c(250,500,1000,1250,1500),
  mtryIncrement = 2,
  cores = nc,
  sampsize = samsize,
  binary = T,
  verbose = T
)

res.rf <- bind_rows(rf.cv$tuning) %>% 
  group_by(numTrees) %>% 
  summarise(cbiTest = mean(cbiTest, na.rm = T), 
            aucTest = mean(aucTest, na.rm = T),
            tssTest = mean(tssTest, na.rm = T),
            msssThold = mean(msssThold, na.rm = T))
sel.rf <- res.rf %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)

mrf <- trainRF(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  numTrees = sel.rf$numTrees,
  mtryIncrement = 2,
  w = cwt,
  cores = nc
)

rfMap <- predictEnmSdm(mrf, aggregate(env, 20))
rfPred <- data_pred %>% mutate(predicted = predictEnmSdm(mrf, as.data.frame(data_pred)))

#### GLM ====
# glm.cv <- trainByCrossValid(
#   data = dat.tmp,
#   resp = "presBg",
#   preds = names(dat.tmp)[2:ncol(dat.tmp)],
#   folds = cv.r,
#   trainFx = trainGLM,
#   scale = T,
#   quadratic = TRUE,
#   interaction = TRUE,
#   interceptOnly = TRUE,
#   presPerTermInitial = 2,
#   presPerTermFinal = 2,
#   maxTerms = 5,
#   cores = nc,
#   verbose = T
# )
# 
# res.glm <- bind_rows(glm.cv$tuning) %>%
#   group_by(k) %>%
#   summarise(cbiTest = unique(cbiTest[AICc == min(AICc)]),
#             aucTest = unique(aucTest[AICc == min(AICc)]),
#             tssTest = unique(tssTest[AICc == min(AICc)]))
# sel.glm <- res.glm %>% ungroup %>%
#   summarise(cbiTest = mean(cbiTest, na.rm = T),
#             aucTest = mean(aucTest, na.rm = T),
#             tssTest = mean(tssTest, na.rm = T))
# 
# mglm <- trainGLM(
#   data = dat.tmp,
#   resp = "presBg",
#   preds = names(dat.tmp)[2:ncol(dat.tmp)],
#   quadratic = T,
#   interaction = T,
#   interceptOnly = F,
#   scale = T,
#   presPerTermInitial = 2,
#   presPerTermFinal = 2,
#   maxTerms = 5,
#   cores = nc
# )
# 
# sel.glm$formula = as.character(mglm$formula)[3]
# 
# glmMap <- predictEnmSdm(mglm, env)
# glmPred <- data_pred %>% mutate(predicted = predictEnmSdm(mglm, as.data.frame(data_pred)))

#### BRT  ====
print(paste("Modelling:", i, "; Boosted regression trees"))

brt.cv <- trainByCrossValid(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  folds = cv.r,
  trainFx = trainBRT,
  learningRate = c(1e-04, 0.001, 0.01),
  treeComplexity = c(1,3,5,7,9),
  bagFraction = c(5:8)/10,
  tries = 5,
  tryBy = c("learningRate", "treeComplexity", "maxTrees", "stepSize"),
  cores = nc,
  verbose = T
)

res.brt <- bind_rows(brt.cv$tuning) %>% 
  group_by(treeComplexity, bagFraction, learningRate) %>% 
  summarise(cbiTest = mean(cbiTest, na.rm = T), 
            aucTest = mean(aucTest, na.rm = T),
            tssTest = mean(tssTest, na.rm = T),
            msssThold = mean(msssThold, na.rm = T))
sel.brt <- res.brt %>% ungroup %>% filter(cbiTest == max(cbiTest, na.rm = T)) %>% sample_n(1)

mbrt <- trainBRT(
  data = dat.tmp,
  resp = "presBg",
  preds = names(dat.tmp)[2:ncol(dat.tmp)],
  learningRate = sel.brt$learningRate,
  treeComplexity = sel.brt$treeComplexity,
  bagFraction = sel.brt$bagFraction,
  anyway = T,
  cores = nc
)

brtMap <- predictEnmSdm(mbrt, env)
brtPred <- data_pred %>% mutate(predicted = predictEnmSdm(mbrt, as.data.frame(data_pred)))

### 6.  Doing  Ensemble  Modelling ====
# ensemble predictions
# preds <- c(mxMap,brtMap,rfMap,glmMap)
preds <- rast(lapply(
  list(mxMap,rfMap,brtMap),
  climateStability::rescale0to1
))
names(preds) <- c("MaxEnt", "RandomForest", "BRT")

weights <- c(sel.mx$cbiTest, sel.rf$cbiTest, sel.brt$cbiTest)
pred.ens <- weighted.mean(preds, weights)
pred.sd <- stdev(preds)


### 8.  Merge evaluation end prediction results ====
# evaluation & settings
results_tmp <- bind_cols(
  Species = i,
  bind_rows(
    cbind.data.frame(
      algorithm = "Maxnet",
      sel.mx %>% dplyr::select(c("cbiTest", "aucTest", "tssTest", "msssThold")),
      sel.mx %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% unite(setting_values),
      settings = sel.mx %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% names %>% paste(collapse = ' ')
    ),
    cbind.data.frame(
      algorithm = "RF",
      sel.rf %>% dplyr::select(c("cbiTest", "aucTest", "tssTest", "msssThold")),
      sel.rf %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% unite(setting_values),
      settings = sel.rf %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% names %>% paste(collapse = ' ')
    ),
    cbind.data.frame(
      algorithm = "BRT",
      sel.brt %>% dplyr::select(c("cbiTest", "aucTest", "tssTest", "msssThold")),
      sel.brt %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% unite(setting_values),
      settings = sel.brt %>% dplyr::select(!c("cbiTest", "aucTest", "tssTest", "msssThold")) %>% names %>% paste(collapse = ' ')
    )
  )
)

# predictions
predictions_linear <- bind_rows(
  MaxEnt = mxPred,
  `Random Forest` = rfPred,
  BRT = brtPred,
  .id = "algorithm"
)


### 7 Find and apply thresholdx ====
# 10% presence
# max TSS
preds.bin_trainSe90 <- c()
preds.bin_maxTSS <- c()
for(mm in 1:3){
  n10.tmp <- ceiling(length(extract(preds[[mm]], occ.tmp)[,2]) * 0.1)
  preds.bin.tmp_trainSe90 <- preds[[mm]] >= sort(extract(preds[[mm]], occ.tmp)[,2])[n10.tmp]
  preds.bin.tmp_maxTSS <- preds[[mm]] >= results_tmp[mm,"msssThold"]
  
  preds.bin_maxTSS <- c(preds.bin_maxTSS, preds.bin.tmp_maxTSS)
  preds.bin_trainSe90 <- c(preds.bin_trainSe90, preds.bin.tmp_trainSe90)
}
preds.bin_maxTSS <- rast(preds.bin_maxTSS)
preds.bin_trainSe90 <- rast(preds.bin_trainSe90)

preds.bin_trainSe90.uncert <- sum(preds.bin_trainSe90)
preds.bin_maxTSS.uncert <- sum(preds.bin_maxTSS)


n10 <- ceiling(length(extract(pred.ens, occ.tmp)[,2]) * 0.1)
pred.ens.bin <- pred.ens > sort(extract(pred.ens, occ.tmp)[,2])[n10]

### 9.  Write everything on disk ====
# rasters
predictions_raster <- c(
  pred.ens,
  pred.sd,
  preds.bin_trainSe90.uncert,
  preds.bin_maxTSS.uncert
)
names(predictions_raster) <- c("EMwmean", "EM_sd", "EM_binary", "EM_ca")
writeRaster(predictions_raster, paste0("../Results/ENM_predictions_",i,".tif"), overwrite = T)

write_csv(results_tmp, paste0("../Results/ENM_results_",i,".csv"))
write_csv(predictions_linear, paste0("../Results/ENM_responseCurves_",i,".csv"))
