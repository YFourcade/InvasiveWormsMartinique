
###############
##### OLD #####
###############


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
      regMult = seq(.5,5,.5),
      forceLinear = F,
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
            aucTest = mean(aucTest, na.rm = T),
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