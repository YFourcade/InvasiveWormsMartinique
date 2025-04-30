# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# DISTRUBTION OF BROMELIADS
# ================================= #
source("functions_sdm.R")

# Load libraries ====
library(tidyverse)
library(readxl)
library(terra)
library(rnaturalearth)
library(tidyterra)
library(flexsdm)
library(sf)
library(parallel)

# Load data ====
## occurrence records ====
# read files
pts_bro <- read_sf(".././Data/EmplacementBromeliacees/Point_Bromeliacees.shp")
trans_bro <- read_sf(".././Data/EmplacementBromeliacees/Transect_Bromeliacees.shp") %>% st_centroid

pts_bro_all <- bind_rows(pts_bro, trans_bro)
pts_bro_all <- st_transform(pts_bro_all, crs = "EPSG:4326")
pts_bro_all <- pts_bro_all %>% st_coordinates() %>% unique


## environmental variables ====
env_all <- rast(".././env.tiff")

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
env_all <- categories(env_all, layer = 3, value = legend_soil)
names(env_all) <- c("Elevation", "Nebulosity", "Soil", "Forest", "Agriculture")


# ENM / SDM ====
## set parameters ====
nb.back.pts <- 10000 # for final models, use ca. 10000 background points

## Model bromeliads distribution  ====
nc = 8 # no. of parallel stuff

occ.tmp <- pts_bro_all %>% as_tibble()
occ.tmp$pr_ab = 1

env <- env_all
categ = "Soil"

### 2.  Formatting  folder ====
# none

### 3. sample background with thickening ====
back.tmp <- sample_background(
  data = occ.tmp,
  x = "X",
  y = "Y",
  n = nb.back.pts,
  method = c("thickening", width = 20000),
  rlayer = env
)

occ_back <- rbind(occ.tmp, back.tmp)

### 3. Prepare stratification into training / test datasets ====
occ_back.folds <- occ_back %>% mutate(
  part = unlist(ENMeval::get.block(
    occs = occ.tmp[,1:2],
    bg = back.tmp[,1:2]
  ))
)

ggplot() +
  geom_spatraster(data = env[[1]], na.rm = T, maxcell = 1E5) +
  geom_point(data = occ_back.folds %>% filter(pr_ab == 0), aes(x = X, y = Y, col = as.factor(part)),
             shape = 16, size = .5) +
  geom_point(data = occ_back.folds %>% filter(pr_ab == 1), aes(x = X, y = Y, col = as.factor(part)),
             shape = 22, size = 3, fill = "white") +
  scale_fill_distiller("Elevation", na.value = "white", palette = "Greys", guide = "none") +
  scale_color_brewer("Evaluation folds", palette = "Set1") +
  theme_void()


### 4. Extract environmental values ====
occ_back_env <- sdm_extract(
  data = occ_back.folds,
  x = "X",
  y = "Y",
  env_layer = env,
  filter_na = T
)

### 5. Tune models ====

# RF
tune_grid.rf <- expand.grid(
  mtry = seq(1, 7, 2), 
  ntree = c(200, 500, 1000, 1500)
)

model.rf.tune <- tune_raf(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  grid = tune_grid.rf,
  thr = "sensitivity",
  n_cores = nc,
  metric = "BOYCE"
)

# pred.rf <- sdm_predict(
#   models =  model.rf.tune,
#   pred = env,
# )$raf
pred.rf <- predict(env, model.rf.tune$model, type = "prob", index = 2)
names(pred.rf) <- "raf"

pred.rf$mtp <- sdm_threshold(pred.rf$raf, occ.tmp[,-3], type = "mtp", binary = T)
pred.rf$p10 <- sdm_threshold(pred.rf$raf, occ.tmp[,-3], type = "p10", binary = T)
pred.rf$raf <- climateStability::rescale0to1(pred.rf$raf)

# MaxEnt
tune_grid.mx <- expand.grid(
  regmult = seq(.5,5,1),
  classes = c("l", "lq", "h", "lqh", "lqhp")
)

model.mx.tune <- tune_max(
  data = occ_back_env %>% filter(pr_ab == 1),
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  background = occ_back_env %>% filter(pr_ab == 0),
  grid = tune_grid.mx,
  thr = "sensitivity",
  n_cores = nc,
  metric = "BOYCE"
)

# pred.mx <- sdm_predict(
#   models =  model.mx.tune,
#   pred = env,
# )$max
pred.mx <- enmSdmX::predictMaxNet(model.mx.tune$model, env, type = "cloglog")
names(pred.mx) <- "max"

pred.mx$mtp <- sdm_threshold(pred.mx$max, occ.tmp[,-3], type = "mtp", binary = T)
pred.mx$p10 <- sdm_threshold(pred.mx$max, occ.tmp[,-3], type = "p10", binary = T)
pred.mx$max <- climateStability::rescale0to1(pred.mx$max)

# BRT
tune_grid.brt <- expand.grid(
  n.trees = c(20, 50, 100, 200, 500),
  shrinkage = c(0.001, 0.01, 0.1, 0.5),
  n.minobsinnode = 10
)

model.brt.tune <- tune_gbm(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  grid = tune_grid.brt,
  thr = "sensitivity",
  n_cores = nc
)

# pred.brt <- sdm_predict(
#   models =  model.brt.tune,
#   pred = env,
# )$gbm
pred.brt <- predict(env, model.brt.tune$model, type = "response")
pred.brt <- mask(pred.brt, env[[1]])
names(pred.brt) <- "gbm"

pred.brt$mtp <- sdm_threshold(pred.brt$gbm, occ.tmp[,-3], type = "mtp", binary = T)
pred.brt$p10 <- sdm_threshold(pred.brt$gbm, occ.tmp[,-3], type = "p10", binary = T)
pred.brt$gbm <- climateStability::rescale0to1(pred.brt$gbm)


# NNeT
tune_grid.net <- expand.grid(
  size = c(2, 4, 6, 8, 10),
  decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10)
)

model.net.tune <- tune_net(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  grid = tune_grid.net,
  thr = "sensitivity",
  n_cores = nc,
  metric = "BOYCE"
)

# pred.net <- sdm_predict(
#   models =  model.net.tune,
#   pred = env
# )$net
pred.net <- predict(env, model.net.tune$model)
names(pred.net) <- "net"

pred.net$mtp <- sdm_threshold(pred.net$net, occ.tmp[,-3], type = "mtp", binary = T)
pred.net$p10 <- sdm_threshold(pred.net$net, occ.tmp[,-3], type = "p10", binary = T)
pred.net$net <- climateStability::rescale0to1(pred.net$net)

# GLM
model.glm.1 <- fit_glm(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  select_pred = F,
  poly = 2,
  inter_order = 1,
  thr = "sensitivity"
)

model.glm.2 <- glm(model.glm.1$model$formula, data = occ_back_env, 
                   family = 'binomial', na.action = na.fail)

clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", nc), type = clusterType))
clusterExport(clust, "occ_back_env")

table.mod.glm <- MuMIn::dredge(
  model.glm.2,
  cluster = clust,
  m.lim= c(
    1, 
    max(2,6)
  ),
  subset = dc("Elevation", "I(Elevation^2)") && dc("Forest", "I(Forest^2)") &&
    dc("Bromeliads", "I(Bromeliads^2)") && dc("Agriculture", "I(Agriculture^2)")  && 
    dc("Nebulosity", "I(Nebulosity^2)")
)

stopCluster(clust)

model.glm.3 <- MuMIn::get.models(table.mod.glm, subset = 1)[[1]]

model.glm.4 <- fit_glm(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  fit_formula = model.glm.3$formula,
  thr = "sensitivity"
)

# pred.glm <- sdm_predict(
#   models =  model.glm.4,
#   pred = env,
# )$glm
pred.glm <- predict(env, model.glm.4$model, type = 'response')
names(pred.glm) <- "glm"

pred.glm$mtp <- sdm_threshold(pred.glm$glm, occ.tmp[,-3], type = "mtp", binary = T)
pred.glm$p10 <- sdm_threshold(pred.glm$glm, occ.tmp[,-3], type = "p10", binary = T)
pred.glm$glm <- climateStability::rescale0to1(pred.glm$glm)


# GAM
model.gam <- fit_gam(
  data = occ_back_env,
  response = "pr_ab",
  predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
  predictors_f = categ,
  partition = "part",
  select_pred = F,
  k=4,
  thr = "sensitivity"
)

model.gam.2 <- mgcv::bam(model.gam$model$formula, data = occ_back_env, 
                         family = 'binomial', discrete = TRUE, na.action = na.fail)

table.mod.gam <- MuMIn::dredge(
  model.gam.2,
  m.lim= c(
    1, 
    max(2,floor(nrow(occ_back_env[occ_back_env$pr_ab==1,"pr_ab"])/3))
  )
)
model.gam.3 <- MuMIn::get.models(table.mod.gam, subset = 1)[[1]]

try(
  model.gam.4 <- fit_gam(
    data = occ_back_env,
    response = "pr_ab",
    predictors = gratia::term_names(model.gam.3),
    predictors_f = categ,
    partition = "part",
    select_pred = F,
    thr = "sensitivity"
  )
)
if (is.null(model.gam.4)){
  model.gam.4 <- model.gam
}

# pred.gam <- sdm_predict(
#   models =  model.gam,
#   pred = env,
# )$gam
pred.gam <- predict(env, model.gam.4$model, type = 'response')
names(pred.gam) <- "gam"

pred.gam$mtp <- sdm_threshold(pred.gam$gam, occ.tmp[,-3], type = "mtp", binary = T)
pred.gam$p10 <- sdm_threshold(pred.gam$gam, occ.tmp[,-3], type = "p10", binary = T)
pred.gam$gam <- climateStability::rescale0to1(pred.gam$gam)


### 6. summarise evaluations and parameters ====
eval_results <- sdm_summarize(
  list(model.rf.tune, 
       model.mx.tune,
       model.brt.tune, 
       model.net.tune, 
       model.glm.4, 
       model.gam.4)
) %>% 
  select(2:6, BOYCE_mean, BOYCE_sd, AUC_mean, AUC_sd, TSS_mean, TSS_sd)

tune_results <- sdm_summarize(
  list(model.rf.tune, 
       model.mx.tune,
       model.brt.tune, 
       model.net.tune, 
       model.glm.4, 
       model.gam.4)
) %>% 
  select(2, mtry, ntree, regmult,classes, n.trees, shrinkage, n.minobsinnode, size, decay)
tune_results <- tune_results %>% mutate_at(.vars = -1, .funs = as.character) %>% 
  pivot_longer(-1) %>% mutate(parameter = paste(name, "=", value)) %>% na.omit %>% 
  select(-name, -value) %>% 
  bind_rows(
    bind_cols(model = "gam", parameter = paste("formula =",as.character(model.gam$model$formula)[[3]])),
    bind_cols(model = "glm", parameter = paste("formula =",as.character(model.glm.4$model$formula)[[3]])),
  ) %>% 
  group_by(model) %>% 
  summarise(parameters = paste(parameter, collapse = " ; "))

eval_results <- eval_results %>% left_join(tune_results) %>% select(-c(2:5))


### 7. Ensemble modelling ====
preds <- c(pred.rf$raf, pred.mx$max, pred.brt$gbm, pred.net$net, pred.glm$glm, pred.gam$gam)
weights <- eval_results %>% pull(BOYCE_mean)
# weights.std <- weights[weights > 0] * length(weights[weights > 0])/sum(weights[weights > 0])

preds <- preds[[weights > 0]]
pred.ens.mean <- weighted.mean(preds, weights[weights > 0], na.rm = T)

pred.ens.ca.mtp  <- sum(c(pred.rf$mtp,
                          pred.mx$mtp,
                          pred.brt$mtp,
                          pred.net$mtp,
                          pred.glm$mtp,
                          pred.gam$mtp)[[weights > 0]])
pred.ens.ca.p10  <- sum(c(pred.rf$p10,
                          pred.mx$p10,
                          pred.brt$p10,
                          pred.net$p10,
                          pred.glm$p10,
                          pred.gam$p10)[[weights > 0]])

# ens.model <- fit_ensemble(
#   models = list(model.rf.tune, 
#                 model.mx.tune,
#                 model.brt.tune, 
#                 model.net.tune, 
#                 model.glm.4, 
#                 model.gam)[weights > 0],
#   ens_method = "meanw",
#   thr = "sensitivity",
#   thr_model = "sensitivity",
#   metric = "BOYCE"
# )
# 
# ens.model.perf <- ens.model$performance %>% 
#   select(1, thr_value, BOYCE_mean, BOYCE_sd, AUC_mean, AUC_sd, TSS_mean, TSS_sd)

# find and apply thresholds on ensemble model
pred.ens.pa.mtp <- sdm_threshold(pred.ens.mean, occ.tmp[,-3], type = 'mtp', binary = T)
pred.ens.pa.p10 <- sdm_threshold(pred.ens.mean, occ.tmp[,-3], type = 'p10', binary = T)

pred.final <- c(pred.ens.mean, pred.ens.ca.mtp, pred.ens.ca.p10, pred.ens.pa.mtp, pred.ens.pa.p10)
names(pred.final) <- c("Weighted mean", "Committee average (MTP)", "Committee average (P10)",
                       "Binary ensemble (MTP)", "Binary ensemble (P10)")


### 8. Response curves ====
linear_predictions <- c()
list.mods <- list(
  "Random Forest" = model.rf.tune, 
  "MaxEnt" = model.mx.tune,
  "Boosted regression trees" = model.brt.tune, 
  "Neural network" = model.net.tune, 
  "GLM" = model.glm.4, 
  "GAM"= model.gam
)[weights > 0]
for(ii in 1:length(list.mods)){
  p1 <- p_pdp(model = list.mods[[ii]]$model, training_data = occ_back_env)
  p1.1 <- lapply(1:length(p1), function(x){p1[[x]] %>% ggplot_build() %>% .$data})
  p1.2 <- unlist(lapply(1:length(p1), function(x){p1[[x]]$label$x}))
  
  var.list <- c()
  for(jj in p1.2){
    if(jj != "Soil"){
      var.list <- c(var.list, rep(jj, 100))
    }else{
      var.list <- c(var.list, rep(jj, length(unique(occ_back_env$Soil) %>% na.omit)))
    }
  }
  
  p1.dat <- bind_rows(p1.1) %>% mutate(Variable = var.list) %>% 
    select(x,y,Variable)
  linear_predictions <- rbind.data.frame(
    linear_predictions,
    p1.dat %>% 
      mutate(Algorithm = names(list.mods)[[ii]], Species = "Bromeliads")
  )
}


### 9. Write everything ====
# rasters
writeRaster(pred.final, "../Variables/Bromeliads_predictions.tif", overwrite = T)
write_csv(eval_results %>% mutate(Species = "Bromeliads", .before = 1), "../Variables/Bromeliads_results.csv")
write_csv(linear_predictions %>% mutate(Species = "Bromeliads", .before = 1), "../Variables/Bromeliads_responseCurves.csv")

