# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# MODELS
# ================================= #
source("functions_sdm.R")

# Load libraries ====
library(tidyverse)
library(readxl)
library(terra)
library(tidyterra)
library(flexsdm)
library(parallel)

# Load data ====
## occurrence records ====
# read files
occ <- list.files(".././Data", pattern = "OK", full.names = T) %>%
  lapply(.,read_excel) %>% bind_rows() %>% select(nom_espece, coord_X, coord_Y) %>%
  rename(Species = 1, X = coord_Y, Y = coord_X) %>%
  mutate(Species = stringr::word(Species, 1, 2),
         Species = ifelse(Species == "Dichogaster sp11", "Dichogaster sp11A", Species),
         Species = ifelse(Species == "Glossodrilus sp. nov. 3", "Glossodrilus spp", Species),
         Species = ifelse(Species == "Glossodrilus sp.", "Glossodrilus sp1", Species)) %>%
  distinct

occ <- occ %>% filter(X > -61.5)

# read species characteristics
sp_char <- read_excel(".././SELECTION_SP.xlsx", sheet = 2)

# merge
occ <- occ %>% left_join(sp_char)


### Select species ====
occ %>% group_by(Origin, Habitat, Species) %>%
  summarise(n = n()) %>% arrange(n) %>% print(n=100)

occ <- occ %>% left_join(
  occ %>% group_by(Species) %>%
    summarise(n = n())
) %>% filter(n >= 5) %>% dplyr::select(-n)


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

bromeliads <- rast(".././Variables/Bromeliads_predictions.tif")[[1]]
env_all <- c(env_all, bromeliads)

names(env_all) <- c("Elevation", "Nebulosity", "Soil", "Forest", "Agriculture", "Bromeliads")

### VIF ====
# usdm::vif(env_all)

# Variables      VIF
#
# Elevation     5.795715
# Nebulosity    3.351553
# Soil          NA
# Forest        2.251301
# Agriculture   1.777860
# Bromeliads    4.607773


# names(env_all) <- c("Elevation (m)",
#                     "Nebulosity (%)",
#                     "Soil type",
#                     "Distance to forest patch (m)",
#                     "Distance to agriculture patch (m)",
#                     "Bromeliad suitability")
# 
# par(mfrow=c(2,3))
# plot(env_all[[1]], maxcell=1E6, axes = F, col = map.pal("haxby", 100), main = names(env_all)[[1]])
# plot(env_all[[2]], maxcell=1E6, axes = F, col = map.pal("blues", 100), main = names(env_all)[[2]])
# plot(env_all[[3]], maxcell=1E6, axes = F, col = map.pal("roygbiv", 8), main = names(env_all)[[3]])
# plot(env_all[[4]], maxcell=1E6, axes = F, col = map.pal("viridis", 100), main = names(env_all)[[4]])
# plot(env_all[[5]], maxcell=1E6, axes = F, col = map.pal("viridis", 100), main = names(env_all)[[5]])
# plot(env_all[[6]], maxcell=1E6, axes = F, col = map.pal("magma", 100), main = names(env_all)[[6]])
# par(mfrow=c(1,1))


# ENM / SDM ====
## set parameters ====
# env_all <- aggregate(env_all, 20, fun = "modal") #for testing only
nb.back.pts <- 10000 # for final models, use ca. 10000 background points


## launch loop over all species ====
nc = 8 # no. of parallel stuff

done <- list.files("../Results/",pattern = "ENM_results_")
done <- gsub("ENM_results_", "", done)
done <- gsub(".csv", "", done)

for(i in setdiff(occ %>% pull(Species) %>% unique, done)){
  print(i)
  
  ### 1. select occurrence data ====
  occ.tmp <- occ %>% filter(Species == i) %>% dplyr::select(X, Y)
  occ.tmp$pr_ab = 1
  
  ##### select variables depending on habitat ====
  if(sp_char[sp_char$Species == i,"Habitat"] == "Arboreal"){
    env <- env_all[[-3]]
    categ = NULL
  } else if(sp_char[sp_char$Species == i,"Habitat"] == "Soil") {
    env <- env_all[[-6]]
    categ = "Soil"
  } else {
    env <- env_all
    categ = "Soil"
  }
  
  
  ### 2. sample background with thickening ====
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
  
  # ggplot() +
  #   geom_spatraster(data = env[[1]], na.rm = T, maxcell = 1E5) +
  #   geom_point(data = occ_back.folds %>% filter(pr_ab == 0), aes(x = X, y = Y, col = as.factor(part)),
  #              shape = 16, size = .5) +
  #   geom_point(data = occ_back.folds %>% filter(pr_ab == 1), aes(x = X, y = Y, col = as.factor(part)),
  #              shape = 22, size = 3, fill = "white") +
  #   scale_fill_distiller("Elevation", na.value = "white", palette = "Greys", guide = "none") +
  #   scale_color_brewer("Evaluation folds", palette = "Set1") +
  #   theme_void()
  
  
  ### 4. Extract environmental values ====
  occ_back_env <- sdm_extract(
    data = occ_back.folds,
    x = "X",
    y = "Y",
    env_layer = env,
    filter_na = T
  )
  
  ### 5. Tune models ====
  # metrics.used = ifelse(nrow(occ.tmp) > 10, "BOYCE", "TSS")
  metrics.used = "BOYCE"
  
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
    metric = metrics.used
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
    metric = metrics.used
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
    n.minobsinnode = c(1,5,10,20)
  )
  
  model.brt.tune <- tune_gbm(
    data = occ_back_env,
    response = "pr_ab",
    predictors = names(occ_back_env)[!names(occ_back_env) %in% c("X","Y","pr_ab","part",categ)],
    predictors_f = categ,
    partition = "part",
    grid = tune_grid.brt,
    thr = "sensitivity",
    n_cores = nc,
    metric = metrics.used
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
    metric = metrics.used
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
  
  # pred <- sdm_predict(
  #   models = list(model.net.tune, model.mx.tune, model.rf.tune),
  #   pred = env
  # )
  # model.net.tune$hyper_performance %>% select(1,2,3,15:22) %>% print(n= 40)
  # model.net.tune$performance %>% select(1,2,3,20:27)
  # 
  # for(i in 1:4){
  #   
  #   partition = i
  #   
  #   mod <- nnet::nnet(pr_ab ~ Elevation+Nebulosity+Forest+Agriculture+Bromeliads, 
  #                     data = occ_back_env %>% filter(part != partition), 
  #                     size = model.net.tune$performance$size, rang = 0.1, 
  #                     decay = model.net.tune$performance$decay, maxit = 200, 
  #                     trace = FALSE)
  #   
  #   pred.net <- cbind(pred = predict(mod, occ_back_env %>% filter(part == partition))[,1],
  #                     pr_ab = occ_back_env %>% filter(part == partition) %>% pull(pr_ab)) %>% as.data.frame
  #   
  #   
  #   print(pROC::auc(response = pred.net$pr_ab,
  #                   predictor = pred.net$pred,
  #                   plot = T))
  #   print(ecospat::ecospat.boyce(fit = pred.net[pred.net$pr_ab == 0, "pred"],
  #                                obs = pred.net[pred.net$pr_ab == 1, "pred"],
  #                                method = "kendall")$cor
  #   )
  #   
  #   
  # }
  
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
      min(ceiling(nrow(occ.tmp)/3),6)
    ),
    subset = dc("Elevation", "I(Elevation^2)") && dc("Forest", "I(Forest^2)") &&
      dc("Bromeliads", "I(Bromeliads^2)") && dc("Agriculture", "I(Agriculture^2)")  && 
      dc("Nebulosity", "I(Nebulosity^2)")
  )
  
  stopCluster(clust)
  
  model.glm.3 <- MuMIn::get.models(table.mod.glm, subset = 1)[[1]]
  
  # var.glm.3 <- attr(model.glm.3$terms, 'variables')
  # var.glm.3 <- lapply(3:length(var.glm.3),function(x)var.glm.3[[x]])
  # var.glm.3.quad <- var.glm.3[grepl("I",var.glm.3)]
  # var.glm.3.quad <- sapply(var.glm.3.quad, `[[`, 2)
  # var.glm.3.quad <- sapply(var.glm.3.quad, `[[`, 2)
  # var.glm.3.lin <- var.glm.3[!grepl("I",var.glm.3)]
  # missing.lin <- setdiff(var.glm.3.quad, var.glm.3.lin)
  # 
  # if(length(missing.lin)>0){
  #   model.glm.4 <- update(model.glm.3, paste("~.+", missing.lin))
  # }else{
  #   model.glm.4 <- model.glm.3
  # }
  
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
    thr = "sensitivity"
  )
  
  model.gam.2 <- mgcv::bam(model.gam$model$formula,
                           data = occ_back_env,
                           family = 'binomial', discrete = TRUE, na.action = na.fail)
  
  table.mod.gam <- MuMIn::dredge(
    model.gam.2,
    m.lim= c(
      1,
      min(ceiling(nrow(occ.tmp)/3),6)
    )
  )
  model.gam.3 <- MuMIn::get.models(table.mod.gam, subset = 1)[[1]]
  
  try(
    model.gam.4 <- fit_gam(
      data = occ_back_env,
      response = "pr_ab",
      predictors = gratia::term_names(model.gam.3),
      predictors_f = if(any(grepl("Soil", gratia::term_names(model.gam.3)))){categ}else{NULL},
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
      bind_cols(model = "gam", parameter = paste("formula =",as.character(model.gam.4$model$formula)[[3]])),
      bind_cols(model = "glm", parameter = paste("formula =",as.character(model.glm.4$model$formula)[[3]])),
    ) %>% 
    group_by(model) %>% 
    summarise(parameters = paste(parameter, collapse = " ; "))
  
  eval_results <- eval_results %>% left_join(tune_results) %>% select(-c(2:5))
  
  
  ### 7. Ensemble modelling ====
  preds <- c(pred.rf$raf, pred.mx$max, pred.brt$gbm, pred.net$net, pred.glm$glm, pred.gam$gam)
  weights <- eval_results %>% pull(paste0(metrics.used, "_mean"))
  # weights.std <- weights[weights > 0] * length(weights[weights > 0])/sum(weights[weights > 0])
  
  preds <- preds[[weights > 0]]
  pred.ens.mean <- weighted.mean(preds, weights[weights > 0], na.rm = T)
  
  # pred.ens.var <- app(preds, function(x){Hmisc::wtd.var(x, weights.std, normwt=F)})
  # pred.ens.sd <- sqrt(pred.ens.var)
  # pred.ens.sd <- mask(pred.ens.sd, pred.ens.mean)
  
  pred.ens.ca.mtp  <- sum(c(pred.rf$mtp,
                            pred.mx$mtp,
                            pred.brt$mtp,
                            pred.net$mtp,
                            pred.glm$mtp,
                            pred.gam$mtp))
  pred.ens.ca.p10  <- sum(c(pred.rf$p10,
                            pred.mx$p10,
                            pred.brt$p10,
                            pred.net$p10,
                            pred.glm$p10,
                            pred.gam$p10))

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
  pred.final[is.na(pred.final) & !is.na(env[[1]])] <- 0
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
    "GAM" = model.gam.4
  )
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
        mutate(Algorithm = names(list.mods)[[ii]], Species = i)
    )
  }
  
  ### 9. Write everything ====
  # rasters
  writeRaster(pred.final, paste0("../Results/ENM_predictions_",i,".tif"), overwrite = T)
  write_csv(eval_results %>% mutate(Species = i, .before = 1), paste0("../Results/ENM_results_",i,".csv"))
  write_csv(linear_predictions %>% mutate(Species = i, .before = 1), paste0("../Results/ENM_responseCurves_",i,".csv"))
}




## test plot
linear_predictions %>% group_by(Algorithm) %>% mutate(y = range01(y)) %>%
  ggplot(aes(x=x, y = y, col = Algorithm))+
  facet_wrap(~Variable, scale = 'free_x', strip.position = "bottom") +
  geom_line() + 
  scale_color_brewer(palette = "Set2") +
  theme_minimal()

