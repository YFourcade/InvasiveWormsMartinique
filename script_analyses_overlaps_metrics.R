# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# ANALYSES - Metrics of overlap
# ================================= #

# Load libraries ====
library(tidyverse)
library(terra)
library(sf)
library(ecospat)
library(FactoMineR)
library(readxl)
library(factoextra)
library(ENMTools)
library(hypervolume)
library(ggrepel)
# library(patchwork)

# reload all data ====
## rasters ====
# continuous
pred_all <- lapply(
  list.files("../Results/", pattern = "ENM_predictions", full.names = T),
  function(x){
    rast(x)[["Weighted mean"]]
  }
) %>% rast

names(pred_all) <- list.files("../Results/", pattern = "ENM_predictions")
names(pred_all) <- gsub("ENM_predictions_", "", names(pred_all))
names(pred_all) <- gsub(".tif", "", names(pred_all))

## occurrences ====
sp_char <- read_excel(".././SELECTION_SP.xlsx", sheet = 2) %>% filter(Species %in% names(pred_all))
occ <- list.files(".././Data", pattern = "OK", full.names = T) %>%
  lapply(.,read_excel) %>% bind_rows() %>% dplyr::select(nom_espece, coord_X, coord_Y) %>%
  rename(Species = 1, X = coord_Y, Y = coord_X) %>%
  mutate(Species = stringr::word(Species, 1, 2),
         Species = ifelse(Species == "Dichogaster sp11", "Dichogaster sp11A", Species),
         Species = ifelse(Species == "Glossodrilus sp. nov. 3", "Glossodrilus spp", Species),
         Species = ifelse(Species == "Glossodrilus sp.", "Glossodrilus sp1", Species)) %>%
  distinct %>% filter(X > -61.5)%>% 
  left_join(sp_char) %>% 
  group_by(Species) %>%
  mutate(n = n()) %>% filter(n >= 5) %>% dplyr::select(-n)

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



# Metrics of geographical and environmental overlap ====
## extract environmental values ====
### at occurrences ====
occ_env <- bind_cols(
  Species = occ[,c(1,6,7)], 
  terra::extract(env_all, occ[,3:2])[,-1]
) %>% ungroup %>% na.omit

### at background points ====
nb.back.pts <- 10000
# back.env <- spatSample(env_all, na.rm = T, size = nb.back.pts)
# write_csv(back.env, "../back.env.csv")
back.env <- read_csv("../back.env.csv")

### all points - Arboreal ====
allpts.arb <- bind_rows(
  background = back.env %>% dplyr::select(-Soil), 
  occurrence = occ_env %>% filter(grepl("Arb", Habitat)) %>%  dplyr::select(colnames(back.env), -Soil),
  .id = "type"
)

### all points - Soil ====
allpts.soil <- bind_rows(
  background = back.env %>% dplyr::select(-Bromeliads), 
  occurrence = occ_env %>% filter(grepl("Soil", Habitat)) %>%  dplyr::select(colnames(back.env), -Bromeliads),
  .id = "type"
)


## overlap in environmental space ====
### Arboreal ====
PCA.cal_Arb <- PCA(
  X = allpts.arb[,-1],
  ind.sup = (length(allpts.arb$type[allpts.arb$type == "background"])+1):
    nrow(allpts.arb)
)

scores.back_Arb <- PCA.cal_Arb$ind$coord
scores.sp_Arb <- PCA.cal_Arb$ind.sup$coord %>% as_tibble %>% 
  mutate(Species = occ_env %>% filter(grepl("Arb", Habitat)) %>% pull(Species),
         Origin = occ_env %>% filter(grepl("Arb", Habitat)) %>% pull(Origin),
         Habitat = occ_env %>% filter(grepl("Arb", Habitat)) %>% pull(Habitat))


p.pca.var <- fviz_pca_var(PCA.cal_Arb, repel = TRUE, geom = c("arrow", "text"))

p.pca <- ggplot() + 
  geom_point(data = as_tibble(scores.back_Arb), 
             aes(x = Dim.1, y = Dim.2), color = "lightgrey", size = .5) +
  geom_count(data = scores.sp_Arb, 
             aes(x = Dim.1, y = Dim.2, color = Origin, shape = Habitat, size = after_stat(prop)),
             alpha = .8) +
  stat_ellipse(data = scores.sp_Arb, 
               aes(x = Dim.1, y = Dim.2, color = Origin), type = 'norm') +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(paste0("Dim 2 (", round(PCA.cal_Arb$eig[2,2],2), "%)")) +
  scale_x_continuous(paste0("Dim 1 (", round(PCA.cal_Arb$eig[1,2],2), "%)")) +
  scale_shape_manual(values = c(15,16)) +
  theme_bw()


### Soil ====
FAMD.cal_Soil <- FAMD(
  base = allpts.soil[,-1],
  ind.sup = (length(allpts.soil$type[allpts.soil$type == "background"])+1):
    nrow(allpts.soil)
)

scores.back_Soil <- FAMD.cal_Soil$ind$coord
scores.sp_Soil <- FAMD.cal_Soil$ind.sup$coord %>% as_tibble %>% 
  mutate(Species = occ_env %>% filter(grepl("Soil", Habitat)) %>% pull(Species),
         Origin = occ_env %>% filter(grepl("Soil", Habitat)) %>% pull(Origin),
         Habitat = occ_env %>% filter(grepl("Soil", Habitat)) %>% pull(Habitat))

p.famd.var <- fviz_famd_var(FAMD.cal_Soil, repel = TRUE)
fviz_famd_var(FAMD.cal_Soil, "quanti.var", repel = TRUE,
              col.var = "black")
fviz_famd_var(FAMD.cal_Soil, "quali.var", repel = TRUE,
              col.var = "black")

fviz_famd_ind(FAMD.cal_Soil,
              geom = c("point"))


p.famd <- ggplot() + 
  geom_point(data = as_tibble(scores.back_Soil), 
             aes(x = Dim.1, y = Dim.2), color = "lightgrey", size = .5) +
  geom_count(data = scores.sp_Soil,
             aes(x = Dim.1, y = Dim.2, color = Origin, shape = Habitat, size = after_stat(prop)),
             alpha = .8) +
  stat_ellipse(data = scores.sp_Soil, 
               aes(x = Dim.1, y = Dim.2, color = Origin), type = 'norm') +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(paste0("Dim 2 (", round(FAMD.cal_Soil$eig[2,2],2), "%)")) +
  scale_x_continuous(paste0("Dim 1 (", round(FAMD.cal_Soil$eig[1,2],2), "%)")) +
  scale_shape_manual(values = c(17,16)) +
  theme_bw()


cowplot::plot_grid(p.pca, p.famd, align = "hv", axis = "tblr")

ggsave(filename = ".././Figures/FigPCA_raw.pdf", width = 9, height = 3)


## D_env ====
overlap_env <- c()
for(i in sp_char %>% filter(Origin == "Native") %>% pull(Species) %>% unique){
  native.sp.tmp.habitat <- sp_char %>% filter(Species == i) %>% pull(Habitat)
  exotic.sp.tmp <- sp_char %>% 
    filter(Origin == "Exotic", grepl(native.sp.tmp.habitat, Habitat)) %>% 
    pull(Species)
  
  for(j in exotic.sp.tmp){
    if(native.sp.tmp.habitat == "Arboreal"){
      scores.back <- scores.back_Arb
      scores.sp1 <- scores.sp_Arb %>% filter(Species == i)
      scores.sp2 <- scores.sp_Arb %>% filter(Species == j)
    } else {
      scores.back <- scores.back_Soil
      scores.sp1 <- scores.sp_Soil %>% filter(Species == i)
      scores.sp2 <- scores.sp_Soil %>% filter(Species == j)
    }
    
    boolFalse<-F
    tryCatch({
      sp1_ecospat <- ecospat.grid.clim.dyn(
        glob = scores.back[,1:2],
        glob1 = scores.back[,1:2],
        sp = scores.sp1[,1:2],
        R = 100,
        th.sp = NULL, 
        th.env = NULL
      )
      
      sp2_ecospat <- ecospat.grid.clim.dyn(
        glob = scores.back[,1:2],
        glob1 = scores.back[,1:2],
        sp = scores.sp2[,1:2],
        R = 100,
        th.sp = NULL, 
        th.env = NULL
      )
      
      D <- ecospat.niche.overlap(sp1_ecospat, sp2_ecospat, cor = F)$D
      
      # niche.test <- ecospat.niche.similarity.test(
      #   sp1_ecospat, sp2_ecospat, 
      #   rep = 100, rand.type = 1, ncores = 6
      # )
      
      boolFalse<-T
    },error=function(e){
    },finally={})
    
    if(boolFalse == F){
      D <- 0
    }
    
    overlap_env <- rbind.data.frame(
      overlap_env,
      cbind.data.frame(
        Native = i , 
        Exotic = j, 
        D_env = D
      )
    )
  }
}

overlap_env <- overlap_env %>% 
  left_join(sp_char %>% filter(Origin == "Native") %>% 
              dplyr::select(Species, Habitat),
            by = c("Native" = "Species"))


p_d_env <- ggplot(overlap_env, aes(y=Exotic, x=Native, fill=D_env, label = round(D_env,2))) +
  facet_wrap(~ Habitat, scales = "free") +
  geom_tile() +
  geom_text(size = 3, col = "white", fontface = "bold") +
  scale_y_discrete("Exotic species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[env])), low = "chartreuse3", high = "coral") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic"),
        strip.text = element_text(face = 'bold'),
        strip.background = element_blank(),
        axis.text.y = element_text(face = "italic"))


## overlap in geographical space (D_geo) ====

overlap_geo <- c()
for(i in sp_char %>% filter(Origin == "Native") %>% pull(Species) %>% unique){
  native.sp.tmp.habitat <- sp_char %>% filter(Species == i) %>% pull(Habitat)
  exotic.sp.tmp <- sp_char %>% 
    filter(Origin == "Exotic", grepl(native.sp.tmp.habitat, Habitat)) %>% 
    pull(Species)
  
  for(j in exotic.sp.tmp){
    sdm1 <- pred_all[[names(pred_all) == i]]
    sdm2 <- pred_all[[names(pred_all) == j]]
    D <- raster.overlap(sdm1, sdm2)$D
    
    overlap_geo <- rbind.data.frame(
      overlap_geo,
      cbind.data.frame(
        Native = i , 
        Exotic = j, 
        D_geo = D
      ) 
    )
  }
}


overlap_geo <- overlap_geo %>% 
  left_join(sp_char %>% filter(Origin == "Native") %>% 
              dplyr::select(Species, Habitat),
            by = c("Native" = "Species"))


p_d_geo <- ggplot(overlap_geo, aes(y=Exotic, x=Native, fill=D_geo, label = round(D_geo,2))) +
  facet_wrap(~ Habitat, scales = "free") +
  geom_tile() +
  geom_text(size = 3, col = "white", fontface = "bold") +
  scale_y_discrete("Exotic species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[geo])), low = "gold1", high = "royalblue1") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic"),
        strip.text = element_text(face = 'bold'),
        strip.background = element_blank(),
        axis.text.y = element_text(face = "italic"))

## plot ====
cowplot::plot_grid(p_d_env, p_d_geo, nrow = 2, 
                   labels = "AUTO")

ggsave(filename = ".././Figures/Fig3.pdf", width = 9, height = 8)
ggsave(filename = ".././Figures/Fig3.tiff", width = 9, height = 8, dpi = 300)

# merge
overlap <- left_join(overlap_env, overlap_geo)

ggplot(overlap, aes(y = D_geo, x = D_env)) +
  geom_point() +
  facet_grid(~Habitat) +
  scale_y_continuous(bquote(italic("D"[geo]))) +
  scale_x_continuous(bquote(italic("D"[env]))) +
  theme_bw() +
  theme(strip.text = element_text(face = 'bold'),
        strip.background = element_blank(),
        panel.grid = element_blank())

# ggsave(filename = ".././Figures/Fig_niche.breadth.svg", width = 6, height = 4)


overlap %>% group_by(Habitat) %>% nest() %>% 
  mutate(test = map(data, ~ cor.test(.x$D_env, .x$D_geo, method = "spearman")),
         tidied = map(test, broom::tidy)
  ) %>% 
  unnest(tidied)

wilcox.test(overlap$D_env, overlap$D_geo)


overlap %>% summarise(D_geo = mean(D_geo), D_env = mean(D_env))

# Habitat   data              test    estimate statistic    p.value method                          alternative
# 1 Soil     <tibble [30 × 4]> <htest>    0.732      1204 0.00000831 Spearman's rank correlation rho two.sided  
# 2 Arboreal <tibble [15 × 4]> <htest>    0.829        96 0.000188   Spearman's rank correlation rho two.sided  


# Metrics of geographical and environmental niche width ====
allpts.allEnv <- bind_rows(
  background = back.env, 
  occurrence = occ_env %>% filter(Habitat == "Soil / Arboreal") %>% dplyr::select(colnames(back.env)),
  .id = "type"
)

FAMD.cal_All <- FAMD(
  base = allpts.allEnv[,-1],
  ind.sup = (length(allpts.allEnv$type[allpts.allEnv$type == "background"])+1):
    nrow(allpts.allEnv)
)

scores.sp_allEnv <- FAMD.cal_All$ind.sup$coord %>% as_tibble %>% 
  mutate(Species = occ_env %>% filter(Habitat == "Soil / Arboreal") %>% pull(Species),
         Origin = occ_env %>% filter(Habitat == "Soil / Arboreal") %>% pull(Origin),
         Habitat = occ_env %>% filter(Habitat == "Soil / Arboreal") %>% pull(Habitat))


niche.breadth <- c()
for(i in names(pred_all)){
  pred_temp <- rast(paste0(".././Results/ENM_predictions_", i, ".tif"))[[5]]
  sdm_breadth <- as.data.frame(pred_temp) %>% table %>% .["1"] / global(pred_temp, "notNA")[[1]]
  
  if((sp_char %>% filter(Species == i) %>% pull(Habitat)) == "Arboreal"){
    env.tmp <-  scores.sp_Arb %>% filter(Species == i) %>% dplyr::select(1:2)
  } else if((sp_char %>% filter(Species == i) %>% pull(Habitat)) == "Soil"){
    env.tmp <-  scores.sp_Soil %>% filter(Species == i) %>% dplyr::select(1:2)
  } else {
    env.tmp <-  scores.sp_allEnv %>% filter(Species == i) %>% dplyr::select(1:2)
  }
  
  hyper <- hypervolume(env.tmp, method = 'svm')
  env_breadth <- get_volume(hyper)
  
  niche.breadth <- rbind.data.frame(
    niche.breadth,
    cbind.data.frame(
      Species = i , 
      Niche_breadth_geo = sdm_breadth, 
      Niche_breadth_env = env_breadth
    ) 
  )
  
}

niche.breadth %>% left_join(sp_char) %>% ungroup %>% 
  mutate(Niche_breadth_geo = scale(Niche_breadth_geo),
         Niche_breadth_env = scale(Niche_breadth_env)) %>% 
ggplot(aes(y = Niche_breadth_geo, x = Niche_breadth_env)) +
  # geom_smooth(method = 'lm', se = F, color = 'grey50', linewidth = .5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(color = Origin, shape = Habitat)) +
  geom_text_repel(aes(color = Origin, label = Species), fontface = "italic", size = 3, min.segment.length = .50) + 
  scale_y_continuous("Niche breadth in the geographical space (scaled)") +
  scale_x_continuous("Niche breadth in the environmental space (scaled)") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(panel.grid = element_blank(), aspect.ratio=1)

niche.breadth %>%  left_join(sp_char) %>% 
  pivot_longer(cols = c(2,3), names_to = "type", values_to = "breadth") %>% 
  mutate(type = gsub("Niche_breadth_", "", type)) %>% 
  group_by(type) %>% 
  group_modify(~ broom::tidy(wilcox.test(breadth ~ Origin, data = .x)))


niche.breadth %>%  left_join(sp_char) %>% 
  pivot_longer(cols = c(2,3), names_to = "type", values_to = "breadth") %>% 
  mutate(type = gsub("Niche_breadth_", "", type)) %>% 
  group_by(type) %>% 
  group_modify(~ broom::tidy(wilcox.test(breadth ~ Origin, data = .x)))

niche.breadth %>%  left_join(sp_char) %>% 
  pivot_longer(cols = c(2,3), names_to = "type", values_to = "breadth") %>% 
  filter(Habitat %in% c("Soil", "Arboreal")) %>% 
  mutate(type = gsub("Niche_breadth_", "", type)) %>% 
  group_by(type) %>% 
  group_modify(~ broom::tidy(wilcox.test(breadth ~ Habitat, data = .x)))

niche.breadth %>%  left_join(sp_char) %>% 
  pivot_longer(cols = c(2,3), names_to = "type", values_to = "breadth") %>% 
  filter(Habitat %in% c("Soil", "Soil / Arboreal")) %>% 
  mutate(type = gsub("Niche_breadth_", "", type)) %>% 
  group_by(type) %>% 
  group_modify(~ broom::tidy(wilcox.test(breadth ~ Habitat, data = .x)))

niche.breadth %>%  left_join(sp_char) %>% 
  pivot_longer(cols = c(2,3), names_to = "type", values_to = "breadth") %>% 
  filter(Habitat %in% c("Arboreal", "Soil / Arboreal")) %>% 
  group_by(type) %>% 
  group_modify(~ broom::tidy(wilcox.test(breadth ~ Habitat, data = .x)))


niche.breadth %>% left_join(sp_char) %>% 
  group_by(Soil, Arboreal) %>% 
  summarise(Niche_breadth_env = mean(Niche_breadth_env),
            Niche_breadth_geo = mean(Niche_breadth_geo))


cor.test(~ Niche_breadth_geo + Niche_breadth_env, data = niche.breadth, method = "spearman")

lm(Niche_breadth_geo ~ Niche_breadth_env, data = niche.breadth) %>% 
  parameters::parameters(standardize = "refit")

lm(scale(Niche_breadth_geo) ~ scale(Niche_breadth_env), data = niche.breadth) %>% 
  performance::model_performance()

ggsave(filename = ".././Figures/Fig2.pdf", width = 7, height = 6)
ggsave(filename = ".././Figures/Fig2.tiff", width = 7, height = 6, dpi = 300)

