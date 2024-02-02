# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# ANALYSES
# ================================= #

# Load libraries ====
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(ecospat)
library(FactoMineR)
# library(ENMeval)
library(ENMTools)
# library(retry)

# set aggregation parameters ====
max.c = 10000
aggr.fact = 10

# reload all data ====
## rasters ====
pred_all <- rast("../SDM_predictions.tif")
pred.bin_all <- rast("../SDM_predictions.binary.tif")

## eval / importance / curves ====
Var_importance <- read_csv(".././Var_importance.csv")
SDM_evaluations <- read_csv(".././SDM_evaluations.csv")
Var_response <- read_csv(".././Var_response.csv")

## occurrences
occ <- read_csv("../species_records.csv") %>% filter(Y<15, Y>14, X < -60.8) %>% 
  distinct()
sp.char <- occ %>% group_by(Species) %>% 
  summarise(Origin = unique(Origin), Habitat = unique(Habitat), n.occ = n())
sp.char %>% arrange(Habitat, Origin) %>% knitr::kable(.)

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

# Check models ====
## evaluation ====
SDM_evaluations %>% 
  mutate(metric.eval = gsub("BOYCE", "Boyce index", metric.eval),
         metric.eval = gsub("ROC", "AUC", metric.eval),
         metric.eval = factor(metric.eval, levels = c("Boyce index", "TSS", "AUC"))) %>% 
  left_join(sp.char) %>% 
  ggplot(aes(y = mean, x = Origin, fill = Origin)) +
  geom_boxplot() +
  facet_grid(metric.eval ~ Habitat, scales = "free") +
  scale_y_continuous("Evaluation index", expand = c(0,0)) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white")
  )

ggsave(filename = "../eval.pdf", height = 4, width = 5)


## variable importance ====
Var_importance %>% 
  left_join(sp.char) %>% 
  ggplot(aes(y = mean, x = Variable, fill = Origin)) +
  facet_grid(. ~ Habitat, scales = "free") +
  geom_boxplot() +
  scale_fill_brewer(palette="Set1") +
  scale_y_continuous("Variable importance", expand = c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.9,.8),
        legend.background = element_rect(color = "grey"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white")
  )

## response to variables ====
Var_response %>% 
  left_join(sp.char) %>% 
  ggplot(aes(x = expl.val, y = pred.val.mean, col = Algorithm, fill = Algorithm)) +
  facet_grid(Species~expl.name, scale = "free") + 
  geom_line() +
  geom_ribbon(aes(ymin = pred.val.mean - pred.val.se, ymax = pred.val.mean + pred.val.se), 
              alpha = .3, color = NA) +
  theme_classic()


# Plot model predictions ====
## continuous ====
ggplot() +
  geom_spatraster(data = pred_all/1000, maxcell = max.c) +
  facet_wrap(~ lyr) +
  scale_fill_viridis_c("Suitability", na.value = "transparent") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold.italic", color = "white")
  ) +
  coord_sf()

ggsave(filename = "../Suitability_maps.pdf", height = 7, width = 9)
ggsave(filename = "Suitability_maps.png", height = 7, width = 9)

## binary ====
ggplot() +
  geom_spatraster(data = as.factor(pred.bin_all), maxcell = max.c) +
  facet_wrap(~ lyr) +
  scale_fill_manual(
    "Suitability", 
    values = c("lightgrey", "orange3"), 
    labels = c("Unsuitable", "Suitable", ""), 
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold.italic", color = "white")
  ) +
  coord_sf()

ggsave(filename = "Suitability_maps_binary.png", height = 7, width = 9.5)


# Overlap between introduced and native species ====
## all species ====
n.over <- nrow(unique(sum(pred.bin_all)))
ggplot() +
  geom_spatraster(data = as.factor(sum(pred.bin_all)), maxcell = max.c) +
  scale_fill_manual(
    "No. species",
    values = c("grey85", colorRampPalette(c("khaki1", "red4"))(n.over-1)),
    labels = c(0:(n.over-1),""),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
  ) +
  coord_sf()

## introduced and native ====
as.data.frame(pred.bin_all %>% aggregate(.,aggr.fact, fun = max), xy = T) %>% as_tibble %>% 
  pivot_longer(cols = -c(1:2), names_to = "Species", values_to = "Suitability") %>% 
  left_join(sp.char) %>% group_by(x, y, Origin, Habitat) %>% 
  summarise(Overlap = sum(Suitability)) %>% 
  mutate(Overlap = factor(Overlap, levels = 0:max(Overlap))) %>%
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  facet_grid(Habitat~Origin) +
  geom_raster() +
  scale_fill_manual(
    "No. species",
    values = c("grey85", colorRampPalette(c("khaki1", "red4"))(4)),
    labels = c(0:(4),""),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold", color = "white")
  ) +
  coord_sf(crs = st_crs("EPSG:4326"))

ggsave(filename = "Species_richness_separated.png", height = 7, width = 9.5)

## introduced vs. native ====
overlap.plot <- c()
for(i in sp.char %>% filter(Origin == "Native") %>% pull("Species")){
  focal.sp_rast <- subset(pred.bin_all, i)
  corr.sps <- sp.char %>% 
    filter(Origin == "Introduced", Habitat == sp.char[sp.char$Species == i, "Habitat"][[1]]) %>% 
    pull(Species)
  for(j in corr.sps){
    corr.sps_rast <- subset(pred.bin_all, c(i,j))
    overlap.plot.tmp <- as.data.frame(aggregate(corr.sps_rast, aggr.fact, fun = max), xy = T) %>% 
      as_tibble %>% 
      rename(Introduced = 3, Native = 4) %>% 
      mutate(Overlap = ifelse(Native == 1 & Introduced == 1, "Overlap",
                              ifelse(Native == 1 & Introduced == 0, "Native only",
                                     ifelse(Native == 0 & Introduced == 1, "Introduced only", 
                                            "No species")))) %>% 
      mutate(Introduced = j, Native = i)
    overlap.plot <- bind_rows(overlap.plot, overlap.plot.tmp %>% 
                                mutate(Habitat = sp.char[sp.char$Species == i, "Habitat"][[1]]))
  }
}

### plot a matrix of overlap area ====
overlap.plot %>% group_by(Introduced, Native, Habitat) %>% 
  summarise(n.overlap = length(Overlap[Overlap == "Overlap"])) %>% 
  ggplot(aes(x = Native, y = Introduced, 
             fill = n.overlap, 
             label = n.overlap)) +
  facet_wrap(~Habitat, scale = "free", shrink  =F) +
  geom_tile(color = "white") +
  geom_text(size = 3, col = "white") +
  scale_fill_gradient(low = "lightgrey", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic"),
        axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white"),
        axis.line = element_blank(),
        legend.pos = "none")


### plot overal areas ====
#### Soil
overlap.plot %>% filter(Habitat == "Soil") %>% 
  mutate(Overlap = factor(Overlap, levels = c("No species", "Native only", "Introduced only", "Overlap"))) %>%
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  facet_grid(Native ~ Introduced) +
  geom_raster() +
  scale_fill_manual(
    "",
    values = c("lightgrey", "green4", "yellow3", "red3")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold.italic", color = "white")
  ) +
  coord_sf(crs = st_crs("EPSG:4326"))

ggsave(filename = "Pairwise_overlap_soil.png", height = 4, width = 8)

overlap.plot %>% filter(Habitat == "Arboreal") %>% 
  mutate(Overlap = factor(Overlap, levels = c("No species", "Native only", "Introduced only", "Overlap"))) %>%
  ggplot(aes(x = x, y = y, fill = Overlap)) +
  facet_grid(Native ~ Introduced) +
  geom_raster() +
  scale_fill_manual(
    "",
    values = c("lightgrey", "green4", "yellow3", "red3")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(face = "bold.italic", color = "white")
  ) +
  coord_sf(crs = st_crs("EPSG:4326"))

ggsave(filename = "Pairwise_overlap_arboreal.png", height = 6, width = 7)

### synthetic invasion risk ====
overlap.plot %>% group_by(x, y, Habitat) %>% 
  summarise(n.overlap = length(Overlap[Overlap == "Overlap"])) %>% 
  ggplot(aes(x = x, y = y, fill = as.factor(n.overlap))) +
  facet_grid(~Habitat) +
  geom_raster() +
  scale_fill_viridis_d(
    "No. pairs of\noverlapping species",
    option = "magma"
  ) +
  theme_void() +
  theme(strip.text = element_text(face = 'bold')) +
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
nb.back.pts <- 5000
back.env <- spatSample(env, na.rm = T, size = nb.back.pts)

## overlap in environmental space (D_env) ====
PCA.cal <- PCA(
  bind_rows(back.env %>% na.omit, 
            occ_env %>% dplyr::select(colnames(back.env)) %>% na.omit),
  ind.sup = (nrow(back.env %>% na.omit)+1):
    (nrow(back.env %>% na.omit)+nrow(occ_env %>% dplyr::select(colnames(back.env)) %>% na.omit)),
  ncp = 2,
  graph = F
)

scores.back <- PCA.cal$ind$coord[1:nrow(back.env %>% na.omit),]
scores.sp <- PCA.cal$ind.sup$coord %>% as_tibble %>% 
  mutate(Species = occ_env %>% na.omit %>% pull(Species))

ggplot() + 
  geom_point(data = as_tibble(scores.back), aes(x = Dim.1, y = Dim.2), color = "lightgrey") +
  geom_point(data = scores.sp, aes(x = Dim.1, y = Dim.2, color = Species)) +
  theme_bw()

overlap_env <- c()
for(i in sp.char %>% filter(Origin == "Introduced") %>% pull(Species) %>% unique){
  for(j in sp.char %>% filter(Origin == "Native") %>% pull(Species) %>% unique){
    if(i != j){
      scores.sp1 <- scores.sp %>% filter(Species == i)
      scores.sp2 <- scores.sp %>% filter(Species == j)
      
      boolFalse<-F
      tryCatch({
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
        
        boolFalse<-T
      },error=function(e){
      },finally={})
      
      if(boolFalse == F){
        D <- 0
      }
      
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


p_d_env <- overlap_env %>% left_join(sp.char, by = c("Introduced" = "Species")) %>% 
  mutate(Introduced = factor(Introduced, levels = sp.char[order(sp.char$Habitat), "Species"]$Species)) %>% 
  ggplot(aes(y=Introduced, x=Native, fill=D_env, label = round(D_env,2))) +
  geom_tile() +
  geom_text(size = 2.5, col = "white") +
  scale_y_discrete("Introduced species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[env])), low = "green3", high = "dark red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_fixed() +
  facet_grid(~ Habitat) +
  theme(strip.text = element_text(face = 'bold'),
        strip.background = element_blank())


## overlap in geographical space (D_geo) ====

overlap_geo <- c()
for(i in sp.char %>% filter(Origin == "Introduced") %>% pull(Species) %>% unique){
  for(j in sp.char %>% filter(Origin == "Native") %>% pull(Species) %>% unique){
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

p_d_geo <-  overlap_geo %>% left_join(sp.char, by = c("Introduced" = "Species")) %>% 
  mutate(Introduced = factor(Introduced, levels = sp.char[order(sp.char$Habitat), "Species"]$Species)) %>% 
  ggplot(aes(y=Introduced, x=Native, fill=D_geo, label = round(D_geo,2))) +
  geom_tile() +
  geom_text(size = 2.5, col = "white") +
  scale_y_discrete("Introduced species") +
  scale_x_discrete("Native species") +
  scale_fill_gradient(bquote(italic("D"[geo])), low = "blue", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_fixed() +
  facet_grid(~ Habitat) +
  theme(strip.text = element_text(face = 'bold'),
        strip.background = element_blank())

## plot ====
cowplot::plot_grid(p_d_env, p_d_geo, nrow = 2, labels = "AUTO")

ggsave(filename = "overlap_indices.png", width = 8, height = 10)

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
