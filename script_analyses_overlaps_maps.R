# ================================= #
# Ecological modelling of native
# and introduced earthworm species
# in Martinique
#
# ANALYSES
# ================================= #
source("functions_sdm.R")

# Load libraries ====
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(readxl)
library(rnaturalearth)
library(patchwork)
library(lmerTest)

# set aggregation parameters ====
max.c = 1E6

# reload all data ====
## rasters ====
# binary
pred_all_bin <- lapply(
  list.files("../Results/", pattern = "ENM_predictions", full.names = T),
  function(x){
    rast(x)[["Binary ensemble (P10)"]]
  }
) %>% rast

names(pred_all_bin) <- list.files("../Results/", pattern = "ENM_predictions")
names(pred_all_bin) <- gsub("ENM_predictions_", "", names(pred_all_bin))
names(pred_all_bin) <- gsub(".tif", "", names(pred_all_bin))

## eval / importance / curves ====
SDM_evaluations <- lapply(
  list.files("../Results/", pattern = "ENM_results", full.names = T),
  read_csv
) %>% bind_rows %>% 
  left_join(sp_char)

Var_response <- lapply(
  list.files("../Results/", pattern = "response", full.names = T),
  read_csv
) %>% bind_rows %>% 
  mutate(Variable = gsub("Bromeliads", "Bromeliad suitability", Variable),
         Variable = gsub("Agriculture", "Distance to agriculture patch (m)", Variable),
         Variable = gsub("Forest", "Distance to forest patch (m)", Variable),
         Variable = gsub("Nebulosity", "Nebulosity (%)", Variable),
         Variable = gsub("Elevation", "Elevation (m)", Variable))

## occurrences ====
sp_char <- read_excel(".././SELECTION_SP.xlsx", sheet = 2)
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

occ.nb <- occ %>% group_by(Species) %>% 
  summarise(n.occ = n())

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


# Check models ====
## evaluation ====
SDM_evaluations %>% dplyr::select(Species, model, BOYCE_mean, AUC_mean, TSS_mean) %>% 
  pivot_longer(-c(1,2), names_to = "metric.eval") %>% 
  mutate(metric.eval = gsub("BOYCE_mean", "Boyce index", metric.eval),
         metric.eval = gsub("AUC_mean", "AUC", metric.eval),
         metric.eval = gsub("TSS_mean", "TSS", metric.eval),
         metric.eval = factor(metric.eval, levels = c("Boyce index", "TSS", "AUC"))) %>% 
  group_by(Species,metric.eval) %>% 
  summarise(value = mean(value)) %>% 
  left_join(sp_char) %>% 
  ggplot(aes(y = value, x = Origin, fill = Origin)) +
  geom_boxplot() +
  facet_grid(metric.eval ~ Habitat, scales = "free", space = "free_x") +
  scale_y_continuous("Evaluation index") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white")
  )

SDM_evaluations %>% dplyr::select(Species, model, BOYCE_mean, AUC_mean, TSS_mean) %>% 
  pivot_longer(-c(1,2), names_to = "metric.eval") %>% 
  mutate(metric.eval = gsub("BOYCE_mean", "Boyce index", metric.eval),
         metric.eval = gsub("AUC_mean", "AUC", metric.eval),
         metric.eval = gsub("TSS_mean", "TSS", metric.eval),
         metric.eval = factor(metric.eval, levels = c("Boyce index", "TSS", "AUC"))) %>% 
  left_join(sp_char) %>% 
  group_by(Habitat, Origin, metric.eval) %>% 
  summarise(value = mean(value),n=length(unique(Species)))

SDM_evaluations %>% left_join(occ.nb) %>% 
  ggplot(aes(y = BOYCE_mean, x = n.occ)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()


### test differences between origin x habitat ====
m.eval.boyce <- lmer(BOYCE_mean ~ Origin * Habitat + (1|model) + (1|Species),
                     data = SDM_evaluations %>% filter(!Habitat %in% "Soil / Arboreal"))
car::Anova(m.eval.boyce, type = 3)
emmeans::emmeans(m.eval.boyce, ~ Origin | Habitat)
coef(m.eval.boyce)

m.eval.auc <- lmer(AUC_mean ~ Origin * Habitat + (1|model) + (1|Species),
                   data = SDM_evaluations %>% left_join(sp_char))
car::Anova(m.eval.auc, type = 3)

m.eval.tss <- lmer(TSS_mean ~ Origin * Habitat + (1|model) + (1|Species),
                   data = SDM_evaluations %>% left_join(sp_char))
car::Anova(m.eval.tss, type = 3)



SDM_evaluations2 <- SDM_evaluations %>% 
  mutate(model = fct_recode(model, 
                            `Random forest` = "raf",
                            GAM = "gam",
                            `Boosted regression tree` = "gbm",
                            MaxEnt = "max",
                            `Neural network` = "net",
                            GLM = "glm"
  )
  ) %>% 
  rename(
    `Boyce index` = BOYCE_mean,
    `AUC` = AUC_mean,
    `TSS` = TSS_mean,
    Parameters = parameters
  ) %>% 
  dplyr::select(1,2,3,5,7,9)

SDM_evaluations2 %>% group_by(model) %>% 
  summarise_at(2:4, mean)

GGally::ggpairs(SDM_evaluations2 %>% select(3:5))

SDM_evaluations %>% group_by(Origin) %>% summarise_at(c(3,5,7), mean)
SDM_evaluations %>% group_by(Habitat) %>% summarise_at(c(3,5,7), mean)

## variable importance ====
# Var_importance %>% 
#   left_join(sp.char) %>% 
#   ggplot(aes(y = mean, x = Variable, fill = Origin)) +
#   facet_grid(. ~ Habitat, scales = "free") +
#   geom_boxplot() +
#   scale_fill_brewer(palette="Set1") +
#   scale_y_continuous("Variable importance", expand = c(0,0)) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         legend.position = c(.9,.8),
#         legend.background = element_rect(color = "grey"),
#         strip.background = element_rect(fill = "black"),
#         strip.text = element_text(face = "bold", color = "white")
#   )


# Plot model predictions ====
martinique <- ne_countries(scale = 10, type = "map_units", country = 'france', returnclass = "sf") %>%
  filter(name_fr == "Martinique") %>% .$geometry

## by species ====
for(i in unique(occ$Species)){

  ### suitability maps ====
  p1 <- ggplot() +
    geom_sf(data = martinique, fill = NA) +
    geom_point(data = occ %>% filter(Species == i) %>% mutate(what = "Occurrences"),
               aes(x = X, y = Y), size = .5) +
    facet_wrap(~what)
  
  pred_temp <- rast(paste0(".././Results/ENM_predictions_", i, ".tif"))[[1:5]]
  names(pred_temp) <- c(
    "Ensemble\n(weighted mean suitability)",
    "Committee average\n(minimum training threshold)",
    "Committee average\n(10th percentile training threshold)",
    "Ensemble binary\n(minimum training threshold)",
    "Ensemble binary\n(10th percentile training threshold)"
  )
  p2 <- ggplot() +
    geom_spatraster(data = pred_temp[[1]], na.rm = T, maxcell = max.c) +
    scale_fill_viridis_c("Suitability", na.value = "transparent") +
    facet_wrap(~lyr)
  
  p3 <- ggplot() +
    geom_spatraster(data = as.factor(pred_temp[[2:3]]), na.rm = T, maxcell = max.c) +
    facet_wrap(~lyr) +
    scale_fill_viridis_d("No. models\nin agreement", na.value = "transparent", option = "B", na.translate = F,
                         limits = c("0","1","2","3","4","5","6"), drop = F)
  
  p4 <- ggplot() +
    geom_spatraster(data = as.factor(pred_temp[[4:5]]), na.rm = T, maxcell = max.c) +
    facet_wrap(~lyr) +
    scale_fill_manual(
      "", 
      values = c("lightgrey", "orange3"), 
      labels = c("Unsuitable", "Suitable", ""), 
      na.value = "transparent"
    )
  
  p_suit_tot <- ((p1 | p3 ) / (p2  | p4)) & 
    theme_void() & 
    theme(strip.text = element_text(margin = margin(t = 0, r = 0, b = 1, l = 0, unit = "pt"), size = 6),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.width = unit(10, "pt"))
  
  # ggsave(p_suit_tot, filename = paste0(".././Figures/Figure_suitability_",i,".pdf"))
  
  ### response to variables ====
  resp_temp <- Var_response %>% filter(Species == i) %>% 
    complete(Algorithm, Variable) %>% 
    group_by(Algorithm) %>% 
    mutate(y = range01(y)) %>% 
    group_by(Species, Variable, x) %>% 
    mutate(y_mean = mean(y), 
           y_min = quantile(y, p= .05, na.rm = T)[[1]],
           y_max = quantile(y, p= .95, na.rm = T)[[1]])
    
  
  p_var <- vector("list", length(unique(resp_temp$Variable)))
  for(j in 1:length(unique(resp_temp$Variable))){
    
    if(unique(resp_temp$Variable)[j] != "Soil"){
      p_var[[j]] <- resp_temp %>% filter(Variable == unique(resp_temp$Variable)[j]) %>% 
        ggplot() +
        geom_ribbon(aes(x = x, ymin = y_min, ymax = y_max), fill = "grey90") +
        geom_line(aes(x = x, y = y_mean), col = "grey20", linewidth = 1.2) +
        geom_line(aes(x = x, y = y, col = Algorithm), alpha = .3) +
        scale_color_manual(values = c("#B22222", "#5F9EA0", "#8470FF", "#63B8FF", "#FFA54F", "#66CD00"),
                           limits = c(
                             "Random Forest",
                             "GAM",
                             "Boosted regression trees",
                             "MaxEnt",
                             "Neural network",
                             "GLM"
                           )) +
        scale_y_continuous("Suitability", limits = c(0,1)) +
        scale_x_continuous(unique(resp_temp$Variable)[j]) +
        theme_classic() +
        theme(strip.placement = "outside", strip.background = element_blank(),
              axis.title.y = element_text(size = 7),
              axis.title.x = element_text(size = 7),
              axis.text = element_text(size = 6),
              legend.text = element_text(size = 7),
              legend.title = element_text(size = 7, face = "bold"),
              legend.key.height = unit(8, "pt"),
              legend.key.width = unit(5, 'pt')
        )
    } else {
      p_var[[j]] <- ggplot() +
        geom_bar(data = bind_rows(
          resp_temp %>% filter(Variable == unique(resp_temp$Variable)[j]) %>% 
            select(-c(y_mean,y_min,y_max)),
          resp_temp %>% filter(Variable == unique(resp_temp$Variable)[j]) %>%
            group_by(Variable, x, Species) %>% summarise(y = unique(y_mean)) %>% 
          mutate(Algorithm = "Mean")
          ) %>% mutate(x = as.factor(x)),
                 aes(x = x, y = y, fill = Algorithm),
                 stat = "identity", position = position_dodge(width = .9)) +
        geom_text(data = cbind.data.frame(label = tolower(legend_soil[1:8,2]),
                                          y = 0.01,
                                          x = 1:8),
                  aes(x=x, y = y , label = label), 
                  angle = 90, hjust = 0, vjust = 2, size = 2, color = "grey20") +
        scale_fill_manual(values = c("#B22222", "#5F9EA0", "#8470FF", "#63B8FF", "#FFA54F", "#66CD00"),
                           limits = c(
                             "Random Forest",
                             "GAM",
                             "Boosted regression trees",
                             "MaxEnt",
                             "Neural network",
                             "GLM"
                           ),
                          na.value = "grey20") +
        scale_y_continuous("Suitability", limits = c(0,1)) +
        scale_x_discrete(unique(resp_temp$Variable)[j]) +
        theme_classic() +
        theme(strip.placement = "outside", strip.background = element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = 7),
              axis.title.x = element_text(size = 7),
              axis.text = element_text(size = 6)
              )
    }
  }
  
  p_var_tot <- wrap_plots(p_var) +
    plot_layout(guides = "collect")
  
  ### evaluation table ====
  tbl_plot <- ggpubr::ggtexttable(
    SDM_evaluations2 %>% 
      filter(Species == i) %>%
      left_join(occ.nb) %>% 
      mutate(`Boyce index` = `Boyce index`) %>% 
      select_if(~ !any(is.na(.))) %>% 
      dplyr::select(-Species, - Parameters, -n.occ) %>% 
      column_to_rownames("model") %>% 
      mutate_all(.funs = function(x)format(round(x,2), nsmall = 2)),
    theme = ggpubr::ttheme(
      "light",
      base_size = 6,
      padding = unit(c(2, 2), "mm")
    )
  )
  
  
  tit <- title <- cowplot::ggdraw() + 
    cowplot::draw_label(i,
                        size = 10,
                        y= .8, x = 0, hjust = 0, fontface = "bold.italic") +
    cowplot::draw_label(paste("Number of occurrences:", nrow(occ %>% filter(Species == i))),
                        size = 8,
                        y = .4, x = 0, hjust = 0) +
    cowplot::draw_label(paste("Origin:", occ %>% filter(Species == i) %>% pull(Origin) %>% unique),
                        size = 8,
                        y = .2, x = 0, hjust = 0) +
    cowplot::draw_label(paste("Habitat:", occ %>% filter(Species == i) %>% pull(Habitat) %>% unique),
                        size = 8,
                        y = 0, x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 10, 5))
  
  
  p.final.sp <- cowplot::plot_grid(
    tit, 
    tbl_plot + theme(plot.margin = margin(0,0,10,-190)), 
    p_suit_tot,
    p_var_tot, 
    rel_heights = c(.2,.3,1,.8), 
    nrow = 4, ncol = 1,
    align = "v", axis = "r",
    labels = c("","A","B","C")
  )
  
  
  # Enregistrer
  ggsave(p.final.sp, filename = paste0(".././Figures/Appendix by species/",i,".png"), width = 7.5, height = 10.5)
}



# Overlap between introduced and native species ====
## all species ====
n.over <- nrow(unique(sum(pred_all_bin)))
ggplot() +
  geom_spatraster(data = as.factor(sum(pred_all_bin)), maxcell = max.c) +
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

## introduced and native, by habitat ====
sp_char2 <- bind_rows(
  sp_char %>% filter(Soil == 0) %>% mutate(Habitat2 = "Arboreal"),
  sp_char %>% filter(Soil == 1) %>% mutate(Habitat2 = "Soil"),
  sp_char %>% filter(Arboreal == 0) %>% mutate(Habitat2 = "Soil"),
  sp_char %>% filter(Arboreal == 1) %>% mutate(Habitat2 = "Arboreal")
) %>% unique

pred_bin_NatArbo <- subset(pred_all_bin, 
                           sp_char2 %>% filter(Habitat2 == "Arboreal", Origin == "Native") %>% 
                             pull(Species) %>% intersect(names(pred_all_bin))) %>% sum
pred_bin_NatSoil <- subset(pred_all_bin, 
                           sp_char2 %>% filter(Habitat2 == "Soil", Origin == "Native") %>% 
                             pull(Species) %>% intersect(names(pred_all_bin))) %>% sum
pred_bin_ExoArbo <- subset(pred_all_bin, 
                           sp_char2 %>% filter(Habitat2 == "Arboreal", Origin == "Exotic") %>% 
                             pull(Species) %>% intersect(names(pred_all_bin))) %>% sum
pred_bin_ExoSoil <- subset(pred_all_bin, 
                           sp_char2 %>% filter(Habitat2 == "Soil", Origin == "Exotic") %>% 
                             pull(Species) %>% intersect(names(pred_all_bin))) %>% sum


pred_bin_Nat <- c(pred_bin_NatArbo, pred_bin_NatSoil)
pred_bin_Exo <- c(pred_bin_ExoArbo, pred_bin_ExoSoil)
names(pred_bin_Nat) <- names(pred_bin_Exo) <- c("Arboreal", "Soil")

p_sr_nat <- ggplot() +
  geom_spatraster(data = as.factor(pred_bin_Nat), maxcell = max.c, na.rm = T) +
  facet_grid("Native species"~lyr, switch = "y") +
  theme_void() +
  theme(
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y = element_text(face = "bold", angle = 90, size = 14),
    legend.position = "none"
  )
p_sr_exo <- ggplot() +
  geom_spatraster(data = as.factor(pred_bin_Exo), maxcell = max.c, na.rm = T) +
  facet_grid("Exotic species"~lyr, switch = "y") +
  theme_void() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 90, size = 14)
  )


(p_sr_nat / p_sr_exo) &
  scale_fill_manual(
    "No. species",
    values = c("grey85", colorRampPalette(c("khaki1", "red4"))(6)),
    labels = c(0:6,""),
    na.value = "transparent"
  ) & plot_layout(guides = 'collect', ncol = 1)


ggsave(filename = ".././Figures/Fig1.pdf", width = 7, height = 7)

## pairwise introduced vs. native ====
sp_char3 <- sp_char2 %>% filter(Species %in% names(pred_all_bin))
overlap.all <- c()
for(i in sp_char3 %>% filter(Origin == "Native") %>% pull("Species")){
  focal.sp_rast <- subset(pred_all_bin, i)
  corr.sps <- sp_char3 %>% 
    filter(Origin == "Exotic", Habitat2 == sp_char3[sp_char3$Species == i, "Habitat2"][[1]]) %>% 
    pull(Species)
  
  corr.sps_rast_nat <- subset(pred_all_bin,i) %>% as.numeric()
  
  corr.sps_rast_over <- c()
  for(j in corr.sps){
    corr.sps_rast_exo <- subset(pred_all_bin, j) %>% as.numeric()
    corr.sps_rast_exo <- subst(corr.sps_rast_exo, 1, 2)
    
    overlap.tmp <- sum(corr.sps_rast_nat, corr.sps_rast_exo)
    levels(overlap.tmp) <- data.frame(id=0:3, cover=c("None", "Native only", "Exotic only", "Overlap"))
    
    corr.sps_rast_over <- c(corr.sps_rast_over, overlap.tmp)
  }
  corr.sps_rast_over <- rast(corr.sps_rast_over)
  names(corr.sps_rast_over) <- paste(i, "-", corr.sps)
  overlap.all <- c(overlap.all, corr.sps_rast_over)
}
names(overlap.all) <- sp_char3 %>% filter(Origin == "Native") %>% pull("Species")

## plot a matrix of overlap area ====
### Arboreal ====
plot.overlap_Arb <- lapply( 
  sp_char3 %>% filter(Origin == "Native", Habitat2 == "Arboreal") %>% pull(Species),
  function(x){
    ggplot() +
      geom_spatraster(data = setNames(overlap.all[[x]], 
                                      unlist(lapply(strsplit(names(overlap.all[[x]]), "-"), function(x)x[[2]]))), 
                      maxcell = max.c, na.rm = T) +
      facet_grid(get(substitute(x))~lyr, switch = "y") +
      scale_fill_manual(
        "",
        values = c("lightgrey", "green4", "yellow3", "red3"),
        na.value = "transparent", na.translate = F
      ) +
      theme_void() +
      theme(
        strip.text.x = element_text(face = "bold.italic", margin = margin(0,1,1,1), size = 7),
        strip.text.y = element_text(face = "bold.italic", angle = 90, margin = margin(0,1,1,1))
      )
  }
)

for(j in 2:length(plot.overlap_Arb)){
  plot.overlap_Arb[[j]] <- plot.overlap_Arb[[j]] + theme(strip.text.x = element_blank())
}

plot.overlap_Arb <- wrap_plots(plot.overlap_Arb) + 
  plot_layout(nrow = length(plot.overlap_Arb), guides = "collect")


### Soil ====
plot.overlap_Soil <- lapply( 
  sp_char3 %>% filter(Origin == "Native", Habitat2 == "Soil") %>% pull(Species),
  function(x){
    ggplot() +
      geom_spatraster(data = setNames(overlap.all[[x]], 
                                      unlist(lapply(strsplit(names(overlap.all[[x]]), "-"), function(x)x[[2]]))), 
                      maxcell = max.c, na.rm = T) +
      facet_grid(get(substitute(x))~lyr, switch = "y") +
      scale_fill_manual(
        "",
        values = c("lightgrey", "green4", "yellow3", "red3"),
        na.value = "transparent", na.translate = F
      ) +
      theme_void() +
      theme(
        strip.text.x = element_text(face = "bold.italic", margin = margin(0,2,2,1), size = 7),
        strip.text.y = element_text(face = "bold.italic", angle = 90, margin = margin(0,1,1,1))
      )
  }
)

for(j in 2:length(plot.overlap_Soil)){
  plot.overlap_Soil[[j]] <- plot.overlap_Soil[[j]] + theme(strip.text.x = element_blank())
}

plot.overlap_Soil <- wrap_plots(plot.overlap_Soil) + 
  plot_layout(nrow = length(plot.overlap_Soil), guides = "collect")



### merge pairwise plots ====
cowplot::plot_grid(plot.overlap_Arb, plot.overlap_Soil, ncol = 1, 
                   rel_heights = c(3,5), rel_width = c(5,6), align = "v", axis = "r") 
ggsave(filename = ".././Figures/Fig_pairwise_overlaps.pdf", width = 10, height = 15)



### overlaps by native species ====
overlap.byNative <- lapply(overlap.all, function(x){
  overlap.tmp <- x %>% as.numeric()
  overlap.tmp <- subst(overlap.tmp, 0:3, c(NA,0,NA,1))
  overlap.tmp <- overlap.tmp %>% sum %>% as.factor()
  return(overlap.tmp)
})

ggplot() +
  geom_spatraster(data = rast(overlap.byNative), maxcell = max.c, na.rm = T) +
  facet_wrap(~lyr, ncol = 4) +
  geom_sf(data = martinique, fill = NA, linewidth = .1) +
  scale_fill_manual(
    "No. overlapping\nexotic species",
    values = c("lightblue4", colorRampPalette(c("khaki1", "red4"))(5)),
    labels = c(0:5,""),
    na.value = "transparent"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold.italic", margin = margin(0,2,2,1), size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.width = unit(10, 'pt')
  ) 

ggsave(filename = ".././Figures/Fig_overlaps_map.pdf", width = 6, height = 3.5)

sum_overlaps <- lapply(
  overlap.all, 
  function(s){
    lapply(names(s), function(n) {
      res <- freq(s[[n]])
      res$layer <- n
      nam <- strsplit(n, " - ")[[1]][1]
      res$Habitat <- sp_char[sp_char$Species == nam, "Habitat"][[1]]
      return(res)
    })
  }
) %>% bind_rows()

sum_overlaps %>% group_by(layer) %>% 
  group_by(Habitat) %>%
  summarise(n = length(value[value == "Overlap"]))

sum_overlaps %>% filter(value == "Overlap") %>% arrange(count)

### synthesis ====
overlap_numbers <- rast(overlap.byNative) %>% as_tibble() %>%
  apply(.,2,table) %>% bind_rows %>% mutate(Species = names(overlap.byNative)) %>% 
  pivot_longer(-Species, names_to = "No. overlapping species", values_to = "No. cells")


overlap_numbers %>% 
  left_join(sp_char2) %>% 
  mutate(`No. cells` = ifelse(is.na(`No. cells`), 0, `No. cells`)) %>% 
  group_by(Species) %>% 
  mutate(tot_cells = sum(`No. cells`)) %>% 
  mutate(Prop_cells = `No. cells`/tot_cells,
         overlapped_range = sum(Prop_cells[`No. overlapping species` != 0])) %>% 
  ungroup %>% 
  mutate(Species = factor(Species)) %>% 
  mutate(Species = fct_reorder(Species, overlapped_range)) %>% 
  mutate(`No. overlapping species` = forcats::fct_rev(`No. overlapping species`)) %>%
  ggplot(aes(x = Species, fill = `No. overlapping species`, y = Prop_cells)) +
  facet_grid(Habitat2 ~ ., scale = "free", space = "free") +
  geom_bar(stat = 'identity') +
  scale_x_discrete("") +
  scale_y_continuous("Proportion of suitable range", labels = scales::percent) +
  scale_fill_manual(
    "No. overlapping\nexotic species",
    values = c(colorRampPalette(c("red4", "khaki1"))(5), "lightblue4")
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic'),
        strip.text = element_text(face = "bold"))

ggsave(filename = ".././Figures/Fig_overlaps_synthesis.pdf", width = 6, height = 3.5)

## synthetic invasion risk ====
### load PAs ====
PNR <- read_sf("../Data/PAs/mtq_pnr2013.shp")
RB <- read_sf("../Data/PAs/N_ENP_RB_S_972.shp")
# RNN <- read_sf("../Data/PAs/mtq_rn2010.shp")
# RNR <- read_sf("../Data/PAs/N_ENP_RNR_S_972.shp")

# overlap.byNative_all <- lapply(overlap.byNative, function(x){x[x>1]<-1;return(x)}) %>% rast %>% as.numeric
overlap.byNative_all <- overlap.byNative %>% rast %>% as.numeric
overlap.byNative_all <- sum(overlap.byNative_all, na.rm = T)

background <- env_all[[1]]; background[!is.na(background)] <- 0
# overlap.byNative_all <- sum(c(background, overlap.byNative_all), na.rm = T)

background <- as.polygons(background > -Inf)


p.synthesis_map <- ggplot() +
  geom_spatraster(data = overlap.byNative_all, na.rm = T, maxcell = max.c) +
  geom_sf(data = PNR, fill = "#f28d271A", color = NA) +
  geom_sf(data = background, fill = NA, color = "black", linewidth = .5) +
  geom_sf(data = RB, fill = NA, color = "blue", linewidth = .5) +
  scale_fill_viridis_c(
    "No. of overlaps\nbetween native\nand exotic species",
    option = "magma",
    na.value = 'transparent'
  ) +
  theme_void() +
  theme(legend.position = "none")


dat.hist <- as.data.frame(overlap.byNative_all) %>% table %>% as_tibble
p.synthesis_hist <- ggplot(data = dat.hist, aes(x = as.numeric(sum), y = n, fill = as.numeric(sum))) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = n), color = "black", vjust = 1, nudge_y = 1100, size = 2) +
  scale_fill_viridis_c(
    "No. of overlaps\nbetween native\nand exotic species",
    option = "magma",
    na.value = 'transparent'
  ) +
  scale_y_continuous("No. grid cells", expand = c(0, NA)) +
  scale_x_continuous("No. of overlaps between native and exotic species",
                     breaks = 0:16, 
                     expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")

cowplot::plot_grid(p.synthesis_map, p.synthesis_hist)

ggsave(filename = ".././Figures/Fig_sum_overlaps_raw.pdf", width = 10, height = 4.5)


# other species====
occ.rare <- list.files(".././Data", pattern = "OK", full.names = T) %>%
  lapply(.,read_excel) %>% bind_rows() %>% dplyr::select(nom_espece, coord_X, coord_Y) %>%
  rename(Species = 1, X = coord_Y, Y = coord_X) %>%
  mutate(Species = stringr::word(Species, 1, 2),
         Species = ifelse(Species == "Dichogaster sp11", "Dichogaster sp11A", Species),
         Species = ifelse(Species == "Glossodrilus sp. nov. 3", "Glossodrilus spp", Species),
         Species = ifelse(Species == "Glossodrilus sp.", "Glossodrilus sp1", Species)) %>%
  distinct %>% filter(X > -61.5)%>% 
  left_join(sp_char) %>% 
  group_by(Species) %>%
  mutate(n = n()) %>% filter(n < 5) %>% dplyr::select(-n)


ggplot() +
  geom_spatraster(data = as.factor(pred_bin_ExoSoil), maxcell = max.c) +
  geom_point(data = occ.rare, aes(x = X, y = Y, color = Species)) +
  ggrepel::geom_label_repel(data = occ.rare, aes(x = X, y = Y, label = Species, color = Species), 
                            max.overlaps = 10000,
                            label.size=NA, 
                            fill = NA,
                            fontface ="bold.italic") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_manual(
    "No. exotic\nsoil species",
    values = c("grey85", colorRampPalette(c("khaki1", "red4"))(6)),
    labels = c(0:6,""),
    na.value = "transparent"
  ) +
  theme_void()

ggsave(filename = ".././Figures/Fig_rareSp_raw.pdf", width = 8, height = 6)

occ.rare %>% ungroup %>% mutate(b.exo = extract(pred_bin_ExoSoil, occ.rare[,c("X","Y")])[,2]) %>% 
  group_by(Species, b.exo) %>% 
  summarise(n.occ = n())

