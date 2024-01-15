##################################################
## Project: 
## Script purpose:
## Date:
## Author:
##################################################

## Section: Set up
##################################################

library(fuzzySim)
library(modEvA)
library(mecofun)
library(foreign)
library(here)
library(sf)
library(readxl)
library(tidyverse)
library(corrplot)
library(blorr)
library(gridExtra)
library(broom)
library(flextable)
library(kableExtra)
library(mecofun)
library(landscapemetrics)
library(writexl)

here::i_am('Aquila_adalberti_SDM.Rproj')

## Functions ----
# ################################################## 

source("02_SCRIPTS/functions.R")

## Load data
##################################################

landscape <- terra::rast("01_DATA/INPUT/20240108_PaisajeHabitat/clc18ext.tif")
landscape_buffer <- st_read("01_DATA/INPUT/20240108_PaisajeHabitat/AQUADA_2019_2021_PS_Habitat_Random_BuffExt.shp")

landscape <- terra::project(landscape, landscape_buffer)

landscape_mask <- terra::mask(x = landscape, mask = landscape_buffer)

metadata_predictors <- readr::read_delim('01_DATA/OUTPUT/metadata_predictors.csv') %>%
  filter(!is.na(B_Habitat)) %>% select(predictors_names)

xy <- read_excel(path = here::here('01_DATA/INPUT/VARIABLES_Habitat_AQUADA.xlsx'), sheet = "medias") %>% 
  mutate(presence= if_else(grepl('^Habitat', ID), 1, 0)) %>% select(-ID) %>%
  # select(metadata_predictors) %>%
  select(-c(AltRan, Oeste180, Sur180, TabsMax1, TabsMax7, TabsMin1, TabsMin7,
            RainMax1, RainMax7, RainDay1, RainDay7, areaExt_m2
            
  ))


##################################################
## Section: configure
##################################################

library(yaml)

config = yaml.load_file("config.yml")

model_name <- config$configuration$folder_name

paths <- create_folder(folder_name = config$configuration$folder_name) # create folder for saving output


##################################################
## Section: landscape metrics
##################################################

check_landscape(landscape_mask)

terra::plot(landscape_mask)

# lsm_metrics <- calculate_lsm(landscape_mask, type="aggregation metric",
#                              level = c("landscape"),
#                              what = c("lsm_l_np", "lsm_l_pd", "lsm_l_lpi", "lsm_l_lsi", 
#                                       "lsm_l_area_mn", "lsm_l_frac_mn", "lsm_l_contag", "lsm_l_shdi", "lsm_l_shei"))


if(file.exists("01_DATA/OUTPUT/habitat_fragsat.csv")) {
  
  lsm_metrics_sample_long <- read_csv( "01_DATA/OUTPUT/habitat_fragsat.csv")
  
} else {
  
  lsm_metrics_sample <-  sample_lsm(landscape, y = landscape_buffer, level = c("landscape"),
                                    what = c("lsm_l_np", "lsm_l_pd", "lsm_l_lpi", "lsm_l_lsi", 
                                             "lsm_l_area_mn", "lsm_l_frac_mn", "lsm_l_contag", "lsm_l_shdi", "lsm_l_shei"))
  
  
  lsm_metrics_sample_long <- lsm_metrics_sample %>% 
    pivot_wider(names_from = metric, values_from = value)
  
  write_csv(lsm_metrics_sample_long, file = "01_DATA/OUTPUT/habitat_fragsat.csv")
  write_xlsx(lsm_metrics_sample_long,  path  = "01_DATA/OUTPUT/habitat_fragsat.xlsx")
  
}



# check raster

lsm_metrics_sample_r <-  sample_lsm(landscape, y = landscape_buffer, level = c("landscape"),
                                  what = c("lsm_l_np", "lsm_l_pd", "lsm_l_lpi", "lsm_l_lsi", 
                                           "lsm_l_area_mn", "lsm_l_frac_mn", "lsm_l_contag", "lsm_l_shdi", "lsm_l_shei"), return_raster=TRUE)

terra::plot(lsm_metrics_sample_r$raster_sample_plots[[1]])


lsm_metrics_sample_long <- lsm_metrics_sample_long %>% select(area_mn , contag, frac_mn, lpi, lsi, np, shdi, shei)


##################################################
## Section: join xy and landscape metrics
##################################################

xy <-  cbind(lsm_metrics_sample_long, xy) 


# ## Collinearity: VIF
# ##################################################
# 
# #https://quantifyinghealth.com/vif-threshold/
# 
# mult1 <- multicol(subset(xy, select = 1:(ncol(xy)-1)))
# mult1  # muestra la multicolinealidad
# 
# 
# VIF_threshold <- config$configuration$vif_threshold
# 
# xy.vif <- subset(mult1, VIF < VIF_threshold)
# 

## FDR - Correlation
##################################################

xy.fdr <- FDR(data = xy, sp.cols = length(xy), var.cols = 1:(length(xy)-1), q = 0.05)

xy.fdr$exclude

saveRDS(xy.fdr, file =  file.path(paths$model_path, "habitat_fdr.rds"))

variables_to_model <- paste(rownames(xy.fdr$select)[1:length(rownames(xy.fdr$select))], collapse = " + ")

## Fit model
##################################################

null.model <- glm(presence ~ 1, data=xy, family = binomial)

habitat_model.glm <- step(null.model, direction='forward',
                  keep =  function(model, aic) list(model = model, aic = aic),
                  scope = paste0("~", variables_to_model))

summary(habitat_model.glm)

saveRDS(object = habitat_model.glm, file =  file.path(paths$model_path, "habitat_model_steps.rds"))

habitat_model.glm <- modelTrim(method = "summary", habitat_model.glm) # using anova to remove non signifincant

fav <- data.frame(fav=Fav(habitat_model.glm))
colnames(fav) <- "fav"

saveRDS(object = habitat_model.glm, file =  file.path(paths$model_path, "habitat_model.rds"))

ggplot(data = fav, mapping = aes(fav))+
  geom_histogram(bins = 10, col="white") + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  # scale_x_continuous(labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4","0.4-0.5", "0.5-0-6", "0.6-0.7", "0.7-0.8","0.8-0.9", "0.9-1")) +
  theme_minimal()

fav$interval <- cut(fav$fav,10, dig.lab = 1)

fav$interval_2 <- cut(fav$fav,breaks=c(0,0.2,0.8,1))

xy$fav <- fav$fav
xy$interval <- as.factor(fav$interval)
xy$interval_2 <- as.factor(fav$interval_2)
xy.summarize <- xy %>% filter(presence==1) %>% group_by(presence, interval) %>% summarise(n=n())
xy.summarize_2 <- xy %>% filter(presence==1) %>% group_by(presence, interval_2) %>% summarise(n=n())

p <- ggplot(data = xy.summarize , mapping = aes(x = interval, y=n))+
  geom_bar( fill="#70AD47", stat = "identity") + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  scale_x_discrete(labels = c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4","0.4 -0.5 ", "0.5 - 0-6", "0.6 - 0.7", "0.7 - 0.8","0.8 - 0.9", "0.9 - 1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

p

ggsave(filename = file.path(paths$figure_path, 'habitat_barplot_10.png'), 
       plot = p, dpi = "retina")


q <- ggplot(data = xy.summarize_2 , mapping = aes(x = interval_2, y=n))+
  geom_col(aes(fill=interval_2)) + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  scale_x_discrete(labels = c("0 - 0.2", "0.2 - 0.8", "0.8 - 1")) +
  scale_fill_manual( values = c("red", "yellow",  "#70AD47")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none")

q

ggsave(filename = file.path(paths$figure_path, 'habitat_barplot_3.png'), 
       plot = q, dpi = "retina")


## Coefficients ----
##################################################


tidy_coef <- build_coef_table(habitat_model.glm)

write_rds(tidy_coef, file = file.path(paths$model_path, "habitat_tabla_coef.rds"))

## Metrics
##################################################

blr <- blr_test_hosmer_lemeshow(habitat_model.glm)

crosspred_glm <- mecofun::crossvalSDM(habitat_model.glm, traindat= xy,
                                      colname_pred=colnames(xy)[1:(dim(xy)[2]-1)], 
                                      colname_species = "presence", kfold= 10)

write_rds(blr, file = file.path(paths$model_path, "habitat_blr.rds"))


evalsdm <-  mecofun::evalSDM(xy$presence, crosspred_glm)

tm <- threshMeasures(model = habitat_model.glm, ylim = c(0, 1), thresh = 0.5)

c_m <-  as.data.frame(tm$ConfusionMatrix)

rownames(c_m) <- c("Favorable", "Desfavorable")
colnames(c_m) <- c("Presencia", "Ausencia")
c_m


model.metrics <- data.frame(AUC=evalsdm[c("AUC")], 
                            UPR=tm$ThreshMeasures[c("UPR"),],
                            OPR=tm$ThreshMeasures[c("OPR"),],
                            HyL=blr$pvalue)
model.metrics

saveRDS(object = model.metrics, file =  file.path(paths$model_path, "habitat_model_metrics.rds"))
