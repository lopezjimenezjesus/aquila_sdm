##################################################
## Project: Aquila  adalberti monograph
## Script purpose: Model future favourability
## Date: `sys.Date()`
## Author: Jesús Jiménez López
##################################################

##################################################
## Section: SET UP ----
##################################################

library(fuzzySim)
library(modEvA)
library(blorr)
library(foreign)
library(here)
library(sf)
library(readxl)
library(tidyverse)
library(corrplot)
library(GGally)
library(geodata)
library(mapSpain)
library(exactextractr)
library(caret)
library(gridExtra)
library(kableExtra)

here::i_am('Aquila_adalberti_SDM.Rproj')


# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

## Functions ----
# ################################################## 

source("02_SCRIPTS/functions.R")

##################################################
## Section: configure
##################################################

library(yaml)

config = yaml.load_file("config.yml")

model_name <- config$configuration$folder_name

paths <- create_folder(folder_name = config$configuration$folder_name) # create folder for saving output


##################################################
## Section: LOAD DATA  ----
##################################################


## Section: load sp and historical bioclimate vars
##################################################

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

utm1 <- st_read('01_DATA/INPUT/CUTM1_aquada/CUTM1x1Ext_AQUADA.shp') %>% dplyr::select("UTMCODE1X1")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariablesWldClimPresente.xlsx'))
  
y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x %>% dplyr::select("CUADRICULA", starts_with("Bio")), y) %>% 
  dplyr::select(-"CUADRICULA") %>% dplyr::rename("sp"="AQUADA")

aa <- data.frame(CUADRICULA=x$CUADRICULA, sp=xy$sp)

aa.sf <- inner_join(utm10, aa, by="CUADRICULA")


ext.raster <- rast(x = '01_DATA/INPUT/extremadura_10.tiff')  %>% 
  scale(center=TRUE, scale=TRUE) # present

names(ext.raster)  <- c(paste0("Bio0", 1:9), paste0("Bio", 10:19))


## Section: load future bioclimate variables
##################################################

path_to_tif <- "01_DATA/INPUT/bioclim_future/ext10KM/"

tiff_files_2021_2040 <- load_list_of_tif(path_to_tif, "2021-2040")
tiff_files_2041_2060 <- load_list_of_tif(path_to_tif, "2041-2060")
tiff_files_2061_2080 <- load_list_of_tif(path_to_tif, "2061-2080")
tiff_files_2081_2100 <- load_list_of_tif(path_to_tif, "2081-2100")

resample_and_scale_tif <- function(tif_file_path, from_raster) {
  r <- terra::rast(tif_file_path) # load raster
  names(r)  <- c(paste0("Bio0", 1:9), paste0("Bio", 10:19))
  # r <- terra::resample(r, from_raster)  %>% 
  #   scale(center=TRUE, scale=TRUE) # empty cells!
  r <-r  %>% 
    scale(center=TRUE, scale=TRUE)
  return(r)
}


rename_list_of_tif <- function(tiff_list) {
  names(tiff_list) <- c("126", "245", "370", "585")
  return(tiff_list)
}

tiff_files_2021_2040_scaled <- rename_list_of_tif(lapply(tiff_files_2021_2040, resample_and_scale_tif, ext.raster))
tiff_files_2041_2060_scaled <- rename_list_of_tif(lapply(tiff_files_2041_2060, resample_and_scale_tif, ext.raster))
tiff_files_2061_2080_scaled <- rename_list_of_tif(lapply(tiff_files_2061_2080, resample_and_scale_tif, ext.raster))
tiff_files_2081_2100_scaled <- rename_list_of_tif(lapply(tiff_files_2081_2100, resample_and_scale_tif, ext.raster))

## Get predictors from presence / absence data
##################################################

presence_predictors <- extract(ext.raster, aa.sf %>% st_transform(st_crs(ext.raster)) %>% st_centroid() %>% st_coordinates())  %>% scale()


## Add slope and alt
##################################################


utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

x_2 <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y_2 <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy_2 <- right_join(x_2,y_2, by="CUADRICULA")


paste(names(xy_2)[1:dim(xy_2)[2]], collapse = " + ")

metadata <- readr::read_delim(file = "01_DATA/OUTPUT/metadata_predictors.csv")

paisaje <- read_csv(file = "01_DATA/OUTPUT/paisaje_predictors.csv")

xy_2 <- xy_2 %>% dplyr::select(all_of(c(paisaje$predictors_names, "AQUADA")))

presence_predictors <- cbind(presence_predictors, alt=xy_2$ALT_MED, alt2=(xy_2$ALT_MED)*(xy_2$ALT_MED), slope=xy_2$Slope_mean)


xy.scaled <- cbind(presence_predictors, sp=aa.sf$sp) %>% as.data.frame()

##################################################
## Section: FUNCTIONS  ----
##################################################

## Functions
##################################################

historical_prediction <- function(bioclim_predictor, model) {
  pred <-terra::predict(object = bioclim_predictor  ,model=model, type="response",  na.rm=TRUE)
  return(pred)
}

historical_favourability <- function(historical_prediction, presence_absence_data, projection=TRUE, crs="4326") {
  f <- Fav(obs = presence_absence_data, pred = historical_prediction)
  if(projection) {
    f <- terra::project(f, crs)
  }
  
  return(f)
}

extract_fav_from_grid <- function(favourability_raster, umt_grid,  plot_points=FALSE) {
  
  fav.rast <-  exact_extract(favourability_raster, umt_grid, 'mean', append_cols  =TRUE)
  
  utm10.fav <- inner_join(umt_grid, fav.rast, by="CUADRICULA")
  
  utm10.fav$category <- cut(utm10.fav$mean, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  
  ggplot(data=utm10.fav) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) + {
      if(plot_points) geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(sp==1), size=0.5)} +
    theme_void() + theme(legend.position="none")
}


extract_fav_from_grid_2 <- function(favourability_raster, umt_grid, plot_points=FALSE) {
  
  fav.rast <-  exact_extract(favourability_raster, umt_grid, 'mean', append_cols  =TRUE)
  
  utm10.fav <- inner_join(umt_grid, fav.rast, by="CUADRICULA")
  
  utm10.fav$category <- cut(utm10.fav$mean, breaks=c(-Inf, 0.5, Inf), labels=c("low","high"))
  
  ggplot(data=utm10.fav) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = c("low" = "red", "high" = "green")) + {
      if(plot_points) geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(sp==1), size=0.5)} +
    theme_void() + theme(legend.position="none")
  
}

## More Functions
##################################################

future_prediction <- function(bioclim_predictor, model) {
  pre <- predict(object = bioclim_predictor, model=model,
                 type="response")
  return(pre)
}

future_favourability<- function(future_prediction, presence_absence ) {
  f <- Fav(obs = presence_absence, pred = future_prediction)
  return(f)
}


favourability_to_vector <- function(favourability_raster, breaks= c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high")) {
  fav.rast <-  exact_extract(favourability_raster, utm10, 'mean', append_cols  =TRUE)
  utm10.fav <- inner_join(utm10, fav.rast, by="CUADRICULA")
  utm10.fav$category <- cut(utm10.fav$mean, breaks = breaks, labels=labels)
  #utm10.fav <- utm10.fav %>% fill( category, .direction = 'updown' ) # na values due to grid resolution
  return(utm10.fav)
}

plot_favourability <- function(favourability_vector, name="", values =c("low" = "red", "middle" = "yellow", "high" = "green")) {
  p <- ggplot(data=favourability_vector) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = values) +
    #geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(sp==1)) +
    labs(title=name) +
    theme_void() + 
    theme(legend.position ="bottom") 
  return(p)
}

##################################################
## Section: MODELS  ----
##################################################

## Model 1: All variables
##################################################

null.model <- glm(formula = sp ~ 1, family = binomial, data=xy.scaled)

bmod_aa_biohist <- step(null.model, direction = "forward", 
                        keep = function(model, aic) list(model = model, aic = aic),
                        scope = (~Bio01 + Bio02 + Bio03 + Bio04 + Bio05 + Bio06 +
                                   Bio07 + Bio08 + Bio09 + Bio10 + 
                                   Bio11 + Bio12 + Bio13 + Bio14 +
                                   Bio15 + Bio16 + Bio17 + Bio18 + Bio19))

bmod_aa_biohist.trimmmed <- modelTrim(bmod_aa_biohist)

# 
# foo <-  glm(formula = sp ~ Bio09 + Bio03 + Bio14 + Bio15 + Bio08 + Bio13 , family = binomial, data=xy.scaled)
# 
# 

## Model 2: by Feature selection
##################################################

# make_model <- function(nm) lm(xy.scaled[c("sp", nm)])
# fits <- Map(make_model, xy.scaled)
# 
# glance_tidy <- function(x) c(unlist(glance(x)), unlist(tidy(x)[, -1]))
# out <- sapply(fits, glance_tidy)

# FDR

var_cols = dim(xy.scaled)[2] - 1

xy.scaled.selection <- FDR(data=xy.scaled,sp.cols = dim(xy.scaled)[2], var.cols = 1:var_cols)

bioclim_predictors <- sort(rownames(xy.scaled.selection$select))

xy.scaled.selection <- (xy.scaled  %>% dplyr::select(all_of(bioclim_predictors),  "sp"))

# Correlation

cor_m <- cor(xy.scaled.selection %>% dplyr::select(!all_of("sp")))
index <- sort(findCorrelation(cor_m, cutoff = 0.9))

xy.scaled.selection <- xy.scaled.selection[,  -index] # remove correlated variables

# xy.scaled.selection <- cbind(xy.scaled.selection, sp=xy$sp)


variables_to_model <- paste(colnames(xy.scaled.selection)[-length(colnames(xy.scaled.selection))], collapse = " + ")


null.model <- glm(formula = sp ~ 1, family = binomial, data=xy.scaled.selection)


fmod_aa_biohist <- step(null.model, direction = "forward", 
                        keep = function(model, aic) list(model = model, aic = aic),
                        scope =paste0("~", variables_to_model))


summary(fmod_aa_biohist)

model_n_steps <- length(fmod_aa_biohist$keep) / 2

saveRDS(object = model_n_steps, file = file.path(paths$model_path, "future_n_steps_model.rds"))

fmod_aa_biohist.trimmed <- modelTrim(fmod_aa_biohist)

# 
# nest_model.glm  <- fuzzySim::stepwise(xy.scaled.selection, length(xy.scaled.selection), 1:(length(xy.scaled.selection) -1), family = binomial(link="logit"), 
#                                       simplif=FALSE, direction="forward", trace=2,preds=TRUE, Favourability=FALSE,Wald=TRUE)
# 
# 

## Plot both models
##################################################

# plot(historical_prediction(ext.raster, bmod_aa_biohist))
# plot(historical_prediction(ext.raster, fmod_aa_biohist))
# plot(historical_prediction(ext.raster, fmod_aa_biohist.trimmed))

##################################################
## Section: HISTORICAL - FUTURE FAVOURABILITY  ----
##################################################

##################################################
## MODEL TO EXPLORE ----
##################################################

model_to_explore <- fmod_aa_biohist.trimmed

write_rds(model_to_explore, file = file.path(paths$model_path, "future_model_selected.rds"))

## Historical ----
##################################################

epsg_code <- paste0("epsg:",st_crs(utm10)$epsg)

plot(historical_favourability(historical_prediction(ext.raster, model_to_explore), 
                              xy.scaled$sp, crs=epsg_code))

historical_cat_3 <- extract_fav_from_grid(historical_favourability(historical_prediction(ext.raster, model_to_explore), 
                                               xy.scaled$sp, crs=epsg_code),utm10)

historical_cat_2 <- extract_fav_from_grid_2(historical_favourability(historical_prediction(ext.raster, model_to_explore), 
                                                 xy.scaled$sp, crs=epsg_code),utm10)

saveRDS(historical_cat_3, file=file.path(paths$model_path, "historical_cat_3.rds"))
saveRDS(historical_cat_2, file=file.path(paths$model_path, "historical_cat_2.rds"))

## Future  ----
##################################################

## 2021 - 2040
##################################################

# 3 categories

prediction_2021_2040 <-  lapply(tiff_files_2021_2040_scaled, future_prediction, model_to_explore)

favourability_2021_2040 <-  lapply(prediction_2021_2040, future_favourability, xy.scaled$sp)

favourability_vector_2021_2040 <- lapply(favourability_2021_2040, favourability_to_vector)

favourability_2021_2040_plot <- lapply(favourability_vector_2021_2040, plot_favourability)

favourability_2021_2040_plot_ <-lapply(
  names(favourability_2021_2040_plot),
  function(i) favourability_2021_2040_plot[[i]] + 
    labs(title=paste0(i, " (2021-2040)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)


do.call("grid.arrange", c(favourability_2021_2040_plot_, ncol = 2))


# 2 categories

favourability_vector_2021_2040 <- lapply(favourability_2021_2040, 
                                         favourability_to_vector, 
                                         breaks= c(-Inf, 0.5, Inf), labels=c("low","high"))

favourability_2021_2040_plot <- lapply(favourability_vector_2021_2040, 
                                       plot_favourability, 
                                       values =c("low" = "red", "high" = "green"))

favourability_2021_2040_plot_2 <-lapply(
  names(favourability_2021_2040_plot),
  function(i) favourability_2021_2040_plot[[i]] + 
    labs(title=paste0(i, " (2021-2040)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)


do.call("grid.arrange", c(favourability_2021_2040_plot_2, ncol = 2))

## 2040 - 2061
##################################################

# 3 categories

prediction_2041_2060 <-  lapply(tiff_files_2041_2060_scaled, future_prediction, model_to_explore)

favourability_2041_2060 <-  lapply(prediction_2041_2060, future_favourability, xy.scaled$sp)

favourability_vector_2041_2060 <- lapply(favourability_2021_2040, favourability_to_vector)

favourability_2041_2060_plot <- lapply(favourability_vector_2041_2060, plot_favourability)

favourability_2041_2060_plot_ <-lapply(
  names(favourability_2021_2040_plot),
  function(i) favourability_2041_2060_plot[[i]] + 
    labs(title=paste0(i, " (2041-2060)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)


do.call("grid.arrange", c(favourability_2041_2060_plot_, ncol = 2))

# 2 categories

favourability_vector_2041_2060 <- lapply(favourability_2041_2060, 
                                         favourability_to_vector, 
                                         breaks= c(-Inf, 0.5, Inf), labels=c("low","high"))

favourability_2041_2060_plot <- lapply(favourability_vector_2041_2060, 
                                       plot_favourability, 
                                       values =c("low" = "red", "high" = "green"))

favourability_2041_2060_plot_2 <-lapply(
  names(favourability_2041_2060_plot),
  function(i) favourability_2041_2060_plot[[i]] + 
    labs(title=paste0(i, " (2041-2060)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)


do.call("grid.arrange", c(favourability_2041_2060_plot_2, ncol = 2))


## 2061 - 2080
##################################################

# 3 categories

prediction_2061_2080 <-  lapply(tiff_files_2061_2080_scaled, future_prediction, model_to_explore)

favourability_2061_2080 <-  lapply(prediction_2061_2080, future_favourability, xy.scaled$sp)

favourability_vector_2061_2080 <- lapply(favourability_2061_2080, favourability_to_vector)

favourability_2061_2080_plot <- lapply(favourability_vector_2061_2080, plot_favourability)

favourability_2061_2080_plot_ <-lapply(
  names(favourability_2061_2080_plot),
  function(i) favourability_2061_2080_plot[[i]] +
    labs(title=paste0(i, " (2061-2080)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)

do.call("grid.arrange", c(favourability_2061_2080_plot_, ncol = 2))

# 2 categories

favourability_vector_2061_2080 <- lapply(favourability_2061_2080, 
                                         favourability_to_vector, 
                                         breaks= c(-Inf, 0.5, Inf), labels=c("low","high"))

favourability_2061_2080_plot <- lapply(favourability_vector_2061_2080, 
                                       plot_favourability, 
                                       values =c("low" = "red", "high" = "green"))

favourability_2061_2080_plot_2 <-lapply(
  names(favourability_2061_2080_plot),
  function(i) favourability_2061_2080_plot[[i]] + 
    labs(title=paste0(i, " (2061-2080)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)


do.call("grid.arrange", c(favourability_2061_2080_plot_2, ncol = 2))



## 2081 - 2100
##################################################

# 3 categories

prediction_2081_2100 <-  lapply(tiff_files_2081_2100_scaled, future_prediction, model_to_explore)

favourability_2081_2100 <-  lapply(prediction_2081_2100, future_favourability, xy.scaled$sp)

favourability_vector_2081_2100 <- lapply(favourability_2081_2100, favourability_to_vector)

favourability_2081_2100_plot <- lapply(favourability_vector_2081_2100, plot_favourability)

favourability_2081_2100_plot_ <-lapply(
  names(favourability_2081_2100_plot),
  function(i) favourability_2081_2100_plot[[i]] +
    labs(title=paste0(i, " (2081-2100)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
)

do.call("grid.arrange", c(favourability_2081_2100_plot_, ncol = 2))

# 2 categories

favourability_vector_2081_2100 <- lapply(favourability_2081_2100, 
                                         favourability_to_vector, 
                                         breaks= c(-Inf, 0.5, Inf), labels=c("low","high"))

favourability_2081_2100_plot <- lapply(favourability_vector_2081_2100, 
                                       plot_favourability, 
                                       values =c("low" = "red", "high" = "green"))

favourability_2081_2100_plot_2 <-lapply(
  names(favourability_2081_2100_plot),
  function(i) favourability_2081_2100_plot[[i]] + 
    labs(title=paste0(i, " (2081-2100)")) + 
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") 
)


do.call("grid.arrange", c(favourability_2081_2100_plot_2, ncol = 2))



## Tendencias ---
##################################################

## Layout all
##################################################

x1 <- c(rbind(favourability_2021_2040_plot_,favourability_2021_2040_plot_2))
x2 <- c(rbind(favourability_2041_2060_plot_,favourability_2041_2060_plot_2))
x3 <- c(rbind(favourability_2061_2080_plot_,favourability_2061_2080_plot_2))
x4 <- c(rbind(favourability_2081_2100_plot_,favourability_2081_2100_plot_2))

x <-  c(rbind(x1,x2, x3, x4))

# Get number of cells in each category (low, middle, high) for all models

# Patrón de expresión regular
patron <- "^(\\d+) \\((\\d+-\\d+)\\)$"

cells_in_each_category <- sapply(x, function(x){
  # Aplica la expresión regular y captura los grupos
  resultados <- str_match(x$label$title, patron)
  
  df <- data.frame(rbind(table(x$data$category)), title= x$label$title, escenario=resultados[2], periodo=resultados[3],
                   type=if_else(length(names(table(x$data$category))) == 2, "cat_2", "cat_3"))
  return(df)
})


cells_in_each_category <-bind_rows(cells_in_each_category, .id = "column_label")

cells_in_each_category <- pivot_longer(cells_in_each_category, c(low, middle, high), values_to = "Value", names_to = "category")


cells_in_each_category <- cells_in_each_category %>% group_by(periodo, escenario) %>%  mutate(percent=Value*100 / 516)


cells_historical_3 <- historical_cat_3$data %>% 
  group_by(category) %>% 
  summarise(Value=n()) %>% 
  st_drop_geometry() %>% 
  mutate(percent=Value*100 / 516) %>% 
  mutate(column_label="", title="historical", escenario="0", periodo="0", type="cat_3", variacion=NA) %>%
  select(column_label, title, escenario, periodo, type, category, Value, variacion)
  

cells_historical_2 <- historical_cat_2$data %>% 
  group_by(category) %>% 
  summarise(Value=n()) %>% 
  st_drop_geometry() %>% 
  mutate(percent=Value*100 / 516) %>%
  mutate(column_label="", title="historical", escenario="0", periodo="0", type="cat_2", variacion=NA) %>%
  select(column_label, title, escenario, periodo, type, category, Value, variacion)



high_f_cat2 <- cells_historical_2[cells_historical_2$category=="high",]$Value
high_f_cat3 <- cells_historical_3[cells_historical_3$category=="high",]$Value
# Variación = [Periodo final – Periodo inicial] x 100
# 
# Periodo Inicial


cells_in_each_category <- cells_in_each_category %>% mutate(variacion= if_else(type=="cat_2", ((Value- high_f_cat2) / high_f_cat2)*100, ((Value- high_f_cat3) / high_f_cat3)*100))


# Crea el gráfico de tendencias

ggplot(cells_in_each_category, aes(x = periodo, y = Value, color = escenario, group = interaction(escenario, category))) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, aes(group = interaction(type, category)), color = "black", linetype = "dashed") +
  facet_wrap(~type+category , scales = "free_y") +
  labs(title = "Tendencias por Escenario",
       x = "Período",
       y = "Valor") +
  theme_minimal()

# add historical data

cells_in_each_category <- rbind(cells_in_each_category, cells_historical_2) %>% rbind(cells_historical_3)

# 
# ggplot(data=cells_summary) + 
#   geom_line(aes(x = ))

write_rds(cells_in_each_category,  file.path(paths$model_path, "cells_in_each_category.rds"))

## Arrange x plots

xx <- do.call("grid.arrange", c(x , ncol = 4))

lapply(x, function(p) {
  file_name <-  paste0(p$labels$title, '.png')
  ggsave(file_name, p  + labs(title = ""), path = paths$figure_path,device = "png")
})


ggsave(filename=file.path(paths$figure_path, "future_prediction.pdf"), plot=xx)

ggsave(filename = file.path(paths$figure_path, "future_prediction.png"), 
       plot = xx, dpi = "retina", width = 210, height = 297, units = "mm")


## Plot all in a grid
##################################################

library(patchwork)
library(grid)

plot_grid <- function() {

p126_2cat <- x1[[2]] +labs(title="2021-2040") + theme(plot.title = element_text(size=10)) +
  ylab("126")  + 
  theme(axis.title.y = element_text(angle = 90, hjust = -0.325, size = 10)) + 
  x2[[2]] +labs(title="2041-2060")  + theme(plot.title = element_text(size=10)) +
  x3[[2]] +labs(title="2061-2080")  + theme(plot.title = element_text(size=10)) +
  x4[[2]] +labs(title="2081-2100")  + theme(plot.title = element_text(size=10)) +
  plot_layout(nrow = 1, byrow = FALSE)

p126_3cat <- x1[[1]] + labs(title=NULL) + x2[[1]] + labs(title=NULL) + x3[[1]] +
  labs(title=NULL) + x4[[1]] + labs(title=NULL)  + plot_layout(nrow = 1, byrow = FALSE)

p245_2cat <- x1[[4]] + labs(title=NULL) +    
  ylab("245")  + 
  theme(axis.title.y = element_text(angle = 90, hjust = -0.325, size = 10)) + 
  x2[[4]] + labs(title=NULL) + 
  x3[[4]] + labs(title=NULL) +   
  x4[[4]] + labs(title=NULL) + 
  plot_layout(nrow = 1, byrow = FALSE) 

p245_3cat <- x1[[3]] + labs(title=NULL) + x2[[3]] + labs(title=NULL)  + x3[[3]] + 
  labs(title=NULL) + x4[[3]] + labs(title=NULL)  +
  plot_layout(nrow = 1, byrow = FALSE)

p370_2cat <- x1[[6]]  + labs(title=NULL) + 
  ylab("370")  + 
  theme(axis.title.y = element_text(angle = 90, hjust = -0.325, size = 10)) + 
  x2[[6]]  + labs(title=NULL) +  
  x3[[6]]  + labs(title=NULL) + 
  x4[[6]]  + labs(title=NULL) + 
  plot_layout(nrow = 1, byrow = FALSE) 

p370_3cat <- x1[[5]] + labs(title=NULL)  + x2[[5]] + labs(title=NULL)  + x3[[5]] + labs(title=NULL)  +
  x4[[5]] + labs(title=NULL)  + plot_layout(nrow = 1, byrow = FALSE)


p585_2cat <- x1[[8]]  + labs(title=NULL) + 
  ylab("585")  + 
  theme(axis.title.y = element_text(angle = 90, hjust = -0.325, size = 10)) + 
  x2[[8]]  + labs(title=NULL) +  
  x3[[8]]  + labs(title=NULL) +   
  x4[[8]]  + labs(title=NULL) + 
  plot_layout(nrow = 1, byrow = FALSE) 

p585_3cat <- x1[[7]] + labs(title=NULL) + x2[[7]] + labs(title=NULL)+ x3[[7]] + labs(title=NULL)  +
  x4[[7]] + labs(title=NULL)  + plot_layout(nrow = 1, byrow = FALSE)


h_patch <- p126_2cat / p126_3cat + plot_layout(nrow = 2, byrow = FALSE) &  theme(plot.margin = unit(c(0,0,0,0), "cm"))

h_patch2 <- p245_2cat / p245_3cat + plot_layout(nrow = 2, byrow = FALSE) &  theme(plot.margin =unit(c(0,0,0,0), "cm"))

h_patch3 <- p370_2cat / p370_3cat + plot_layout(nrow = 2, byrow = FALSE)  &  theme(plot.margin =unit(c(0,0,0,0), "cm"))

h_patch4 <- p585_2cat / p585_3cat + plot_layout(nrow = 2, byrow = FALSE)  &  theme(plot.margin = unit(c(0,0,0,0), "cm"))



future_layout <- (wrap_elements(h_patch) +
  # labs(tag = "126") +
  theme(
    plot.tag = element_text(size = rel(0.5), angle = 90),
    plot.tag.position = "left"
  )) /  plot_spacer() / 
  
  (wrap_elements(h_patch2)  +
     # labs(tag = "245") +
     theme(
       plot.tag = element_text(size = rel(0.5), angle = 90),
       plot.tag.position = "left"
     )) /  plot_spacer() / 
  
  (wrap_elements(h_patch3)  +
     # labs(tag = "245") +
     theme(
       plot.tag = element_text(size = rel(0.5), angle = 90),
       plot.tag.position = "left"
     ))  /  plot_spacer() / 
  
  (wrap_elements(h_patch4)  +
     # labs(tag = "245") +
     theme(
       plot.tag = element_text(size = rel(0.5), angle = 90),
       plot.tag.position = "left"
     )) +  
  plot_layout(heights = c(8, 0, 8, 0, 8, 0, 8), guides = "collect")

  future_layout <-  p126_2cat / p126_3cat /  p245_2cat / p245_3cat / p370_2cat / p370_3cat / p585_2cat / p585_3cat


  return(future_layout)
}

future_layout <- plot_grid() +  plot_annotation(theme = theme(plot.margin = margin(r=-1, l=-1)))

ggsave(filename = file.path(paths$figure_path, "future_prediction_2.png"), 
       plot = future_layout,  dpi = "retina", width = 210, height = 260, units = "mm")


knitr::plot_crop( file.path(paths$figure_path, "future_prediction_2.png"))

##################################################
## Section: Metrics ----
##################################################

## model_to_explore

# hosmer BLR ----

blr <- blr_test_hosmer_lemeshow(model_to_explore)

blr


saveRDS(object = blr, file = file.path(paths$model_path, "future_blr.rds"))


# percent cells 

fav_3_cat <- extract_fav_from_grid(historical_favourability(historical_prediction(ext.raster, model_to_explore), 
                                               xy.scaled$sp, crs=epsg_code),utm10)

summary(fav_3_cat$data)

fav_2_cat <- extract_fav_from_grid_2(historical_favourability(historical_prediction(ext.raster, model_to_explore), 
                                                 xy.scaled$sp, crs=epsg_code),utm10)
summary(fav_2_cat$data)


# Plot 

grid.arrange(fav_2_cat + theme(legend.position = "none"), fav_3_cat + theme(legend.position = "none"), ncol=2)

# table of coefficients


build_coefficient_table <- function(model) {
  
  tidy_coef <- broom::tidy(model)
  
  odd.ratio <- exp(coef(model))
  
  wald_coef <- sapply(1:length(coef(model)), function(x){
    temp <- aod::wald.test(b=coef(model), Sigma=vcov(model), Terms = x)
    temp$result$chi2[1]
  })
  
  
  tidy_coef <- cbind(tidy_coef, odd.ratio) %>% cbind(wald=wald_coef)
  
  colnames(tidy_coef) <- c("Variables", "β", "ET", "statistic", "Sig.", "Exp(B)", "Wald")
  
  tidy_coef = tidy_coef[-1,] # remove intersect
  
  
  tidy_coef <- tidy_coef %>% dplyr::select(c("Variables", "β", "ET","Wald", "Sig.", "Exp(B)"))
  return(tidy_coef)
  
}

tidy_coef <- build_coefficient_table(model_to_explore)


write_rds(tidy_coef, file = file.path(paths$model_path, "tabla_coef_future.rds"))



# 
# knitr::kable(tidy_coef, "html") %>%
#   kable_styling("striped") %>% 
#   row_spec(0, background = "#F7CAAC") %>%
#   row_spec(seq(1,dim(tidy_coef)[1],2), background = "#FBE4D5")
# 
# 
