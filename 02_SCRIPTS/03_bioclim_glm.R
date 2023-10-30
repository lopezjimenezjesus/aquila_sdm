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
library("gridExtra")    

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


##################################################
## Section: load data
##################################################


## Section: load sp and historical bioclimate vars
##################################################

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% select("CUADRICULA")

utm1 <- st_read('01_DATA/INPUT/CUTM1_aquada/CUTM1x1Ext_AQUADA.shp') %>% select("UTMCODE1X1")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariablesWldClimPresente.xlsx'))
  
y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x %>% select("CUADRICULA", starts_with("Bio")), y)

xy <- xy %>% select(-"CUADRICULA") %>% rename("sp"="AQUADA")

aa <- data.frame(CUADRICULA=x$CUADRICULA, sp=xy$sp)

aa.sf <- inner_join(utm10, aa, by="CUADRICULA")



ext.raster <- rast(x = '01_DATA/INPUT/extremadura_10.tiff')  %>% 
  scale(center=TRUE, scale=TRUE) # present

names(ext.raster)  <- c(paste0("Bio0", 1:9), paste0("Bio", 10:19))


## Section: load future bioclimate variables
##################################################

load_list_of_tif <- function (path, period) {
  x <- list.files(path = here::here(paste0(path_to_tif, period)), 
             pattern = "*.tif", full.names = TRUE)
  return(x)
}

path_to_tif <- "01_DATA/INPUT/bioclim_future/ext10KM/"

tiff_files_2021_2040 <- load_list_of_tif(path_to_tif, "2021-2040")
tiff_files_2041_2060 <- load_list_of_tif(path_to_tif, "2041-2060")
tiff_files_2061_2080 <- load_list_of_tif(path_to_tif, "2061-2080")
tiff_files_2081_2100 <- load_list_of_tif(path_to_tif, "2081-2100")

resample_and_scale_tif <- function(tif_file_path, from_raster) {
  r <- terra::rast(tif_file_path) # load raster
  names(r)  <- c(paste0("Bio0", 1:9), paste0("Bio", 10:19))
  r <- terra::resample(r, from_raster)  %>% 
    scale(center=TRUE, scale=TRUE)
  return(r)
}

tiff_files_2021_2040_scaled <- lapply(tiff_files_2021_2040, resample_and_scale_tif, ext.raster)
tiff_files_2041_2060_scaled <- lapply(tiff_files_2041_2060, resample_and_scale_tif, ext.raster)
tiff_files_2061_2080_scaled <- lapply(tiff_files_2061_2080, resample_and_scale_tif, ext.raster)
tiff_files_2081_2100_scaled <- lapply(tiff_files_2081_2100, resample_and_scale_tif, ext.raster)


## Get preditors from presence / abssence data
##################################################

presvals <- extract(ext.raster, aa.sf %>% st_transform(st_crs(ext.raster)) %>% st_centroid() %>% st_coordinates())

xy.scaled <- cbind(presvals, sp=aa.sf$sp)


##################################################
## Section: Model
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


## Model 2: by Feature selection
##################################################

# FDR

xy.scaled.selection <- FDR(data=xy.scaled,sp.cols = 20, var.cols = 1:19)

bioclim_predictors <- sort(rownames(xy.scaled.selection$select))

xy.scaled.selection <- (xy.scaled  %>% dplyr::select(all_of(bioclim_predictors), "sp"))

# Correlation

cor_m <- cor(xy.scaled.selection %>% select(!all_of("sp")))
index <- sort(findCorrelation(cor_m))

xy.scaled.selection <- xy.scaled.selection[,  -index] # remove correlated variables

# xy.scaled.selection <- cbind(xy.scaled.selection, sp=xy$sp)

null.model <- glm(formula = sp ~ 1, family = binomial, data=xy.scaled.selection)


fmod_aa_biohist <- step(null.model, direction = "forward", 
                        keep = function(model, aic) list(model = model, aic = aic),
                        scope = (~Bio03 + Bio04 + Bio05 + Bio06 + Bio14 + Bio16))


summary(fmod_aa_biohist)

##################################################
## Section:  Historical - Prediction and Favourability
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

extract_fav_from_grid <- function(favourability_raster, umt_grid) {
  
  fav.rast <-  exact_extract(favourability_raster, umt_grid, 'mean', append_cols  =TRUE)
  
  utm10.fav <- inner_join(umt_grid, fav.rast, by="CUADRICULA")
  
  utm10.fav$category <- cut(utm10.fav$mean, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  
  ggplot(data=utm10.fav) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
    geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(sp==1), size=0.5)
  
}

plot(historical_prediction(ext.raster, bmod_aa_biohist))
plot(historical_prediction(ext.raster, fmod_aa_biohist))

epsg_code <- paste0("epsg:",st_crs(utm10)$epsg)

plot(historical_favourability(historical_prediction(ext.raster, bmod_aa_biohist), 
                              xy.scaled$sp, crs=epsg_code))

extract_fav_from_grid(historical_favourability(historical_prediction(ext.raster, bmod_aa_biohist), 
                                               xy.scaled$sp, crs=epsg_code),utm10)


#####################################################
## Section: bioclimatic model -  future
#####################################################

future_prediction <- function(bioclim_predictor, model) {
  pre <- predict(object = bioclim_predictor, model=model,
                 type="response")
  return(pre)
}

future_favourability<- function(future_prediction, presence_absence ) {
  f <-f <- Fav(obs = presence_absence, pred = future_prediction)
  f
}


favourability_to_vector <- function(favourability_raster) {
  fav.rast <-  exact_extract(favourability_raster, utm10, 'mean', append_cols  =TRUE)
  utm10.fav <- inner_join(utm10, fav.rast, by="CUADRICULA")
  utm10.fav$category <- cut(utm10.fav$mean, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  utm10.fav <- utm10.fav %>% fill( category, .direction = 'updown' ) # na values due to grid resolution
  return(utm10.fav)
}

plot_favourability <- function(favourability_vector) {
  ggplot(data=favourability_vector) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
    geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(sp==1)) +
    theme_void() + 
    theme(legend.position ="bottom")
}


##################################################
## Section: Future prediction and favourability
##################################################

## 2021 - 2040
##################################################

prediction_2021_2040 <-  lapply(tiff_files_2021_2040_scaled, future_prediction, bmod_aa_biohist)

favourability_2021_2040 <-  lapply(prediction_2021_2040, future_favourability, xy.scaled$sp)

favourability_vector_2021_2040 <- lapply(favourability_2021_2040, favourability_to_vector)

favourability_2021_2040_plot <- lapply(favourability_vector_2021_2040, plot_favourability)

do.call("grid.arrange", c(favourability_2021_2040_plot, ncol = 2))

## 2040 - 2061
##################################################

prediction_2041_2060 <-  lapply(tiff_files_2041_2060_scaled, future_prediction, bmod_aa_biohist)

favourability_2041_2060 <-  lapply(prediction_2041_2060, future_favourability, xy.scaled$sp)

favourability_vector_2041_2060 <- lapply(favourability_2021_2040, favourability_to_vector)

favourability_2041_2060_plot <- lapply(favourability_vector_2041_2060, plot_favourability)

do.call("grid.arrange", c(favourability_2041_2060_plot, ncol = 2))

## 2061 - 2080
##################################################

prediction_2061_2080 <-  lapply(tiff_files_2061_2080_scaled, future_prediction, bmod_aa_biohist)

favourability_2061_2080 <-  lapply(prediction_2061_2080, future_favourability, xy.scaled$sp)

favourability_vector_2061_2080 <- lapply(favourability_2061_2080, favourability_to_vector)

favourability_2061_2080_plot <- lapply(favourability_vector_2061_2080, plot_favourability)

do.call("grid.arrange", c(favourability_2061_2080_plot, ncol = 2))


## 2081 - 2100
##################################################


prediction_2081_2100 <-  lapply(tiff_files_2081_2100_scaled, future_prediction, bmod_aa_biohist)

favourability_2081_2100 <-  lapply(prediction_2081_2100, future_favourability, xy.scaled$sp)

favourability_2081_2100_plot <- lapply(favourability_2081_2100, plot_favourability)


do.call("grid.arrange", c(favourability_2081_2100_plot, ncol = 2))



