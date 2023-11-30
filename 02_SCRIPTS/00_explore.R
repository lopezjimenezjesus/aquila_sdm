##################################################
## Project: Aquila adalberti SDM 
## Script purpose: Explore variables
## Date:
## Author: Jesús Jiménez López
##################################################


##################################################
## Section: SET UP ----
##################################################


## Section: Libraries
##################################################


library(fuzzySim)
library(modEvA)
library(foreign)
library(here)
library(sf)
library(readxl)
library(tidyverse)
library(skimr)
library(corrplot)
library(GGally)
library(gridExtra)

here::i_am('Aquila_adalberti_SDM.Rproj')


##################################################
## Section: DATA ----
##################################################

## Section: Load
##################################################


utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x,y, by="CUADRICULA")

head(xy)

colnames(xy)[dim(xy)[2]] <- "occ" # occurrence

glimpse(xy)

skim

paste(names(xy)[1:107], collapse = " + ")

## Section: Create a metadata table with predictors 
##################################################

## Join metadata table with predictors from dataset
## useful to retrieve predictors based on Paisaje, Habitat and NestSite

# read predictors metadata from csv 

metadata_predictors <- readr::read_delim('01_DATA/INPUT/metadatos_variables_predictoras.csv', delim = ";", locale=locale(encoding="latin1")) 

colnames(metadata_predictors) <- c("tipologia", "nombre", "descripcion", "A_Paisaje", "B_Habitat", "C_NestSite")

metadata_predictors <- dplyr::add_row(
  metadata_predictors,
  tipologia = "VEGETACIÓN",
  nombre = "Pino",
  descripcion = "% de FCC de Pinus (8)",
  A_Paisaje = "x",
  B_Habitat = "x",
  C_NestSite = "x",
  .before = 91
)


# get predictors names from dataset and add missing predictors names

predictors_names <- c(colnames(xy)[5:(length(colnames(xy))-1)], c("AlturaSup", "AlturaInf"))

predictors_names <- c("ALT_MAX", "ALT_MIN", "ALT_MED", "ALT_RANGE", "AlturaSup", "AlturaInf",
                      "Tri_mean", "Slope_mean", "Ori_W", "Ori_S", "E", "N", "NE", "NW", 
                      "S", "SE", "SW", "W", "La", "Lo", "La2", "Lo2", "LaLoVaria", 
                      "Tmed", "Rmtd", "Isot", "Test", "Tmax7", "Tmax1", "Tran", "Ptot", "Pvar", 
                      "P_prim", "Taut", "TSpr", "Tsum", "Twin", "Pwin", "Pene", "Pjul", 
                      "Paut", "PSpr", "Psum", "RadSol", "P_DIAS", "TABSMAX1", "TABSMAX7", 
                      "TABSMIN1", "TABSMIN7", "FROSTDAY", "RAINMAX1", "RAINMAX7", "RAINDAY1", "RAINDAY7", 
                      "DenPobla", "DistPobla", "NumPoblaD", "DistCarre", "LongCarrD", 
                      "DistCamin", "LongCaminD", "DistElect", "LongElectD", 
                      "Arroz", "Cul_sec", "Cul_len", "Prad", "Cul_het", "Deh", "Bosq", 
                      "Mat", "Agu_cont", "Pas_nat", "Reg", "Sup_arti", "SinVeg", "NP", 
                      "PD", "LPI", "LSI", "AREA_MN", "FRAC_AM", "CONTAG", "SHDI", "SHEI", 
                      "QFAGPY", "QUESUR", "QUEILE", "CASSAT", "EUCSPP", "PINO", "AltVeg", 
                      "DenCap18", "DenOvi19", "DenPor19", "DenVac19", "CazaMa", "CazaMe", 
                      "Conejo", "Perdiz", "Ciervo", "Jabali", "DenArena", "DenCaliza", "DenGranito", 
                      "DenPizarra", "DenHidro", "DenEmbalse"
)

# bind both tables

predictors_table <- cbind(metadata_predictors, predictors_names)

predictors_table <-predictors_table %>%
  mutate(
    factores = ifelse(
      tipologia %in% c("TOPOGRÁFICAS", "ESPACIALES", "CLIMÁTICAS", "ÍNDICES DEL PAISAJE","LITOLOGÍA", "HIDROLOGÍA"), 
      "Factores Ambientales Abióticos", 
      ifelse(
        tipologia %in% c("ACTIVIDAD HUMANA"), "Factores Antrópicos",
        ifelse(
          tipologia %in% c("USOS DEL SUELO", "VEGETACIÓN", "GANADO Y CAZA"), "Factores Ambientales Bióticos", ""
          )
        )
      )
    )

readr::write_csv(predictors_table, file = "01_DATA/OUTPUT/metadata_predictors.csv")

skim(xy)


# convert climatic variables v/10

reduce_scaling <- function(x, na.rm = FALSE) (x/10)


predictors_to_scale <- predictors_table[predictors_table$tipologia=="CLIMÁTICAS","predictors_names"]

predictors_to_scale <- predictors_to_scale[! predictors_to_scale %in% c('RadSol', 'P_DIAS', "FROSTDAY", "RAINDAY1", "RAINDAY7")]

xy <- xy %>% mutate_at(predictors_to_scale, .funs = reduce_scaling)

##################################################
## Section: FUNCTIONS ----
##################################################


filter_predictors_by_type_analysis <- function(type_analisys="A_Paisaje", name_of_columns=c("tipologia","predictors_names")) {
  
  predictors_table[!is.na(predictors_table[type_analisys]),c(name_of_columns)] 
}


select_variables_by_type__class <- function(type_class="TOPOGRÁFICAS") {
  xy %>% dplyr::select(a_paisaje[a_paisaje$tipologia==type_class,]$predictors_names, "occ")
}




plot_density <- function(dataset, predictor_variable, presence_variable) {
  
  
  ggplot(dataset, aes(x = .data[[predictor_variable]])) +
    geom_density(alpha = .2, fill = "#FF6666") +
    geom_point(aes(col=as.factor(.data[[presence_variable]]), y=0)) + 
    labs(title="", x=predictor_variable, color='Presence') +
    facet_grid(.data[[presence_variable]]~.) + 
    theme_light()
  
}


plot_all_density <- function(dataset) {
  myplots <- vector('list', ncol(dataset))
  
  for (i in seq(from=1, to=dim(dataset)[2]-1)) {
    message(i)
    myplots[[i]] <- local({
      i <- i
      p1 <- plot_density(dataset, colnames(dataset)[i], "occ")
      print(p1)
    })
  }
  
  
  # Grid 
  
  n <- length(myplots)
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(myplots, ncol=nCol))
  
  
}


fdr_var_selection <- function(dataset) {
  
  FDR(dataset, sp.cols = dim(dataset)[2], var.cols = 1:(dim(dataset)[2]-1))
}


fdr_var_subset  <- function(fdr_dataset) {
  
  selected <-  rownames(fdr_dataset$select)
  
  xy %>% dplyr::select(all_of(selected))
  
}



fdr_var_subset  <- function(fdr_dataset) {
  
  selected <-  rownames(fdr_dataset[["select"]])
  
  xy.selected <- xy %>% dplyr::select(all_of(selected))
  
  xy.selected
}


get_model_variables <- function(dataset, init, end) {
  paste(names(dataset)[init:end], collapse = " + ")
}

##################################################
## Section: PREDICTOR VARIABLES SELECTION 1
##################################################


##################################################
## Section: A-Paisaje
##################################################

a_paisaje <- filter_predictors_by_type_analysis()

write_csv(a_paisaje, file = "01_DATA/OUTPUT/pasiaje_predictors.csv")


## TOPOGRÁFICAS
##################################################

xy.topo_paisaje <- select_variables_by_type__class()

# Density plot

plot_all_density(xy.topo_paisaje)

# Corplot

xy.topo_paisaje.cor <- cor(xy.topo_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.topo_paisaje.cor)

ggcorr(xy.topo_paisaje.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.topo_paisaje.fdr <- fdr_var_selection(xy.topo_paisaje)

## ESPACIALES
##################################################

xy.espaciales_paisaje <- select_variables_by_type__class(type_class = "ESPACIALES")

# Density plot

plot_all_density(xy.espaciales_paisaje)

# Corplot

xy.espaciales.cor <- cor(xy.espaciales_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.espaciales.cor)

ggcorr(xy.espaciales.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.espaciales_paisaje.fdr <- fdr_var_selection(xy.espaciales_paisaje)

## CLIMÁTICAS
##################################################

xy.climaticas_paisaje <- select_variables_by_type__class(type_class = "CLIMÁTICAS")

# Density plot

plot_all_density(xy.climaticas_paisaje)

# Corplot

xy.climaticas.cor <- cor(xy.climaticas_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.climaticas.cor)

ggcorr(xy.climaticas.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")



# fdr 

xy.climaticas.fdr <- fdr_var_selection(xy.climaticas_paisaje)



## ACTIVIDAD HUMANA
##################################################

xy.humana_paisaje <- select_variables_by_type__class(type_class = "ACTIVIDAD HUMANA")

# Density plot

plot_all_density(xy.humana_paisaje)

# Corplot

xy.act_humana.cor <- cor(xy.humana_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.act_humana.cor)

ggcorr(xy.act_humana.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")



# fdr 

xy.act_humana.fdr <- fdr_var_selection(xy.humana_paisaje)



## USOS DEL SUELO
##################################################


xy.usos_paisaje <- select_variables_by_type__class(type_class = "USOS DEL SUELO")


# Density plot

plot_all_density(xy.usos_paisaje)

# Corplot

xy.usos_suelo.cor <- cor(xy.usos_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.usos_suelo.cor)

ggcorr(xy.usos_suelo.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.usos_suelo.fdr <- fdr_var_selection(xy.usos_paisaje)


## ÍNDICES DEL PAISAJE
##################################################

xy.indices_paisaje <- select_variables_by_type__class(type_class = "ÍNDICES DEL PAISAJE")


# Density plot

plot_all_density(xy.indices_paisaje)


# Corplot

xy.paisaje.cor <- cor(xy.indices_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.paisaje.cor)

ggcorr(xy.paisaje.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.paisaje.fdr <- fdr_var_selection(xy.indices_paisaje)



## VEGETACIÓN
##################################################


xy.vegetacion_paisaje <- select_variables_by_type__class(type_class = "VEGETACIÓN")

# Density plot

plot_all_density(xy.vegetacion_paisaje)


# Corplot

xy.vegetacion.cor <- cor(xy.vegetacion_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.vegetacion.cor)

ggcorr(xy.vegetacion_paisaje,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.vegetacion.fdr <- fdr_var_selection(xy.vegetacion_paisaje)



## GANADO Y CAZA
##################################################

xy.ganado_caza_paisaje <- select_variables_by_type__class(type_class = "GANADO Y CAZA")

# Density plot

plot_all_density(xy.ganado_caza_paisaje)


# Corplot

xy.ganado_caza.cor <- cor(xy.ganado_caza_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.ganado_caza.cor)

ggcorr(xy.ganado_caza.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.ganado_caza.fdr <- fdr_var_selection(xy.ganado_caza_paisaje)


## LITOLOGÍA
##################################################

# No hay en paisaje

## HIDROLOGÍA
##################################################

# No hay en paisaje


##################################################
## Section: PREDICTOR VARIABLE SELECTION 2
##################################################


## PAISAJE
##################################################

# exclude

# xy.litologia.fdr 
# xy.espaciales.fdr
# xy.hidrologia.fdr
# xy.litologia.fdr

# new set

l_fdr_paisaje <- list(xy.topo_paisaje.fdr, xy.climaticas.fdr, xy.act_humana.fdr,
              xy.usos_suelo.fdr, xy.paisaje.fdr, xy.vegetacion.fdr, xy.ganado_caza.fdr)

l_fdr_paisaje_selected <- lapply(l_fdr_paisaje, fdr_var_subset)

paisaje_variables <- do.call("cbind", l_fdr_paisaje_selected)

xy_paisaje <- data.frame(paisaje_variables, occ= xy$occ)

xy_paisaje.fdr <- fdr_var_selection(xy_paisaje)


paisaje_variables <- xy_paisaje %>% dplyr::select(rownames(xy_paisaje.fdr$select))


ggcorr(cor(paisaje_variables),
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


##################################################
## Section: EXPLORE NAIVE MODELS
##################################################

## Section: Multicolinearity
##################################################

mult1 <- multicol(subset(xy_paisaje, select = 1:(dim(xy_paisaje)[2] -1)))

subset(mult1, VIF < 5)  # muestra cuales tienen VIF elevado

xy_paisaje <- xy %>% dplyr::select(rownames(subset(mult1, VIF < 5)))
xy_paisaje <- data.frame(xy_paisaje, occ= xy$occ)

xy_paisaje.scale <- as.data.frame(scale(xy_paisaje))
xy_paisaje.scale$occ <- xy_paisaje$occ

null.model = glm(occ ~ 1, data = xy_paisaje, family=binomial)

paste(colnames(xy_paisaje), collapse = " + ")

model <- step(null.model, direction='forward',
     keep =  function(model, aic) list(model = model, aic = aic),
     scope = (~Conejo + Tmax7 + RAINDAY1 + QUEILE + RAINDAY7 + Deh + LongCarrD +
                QFAGPY + SinVeg + DistPobla + Reg + Cul_len + NumPoblaD + CASSAT + 
                DenCap18 + TABSMAX1 + Cul_het + LongElectD + QUESUR + EUCSPP))

# scaled

null.model = glm(occ ~ 1, data = xy_paisaje.scale, family='binomial')

paste(colnames(xy_paisaje.scale), collapse = " + ")


model.scaled <- step(null.model, direction='forward',
              keep =  function(model, aic) list(model = model, aic = aic),
              scope = (~Conejo + Tmax7 + RAINDAY1 + QUEILE + RAINDAY7 + Deh + LongCarrD + 
                         QFAGPY + SinVeg + DistPobla + Reg + Cul_len + NumPoblaD + CASSAT + 
                         DenCap18 + TABSMAX1 + Cul_het + LongElectD + QUESUR + EUCSPP))





##################################################
## Section: B-Hábitat ----
##################################################


filter_predictors_by_type_analysis("B_Habitat")


##################################################
## Section: C-NestSite ----
##################################################

filter_predictors_by_type_analysis("C_NestSite")


