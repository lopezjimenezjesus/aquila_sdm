##################################################
## Project: Aquila adalberti SDM 
## Script purpose: Explore variables
## Date:
## Author: Jesús Jiménez López
##################################################


##################################################
## Section: SET UP
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
library(corrplot)
library(GGally)
library(gridExtra)

here::i_am('Aquila_adalberti_SDM.Rproj')


## Section: Data
##################################################


utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x,y, by="CUADRICULA")

paste(names(xy)[1:107], collapse = " + ")


## Functions
##################################################


plot_density <- function(dataset, predictor_variable, presence_variable) {
  
  
  ggplot(dataset, aes(x = .data[[predictor_variable]])) +
    geom_density(alpha = .2, fill = "#FF6666") +
    geom_point(aes(col=as.factor(.data[[presence_variable]]), y=0)) + 
    labs(title="Density plot", x=predictor_variable, color='Presence') +
    facet_grid(.data[[presence_variable]]~.)
  
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
  
  FDR(dataset, sp.cols = dim(dataset)[2], var.cols = 1:dim(dataset)[2]-1)
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

head(xy)

summary(xy)

colnames(xy)[dim(xy)[2]] <- "occ"

## TOPOGRÁFICAS
##################################################

## ALL

topo <- c("ALT_MAX", "ALT_MIN", "ALT_MED", "ALT_RANGE", "Tri_mean", "Slope_mean", 
          "Ori_W", "Ori_S",  "E", "N", "NE", "NW", 
          "S", "SE", "SW", "W")


xy.topo <- xy %>% dplyr::select(all_of(topo), "occ")

# Density plot

plot_all_density(xy.topo)

# Corplot

xy.topo.cor <- cor(xy.topo %>% dplyr::select(-"occ"))

corrplot(xy.topo.cor)

ggcorr(xy.topo.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.topo.fdr <- fdr_var_selection(xy.topo)


## PAISAJE

topo_paisaje <- c("ALT_MAX", "ALT_MIN", "ALT_MED", "ALT_RANGE", "Tri_mean", "Slope_mean")
xy.topo_paisaje <- xy %>% dplyr::select(all_of(topo_paisaje), "occ")

# Density plot

plot_all_density(xy.topo_paisaje)

# Corplot

xy.topo_paisaje.cor <- cor(xy.topo_paisaje %>% dplyr::select(-"occ"))

corrplot(xy.topo_paisaje.cor)

ggcorr(xy.topo.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.topo_paisaje.fdr <- fdr_var_selection(xy.topo_paisaje.cor)



## ESPACIALES
##################################################

espaciales <- c("La", "Lo", "La2", "Lo2", "LaLoVaria")

xy.espaciales <- xy %>% dplyr::select(all_of(espaciales), "occ")

# Density plot

plot_all_density(xy.espaciales)

# Corplot

xy.espaciales.cor <- cor(xy.espaciales %>% dplyr::select(-"occ"))

corrplot(xy.espaciales.cor)

ggcorr(xy.espaciales.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.espaciales.fdr <- fdr_var_selection(xy.espaciales)


## CLIMÁTICAS
##################################################

climaticas <- c("Tmed", "Rmtd", "Isot", "Test", "Tmax7", "Tmax1",  
                "Tran", "Ptot", "Pvar", "P_prim", "Taut", "TSpr", "Tsum", "Twin", "Pwin", "Pene",
                "Pjul", "Paut", "PSpr", "Psum", "RadSol", "P_DIAS", 
                "TABSMAX1", "TABSMAX7", "TABSMIN1", "TABSMIN7", 
                "FROSTDAY", "RAINDAY1", "RAINDAY7", "RAINMAX1", "RAINMAX7")


xy.climaticas <- xy %>% dplyr::select(all_of(climaticas), "occ")


# Density plot

plot_all_density(xy.climaticas)

# Corplot

xy.climaticas.cor <- cor(xy.climaticas %>% dplyr::select(-"occ"))

corrplot(xy.climaticas.cor)

ggcorr(xy.climaticas.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")



# fdr 

xy.climaticas.fdr <- fdr_var_selection(xy.climaticas)



## ACTIVIDAD HUMANA
##################################################

act_humana <- c( "DenPobla", "DistPobla", "NumPoblaD", "DistCarre", 
             "LongCarrD", "DistCamin", "LongCaminD", "DistElect", "LongElectD")

xy.act_humana <- xy %>% dplyr::select(all_of(act_humana), "occ")


# Density plot

plot_all_density(xy.act_humana)

# Corplot

xy.act_humana.cor <- cor(xy.act_humana %>% dplyr::select(-"occ"))

corrplot(xy.act_humana.cor)

ggcorr(xy.act_humana.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")



# fdr 

xy.act_humana.fdr <- fdr_var_selection(xy.act_humana)



## USOS DEL SUELO
##################################################


usos_suelo <- c("Arroz", "Cul_sec", "Cul_len", "Prad", "Cul_het", "Deh", "Bosq", 
                "Mat", "Agu_cont", "Pas_nat", "Reg", "Sup_arti", "SinVeg")


xy.usos_suelo <- xy %>% dplyr::select(all_of(usos_suelo), "occ")

# Density plot

plot_all_density(xy.usos_suelo)

# Corplot

xy.usos_suelo.cor <- cor(xy.usos_suelo %>% dplyr::select(-"occ"))

corrplot(xy.usos_suelo.cor)

ggcorr(xy.usos_suelo.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr 

xy.usos_suelo.fdr <- fdr_var_selection(xy.usos_suelo)


## ÍNDICES DEL PAISAJE
##################################################

paisaje <- c("NP", "PD", "LPI", "LSI", 
             "AREA_MN", "FRAC_AM", "CONTAG", "SHDI", "SHEI")


xy.paisaje <- xy %>% dplyr::select(all_of(paisaje), "occ")

# Density plot

plot_all_density(xy.paisaje)


# Corplot

xy.paisaje.cor <- cor(xy.paisaje %>% dplyr::select(-"occ"))

corrplot(xy.paisaje.cor)

ggcorr(xy.paisaje.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.paisaje.fdr <- fdr_var_selection(xy.paisaje)



## VEGETACIÓN
##################################################

vegetacion <- c("QFAGPY", "QUESUR", "QUEILE", "CASSAT", "EUCSPP", "PINO", "AltVeg")

xy.vegetacion <- xy %>% dplyr::select(all_of(vegetacion), "occ")


# Density plot

plot_all_density(xy.vegetacion)


# Corplot

xy.vegetacion.cor <- cor(xy.vegetacion %>% dplyr::select(-"occ"))

corrplot(xy.vegetacion.cor)

ggcorr(xy.vegetacion.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.vegetacion.fdr <- fdr_var_selection(xy.vegetacion)



## GANADO Y CAZA
##################################################

ganado_caza <- c("DenCap18", "DenOvi19", "DenPor19", "DenVac19", "CazaMa", 
                 "CazaMe", "Conejo", "Perdiz", "Ciervo", "Jabali")

xy.ganado_caza <- xy %>% dplyr::select(all_of(ganado_caza), "occ")


# Density plot

plot_all_density(xy.ganado_caza)


# Corplot

xy.ganado_caza.cor <- cor(xy.ganado_caza %>% dplyr::select(-"occ"))

corrplot(xy.ganado_caza.cor)

ggcorr(xy.ganado_caza.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.ganado_caza.fdr <- fdr_var_selection(xy.ganado_caza)


## LITOLOGÍA
##################################################


litologia <- c( "DenArena", "DenCaliza", "DenGranito", 
                "DenPizarra")



xy.litologia <- xy %>% dplyr::select(all_of(litologia), "occ")


# Density plot

plot_all_density(xy.litologia)


# Corplot

xy.litologia.cor <- cor(xy.litologia %>% dplyr::select(-"occ"))

corrplot(xy.litologia.cor)

ggcorr(xy.litologia.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")


# fdr

xy.litologia.fdr <- fdr_var_selection(xy.litologia)



## HIDROLOGÍA
##################################################

hidrologia <- c( "DenHidro", "DenEmbalse")

xy.hidrologia <- xy %>% dplyr::select(all_of(hidrologia), "occ")


# Density plot

plot_all_density(xy.hidrologia)


# Corplot

xy.hidrologia.cor <- cor(xy.hidrologia %>% dplyr::select(-"occ"))

corrplot(xy.hidrologia.cor)

ggcorr(xy.hidrologia.cor,
       method = c("pairwise", "spearman"),
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 3,
       color = "grey50")

# fdr

xy.hidrologia.fdr <- fdr_var_selection(xy.hidrologia)


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


## Section: 
##################################################


mult1 <- multicol(subset(xy_paisaje, select = 1:dim(xy_paisaje)[2] -1))  # especificar sÃ³lo las columnas con variables!

subset(mult1, VIF < 5)  # muestra cuÃ¡les tienen VIF elevado

xy_paisaje <- xy %>% dplyr::select(rownames(subset(mult1, VIF < 5)))  # muestra cuÃ¡les tienen VIF elevado
xy_paisaje <- data.frame(xy_paisaje, occ= xy$occ)

logit <- step(glm(occ~., data = xy_paisaje, family = 'binomial'), direction = "forward")

min.model = glm(occ ~ 1, data = xy_paisaje, family='binomial')

predictor_variables <- get_model_variables(xy_paisaje, 1, dim(xy_paisaje)[2]-1)


#statistical model

fwd.model = step(min.model, direction='forward',
                 keep = function(model, aic) list(model = model, aic = aic),
                 scope=(~Conejo + ALT_MIN + QUEILE + RAINDAY1 + Deh + RAINDAY7 + LongCarrD + QFAGPY + SinVeg + Reg + DistPobla + NumPoblaD + Cul_len + CASSAT + DenCap18 + TABSMAX1 + Cul_het + LongElectD + QUESUR + EUCSPP + Arroz))

# ecological model: select make sense variables

fwd.model = step(min.model, direction='forward',
                 keep = function(model, aic) list(model = model, aic = aic),
                 scope=(~Conejo + ALT_MIN + QUEILE +RAINDAY1 + Deh  + RAINDAY7 + 
                          LongCarrD + QFAGPY + SinVeg + Reg +   DistPobla + LongElectD + DenCap18 + QUESUR + EUCSPP))

summary(fwd.model)

# summary(fwd.model$keep[,2]$model)

## Section: Favourability
##################################################

Fav(fwd.model)

xy.fav <- data.frame(id=xy$CUADRICULA, occ=xy$occ, f=Fav(fwd.model))

xy.fav.utm10 <- left_join(utm10, xy.fav,join_by(CUADRICULA  == id))

plot(xy.fav.utm10)


xy.fav.utm10$category <- cut(xy.fav.utm10$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))

## Section: Plot
##################################################

ggplot(data=xy.fav.utm10) +
  geom_sf(aes(fill=f)) +
  scale_fill_stepsn(colours=  c("#FF0000", "yellow", "green"),breaks= c(0.2, 0.8, 1))

ggplot(data=xy.fav.utm10) +
  geom_sf(aes(fill=category)) +
  scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
  geom_sf(data=st_centroid(xy.fav.utm10) %>% dplyr::filter(occ==1))



table_mat_ <- table(xy.fav.utm10$occ, xy.fav.utm10$f > 0.2)
table_mat_

table_mat_ <- table(xy.fav.utm10$occ, xy.fav.utm10$f > 0.5)
table_mat_

## Section: Train test 
##################################################


set.seed(1234)
create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1: total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}
data_train <- create_train_test(xy_paisaje, 0.8, train = TRUE)
data_test <- create_train_test(xy_paisaje, 0.8, train = FALSE)
dim(data_train)
dim(data_test)

formula <- occ~.
logit <- step(glm(formula, data = data_train, family = 'binomial'), direction = "forward")
summary(logit)

predict_ <- predict(fwd.model, data_test, type = 'response')

fav_test <-  Fav(pred = predict_, obs=data_test$occ) # transform to fav

# confusion matrix
table_mat <- table(data_test$occ, fav_test > 0.5)
table_mat

accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test


precision <- function(matrix) {
  # True positive
  tp <- matrix[2, 2]
  # false positive
  fp <- matrix[1, 2]
  return (tp / (tp + fp))
}

recall <- function(matrix) {
  # true positive
  tp <- matrix[2, 2]# false positive
  fn <- matrix[2, 1]
  return (tp / (tp + fn))
}

prec <- precision(table_mat)
prec
rec <- recall(table_mat)
rec

f1 <- 2 * ((prec * rec) / (prec + rec))
f1

library(ROCR)
ROCRpred <- prediction(predict, data_test$sp)
ROCRperf <- performance(ROCRpred, 'tpr', 'fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7))



## check partial responses

# Load the mecofun package
library(mecofun)

# Names of our variables:
pred <- colnames(bio.sp.selected)[1:length(bio.sp.selected)-1]

# We want three panels next to each other:
par(mfrow=c(1,3)) 

# Plot the partial responses
partial_response(logit, predictors = bio.sp.selected[,pred])


## Section:
##################################################


# remove Lat Long variables

xy.sel <- xy %>% dplyr::select(c("Slope_mean", "Conejo", "CazaMa", "Deh", "Bosq", "Mat", "Agu_cont",
                              "QUESUR", "QUEILE", "DistPobla", "EUCSPP", "occ"))

data_train <- create_train_test(xy.sel, 0.8, train = TRUE)
data_test <- create_train_test(xy.sel, 0.8, train = FALSE)
dim(data_train)
dim(data_test)

formula <- occ~.
logit <- glm(formula, data = data_train, family = 'binomial')
summary(logit)



predict_ <- predict(logit, data_test, type = 'response')
# confusion matrix
table_mat <- table(data_test$occ, predict_ > 0.5)
table_mat

accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test


# LITOLOGÍA

# HIDROLOGÍA


ggplot(xy) + 
  geom_point(aes(CazaMe, Conejo)) # remove CazaMe and keep conejo

ggplot(xy) + 
  geom_point(aes(CazaMe, Conejo)) # remove CazaMe and keep conejo


ggplot(xy) + 
  geom_point(aes(sqrt(CazaMe), sqrt(Conejo))) 



ggplot(xy) + 
  geom_point(aes(Tmax1, Tmax7))


ggplot(xy, aes(x = Tmax1, y = Tmax7, col = as.factor(occ))) + 
  geom_point(alpha = .3) + 
  coord_equal(ratio = 1)



ggplot(xy, aes(x = CazaMa, y = CazaMe, col = as.factor(occ))) + 
  geom_point(alpha = .3) + 
  coord_equal(ratio = 0.25)




## Section: Compare bioclimate variables
##################################################


## Section: load sp and historical bioclimate vars
##################################################

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariablesWldClimPresente.xlsx'))

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x %>% select("CUADRICULA", starts_with("Bio")), y)

xy <- xy %>% select(-"CUADRICULA") %>% rename("sp"="AQUADA")


## Section: plot bioclim future - present
##################################################

# example

wc_ssp126_2021_2040 <- read.csv2(file =here::here('01_DATA/INPUT/bioclim_future/ext10KM/2021-2040/wc2.1_30s_bioc_ssp126_2021-2040_ext10KM.csv'),
                                 sep =",") 

colnames(wc_ssp126_2021_2040)  <- c("X", "CUADRICULA", paste0("Bio0", 1:9), paste0("Bio", 10:19))

xy.126_2021_2040 <- data.frame(wc_ssp126_2021_2040 %>% dplyr::select(-c("X", "CUADRICULA")) %>%   mutate(across(1:dim(.)[2], as.numeric)),
                               sp=xy$sp)

present <- cbind.data.frame(xy, period="present")
future <- cbind.data.frame(xy.126_2021_2040, period="126_2021_2040")
all <- rbind.data.frame(present %>% dplyr::select(-"sp"), future %>% dplyr::select(-c("sp")), stringsAsFactors = FALSE)

ggplot(data=all) +
  geom_boxplot(aes(x=period, y=Bio01, fill=period))

library(reshape2)
mm = melt(all, id=c('period'))

foo <- mm %>% group_by(period, variable)  %>% summarise(mean=mean(value))


ggplot(data=mm) +
  geom_boxplot(aes(x=variable, y=value, fill=period))

