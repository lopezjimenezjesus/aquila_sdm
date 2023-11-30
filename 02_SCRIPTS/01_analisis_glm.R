##################################################
## Project: PAISAJE
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


here::i_am('Aquila_adalberti_SDM.Rproj')


# x <- st_read(dsn = here::here('01_DATA/INPUT/shp/CUTM10_AQUADA.shp')) %>% select('CUADRICULA', 'AQUADA') %>% st_drop_geometry()
# head(x)

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")


xy <- right_join(x,y, by="CUADRICULA")

paste(names(xy)[1:107], collapse = " + ")
# remove Lat Long variables


## Section: Inspect
##################################################

head(xy)

summary(xy)

# check relationships 

ggplot(xy) + 
  geom_point(aes(CazaMe, Conejo)) # remove CazaMe and keep conejo


ggplot(xy) + 
  geom_point(aes(sqrt(CazaMe), sqrt(Conejo))) 

ggplot(xy) + 
  geom_point(aes(ALT_MAX, ALT_RANGE))


ggplot(xy) + 
  geom_point(aes(ALT_MAX, ALT_RANGE))


ggplot(xy) + 
  geom_point(aes(ALT_MIN, ALT_RANGE))

ggplot(xy) + 
  geom_point(aes(Tmax1, Tmax7))




## Section: select initial predictors
##################################################



xy.varselected <- xy %>% dplyr::select(-c("CUADRICULA", "Ori_W", "Ori_S",
                              "La", "Lo", "Lo2", "La2", "LaLoVaria",
                              "S", "SE", "SW", "W", "N", "N", "NW", "NE", "E",
                              "DenGranito", "DenArena", "DenCaliza", "DenPizarra", "CazaMe"))

corrplot(cor(xy.varselected %>% dplyr::select(-"AQUADA")))


##################################################
## Section: Model 1 - Benchmark
##################################################


model1.aguila <- glm(AQUADA ~ .,
                     family = binomial, 
                     data = xy)


summary(model1.aguila)

##################################################
## Section: Model 2.1
##################################################

xy.varselected <- xy %>% 
  dplyr::select(-c("CUADRICULA", "TABSMAX1","TABSMAX7", "TABSMIN1", 
                   "TABSMIN7", "RAINMAX1", "RAINMAX7", "RAINDAY1", 
                   "DenPizarra", "DenArena", "DenGranito", "DenCaliza",
                   "CazaMa", "CazaMe","ALT_RANGE", "Tmax1"))


# model

model2.aguila <- glm(AQUADA ~ .,
                     family = binomial, 
                     data = xy.varselected)

summary(model2.aguila)

AUC(model = model2.aguila)

crosspred_glm <- mecofun::crossvalSDM(model2.aguila, traindat= xy.varselected,
                                      colname_pred=colnames(xy.varselected)[1:90], colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy.varselected$AQUADA, crosspred_glm)

##################################################
## Section: Model 2.2
##################################################


## Section: address multicollinearity 
#####################################

mult1 <-multicol(model = model2.aguila)  # multicollinearity based on the fitted model
mult1

model.subset <-  rownames(subset(mult1, VIF < 3)) # remove > 3

xy.varselected <- xy %>% dplyr::select(all_of(c(model.subset, "AQUADA")))

fdr.init <- FDR(xy.varselected, sp.cols = 14, var.cols = 1:13)

xy.varselected <-  xy %>% dplyr::select(all_of(c(rownames(fdr.init$select), "AQUADA")))

corrplot(cor(xy.varselected %>% dplyr::select(-"AQUADA")))

paste(colnames(xy.varselected), collapse = " + ")


model21.aguila <- step(glm(AQUADA ~ Ciervo + QUESUR + DenCap18 +  NumPoblaD + 
                             EUCSPP + DistPobla + LongElectD + CASSAT + 
                             Jabali + DenHidro + Conejo, 
                           family = binomial, data = xy))

summary(model21.aguila)


AUC(model = model21.aguila)


crosspred_glm <- mecofun::crossvalSDM(model21.aguila, traindat= xy.varselected,
                                      colname_pred=colnames(xy.varselected)[1:11], colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy.varselected$AQUADA, crosspred_glm)


## Section: Favourability
##################################################

fav <- Fav(model21.aguila) # note favourability is only applicable to unweighted predictons

aa <- data.frame(CUADRICULA=x$CUADRICULA, xy$AQUADA, f=round(fav, 4))

aa.sf <- left_join(utm10, aa, by="CUADRICULA")

aa.sf$category <- cut(aa.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))

## Section: Plot
##################################################

ggplot(data=aa.sf) +
  geom_sf(aes(fill=f)) +
  scale_fill_stepsn(colours=  c("#FF0000", "yellow", "green"),breaks= c(0.2, 0.8, 1))

ggplot(data=aa.sf) +
  geom_sf(aes(fill=category)) +
  scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
  geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(xy.AQUADA==1))



##################################################
## Section: Model 3
##################################################


## Section: glm scaled
##################################################

#FDR
names(xy) # ver en qué columnas están las especies y variables

xy.scaled <- scale(xy %>% dplyr::select(-"CUADRICULA"))
xy.scaled <- as.data.frame(xy.scaled)

# checked and scaling not required, use xy
xy.scaled.fdr <- FDR(xy.scaled, sp.cols = dim(xy.scaled)[-1], var.cols = 1:(dim(xy.scaled)[2]-1))


paste(rownames(xy.scaled.fdr$select), collapse = " + ")




#Regresión por pasos hacia delante
# modelo.aguila <- step(glm(AQUADA ~  QUESUR + Ciervo + CazaMa + DenPizarra + 
#                             TABSMAX1 + RAINDAY1 + QUEILE + P_DIAS + EUCSPP + 
#                             Deh + PSpr + DenGranito + Cul_len + Ptot + P_prim + 
#                             Cul_het + LongCarrD + DenCap18 + DistPobla + Psum + 
#                             Pwin + Jabali + Reg + Pene + Sup_arti + CazaMe + 
#                             Pvar + Conejo + Pjul + LongElectD + QFAGPY,
#                           family = binomial, data = xy),direction = "backward", keep =  function(model, aic) list(model = model, aic = aic))

model3.aguila <- step(glm(AQUADA ~  QUESUR + Ciervo + CazaMa + DenPizarra + TABSMAX1 +
                            RAINDAY1 + QUEILE + P_DIAS + EUCSPP + Deh + PSpr + DenGranito + 
                            Cul_len + NE + Ptot + P_prim + Cul_het + LongCarrD + DenCap18 +
                            DistPobla + W + Psum + Pwin + Jabali + Reg + Pene + Sup_arti + 
                            CazaMe + Pvar + Conejo + Pjul + LongElectD,
                          family = binomial, data = xy),direction = "backward", 
                      keep =  function(model, aic) list(model = model, aic = aic))


model3.aguila  # muestra el modelo
summary(model3.aguila)

xy.scaled$AQUADA <- xy$AQUADA

xy.varselected <- xy.scaled %>% dplyr::select(dput(rownames(xy.scaled.fdr$select)), "AQUADA")

crosspred_glm <- mecofun::crossvalSDM(model3.aguila, traindat= xy.varselected,
                                      colname_pred=colnames(xy.varselected)[1:32], colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy.varselected$AQUADA, crosspred_glm)



# remove not significative variables

model3.aguila <- modelTrim(model3.aguila)  # elimina las no significativas
model3.aguila

summary(model3.aguila)


crosspred_glm <- mecofun::crossvalSDM(model3.aguila, traindat= xy.varselected,
                                      colname_pred=colnames(xy.varselected)[1:32], colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy.varselected$AQUADA, crosspred_glm)


#Ecuaciones finales del modelo

getModEqn(model3.aguila, digits = 4)  # para redondear a menos d?gitos
getModEqn(model3.aguila, digits = 4, type = "P")
getModEqn(model3.aguila, digits = 4, type = "F")

plotGLM(model3.aguila, main = "Modelo con variables FDR")


# Check model

# Área bajo la curva
AUC(model = model3.aguila)

#Capacidad de clasificación
threshMeasures(model = model3.aguila, ylim = c(0, 1), thresh = "preval")

# Para encontrar el valor umbral que optimiza cada medida de discriminación, se puede usar la función optiThresh:
optiThresh(model = model3.aguila)
# También se puede optimizar el umbral para la mejor combinación de una pareja de medidas de evaluación relacionadas:
optiPair(model = model3.aguila)

# por defecto las dos medidas son sensibilidad y especificidad
# pero se puede especificar otras:
optiPair(model = model3.aguila, measures = c("Omission", "Commission"))
optiPair(model = model3.aguila, measures = c("UPR", "OPR"))

# 4.4. Analizar la calibración de un modelo
# Las medidas de calibración miden el ajuste de las probabilidades predichas por el modelo a las frecuencias observadas en los datos. La función HLfit calcula la bondad de ajuste con el estadístico de Hosmer-Lemeshow.
# La significación se refiere a la diferencia entre observado y predicho; por tanto, cuanto más alta la p, mejor. Esta función aporta también la raiz cuadrada de la media de los cuadrados de los errores (RMSE). 

HLfit(model = model3.aguila, bin.method = "prob.bins")

## Section: Favourability
##################################################

# pred <- predict(model3.aguila, xy.scaled[1:(length(xy.scaled)-1)], type = "response")
# columns <- attr(model3.aguila$terms, "term.labels")
# 
# pred <- predict(model3.aguila, xy.varselected %>% dplyr::select(all_of(columns)), type = "response")


fav <- Fav(model3.aguila) # note favourability is only applicable to unweighted predictons

aa <- data.frame(CUADRICULA=x$CUADRICULA, xy$AQUADA, f=round(fav, 4))

aa.sf <- left_join(utm10, aq, by="CUADRICULA")

aa.sf$category <- cut(aa.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))

## Section: Plot
##################################################

ggplot(data=aa.sf) +
  geom_sf(aes(fill=f)) +
  scale_fill_stepsn(colours=  c("#FF0000", "yellow", "green"),breaks= c(0.2, 0.8, 1))

ggplot(data=aa.sf) +
  geom_sf(aes(fill=category)) +
  scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
  geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(xy.AQUADA==1))

## Section: save model
##################################################


write_rds(model3.aguila, file = here('04_RESULTS/glm_model.rds'))


##################################################
## Section: Model 4
##################################################

#  Same as model 3 but Not scaled

#FDR
names(xy)  # ver en qué columnas están las especies y variables

# checked and scaling not required, use xy
xy.fdr <- FDR(xy, sp.cols = dim(xy)[-1], var.cols = 1:dim(xy)[2]-1)


paste(rownames(xy.fdr$select), collapse = " + ")


#Regresión por pasos hacia delante
modelo.aguila.notsc <- step(glm(AQUADA ~ Cul_len + DenPizarra + TABSMAX1 + DenGranito + CazaMa + Ciervo +
                            PSpr + QUESUR + Reg + Cul_het + Sup_arti + RAINDAY1 + DenCap18 + 
                            QFAGPY + P_DIAS + QUEILE + Ptot + P_prim + NumPoblaD + Deh +
                            LongCarrD + Psum + DenPobla + EUCSPP + Pwin + Arroz + 
                            SinVeg + DistPobla + Pene + Pjul + LongElectD + CASSAT + 
                            Jabali + Pvar + Tmax7 + Paut + RAINDAY7 + CazaMe,
                          family = binomial, data = xy))
modelo.aguila.notsc  # muestra el modelo
summary(modelo.aguila.notsc)

#Ecuaciones finales del modelo
getModEqn(modelo.aguila.notsc, digits = 4)  # para redondear a menos d?gitos
getModEqn(modelo.aguila.notsc, digits = 4, type = "P")
getModEqn(modelo.aguila.notsc, digits = 4, type = "F")

plotGLM(modelo.aguila.notsc, main = "Modelo con variables FDR")


# Check model

#?rea bajo la curva
AUC(model = modelo.aguila.notsc)

#Capacidad de clasificación
threshMeasures(model = modelo.aguila.notsc, ylim = c(0, 1), thresh = "preval")

# Para encontrar el valor umbral que optimiza cada medida de discriminación, se puede usar la función optiThresh:
optiThresh(model = modelo.aguila.notsc)
# También se puede optimizar el umbral para la mejor combinación de una pareja de medidas de evaluación relacionadas:
optiPair(model = modelo.aguila.notsc)

# por defecto las dos medidas son sensibilidad y especificidad
# pero se puede especificar otras:
optiPair(model = modelo.aguila.notsc, measures = c("Omission", "Commission"))
optiPair(model = modelo.aguila.notsc, measures = c("UPR", "OPR"))


## Section: Favourability
##################################################

red <- predict(modelo.aguila.notsc, xy[1:(length(xy)-1)], type = "response")

fav <- Fav(modelo.aguila.notsc) # note favourability is only applicable to unweighted predictons

aquila.adalberti_2 <- data.frame(CUADRICULA=x$CUADRICULA, xy$AQUADA, f=round(fav, 4))

aquila.adalberti_2.sf <- left_join(utm10, aquila.adalberti_2, by="CUADRICULA")

aquila.adalberti_2.sf$category <- cut(aquila.adalberti_2.sf$f, 
                                      breaks=c(-Inf, 0.2, 0.8, Inf), 
                                      labels=c("low","middle","high"))

## Section: Plot
##################################################

ggplot(data=aquila.adalberti_2.sf) +
  geom_sf(aes(fill=f)) +
  scale_fill_stepsn(colours=  c("#FF0000", "yellow", "green"),breaks= c(0.2, 0.8, 1))

ggplot(data=aquila.adalberti_2.sf) +
  geom_sf(aes(fill=category)) +
  scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
  geom_sf(data=st_centroid(aquila.adalberti_2.sf) %>% dplyr::filter(xy.AQUADA==1))


## Section:
##################################################

names(modelo.aguila$coefficients)

crosspred_glm <- mecofun::crossvalSDM(modelo.aguila, traindat= xy.varselected, 
                                      colname_pred=colnames(xy.varselected)[1:11], colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy.varselected$AQUADA, crosspred_glm)


##################################################
## Section: Model 5
##################################################

colnames(xy)
xy.sub <- xy %>% dplyr::select(c("Tsum",
                          "Pjul", "Paut", "PSpr", "Psum", "RadSol", "P_DIAS",
                          "RAINMAX1", "RAINMAX7", "DenPobla", "DistPobla", "NumPoblaD", 
                          "DistCarre", "Reg",
                          "PD", "SHDI", "CazaMa",
                          "QFAGPY", "QUESUR", "QUEILE", "CASSAT", "EUCSPP", "AltVeg", "Conejo", "AQUADA"
))

xy.fdr <- FDR(xy.sub, sp.cols = dim(xy.sub)[2], var.cols = 1:dim(xy.sub)[2]-1)
paste(rownames(xy.fdr$select), collapse = " + ")



# remove correlated variables

#Regresión por pasos hacia delante
model5.aguila <- step(glm(AQUADA ~  PSpr + QUESUR + Reg + QFAGPY + P_DIAS + 
                                   QUEILE + NumPoblaD + Psum + DenPobla + EUCSPP + 
                                   DistPobla + Pjul + CASSAT + Paut + Conejo + RAINMAX7,
                                 family = binomial, data = xy))

summary(model5.aguila)

#?rea bajo la curva
model5.aguila <- modelTrim(model5.aguila) 

AUC(model = model5.aguila)


pred <- predict(model5.aguila, xy.sub[-1], type = "response")
fav <- Fav(model5.aguila) # note favourability is only applicable to unweighted predictons



aa<- data.frame(CUADRICULA=x$CUADRICULA, xy %>% dplyr::select("AQUADA", ), f=round(fav, 4))

aa.sf <- left_join(utm10, aquila.adalberti, by="CUADRICULA")

aa.sf$category <- cut(aa.sf$f, 
                   breaks=c(-Inf, 0.2, 0.8, Inf), 
                   labels=c("low","middle","high"))



ggplot(data=aa.sf) +
  geom_sf(aes(fill=f)) +
  scale_fill_stepsn(colours=  c("#FF0000", "yellow", "green"),breaks= c(0.2, 0.8, 1))

ggplot(data=aa.sf) +
  geom_sf(aes(fill=category)) +
  scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
  geom_sf(data=st_centroid(aa.sf) %>% dplyr::filter(xy.AQUADA==1))
