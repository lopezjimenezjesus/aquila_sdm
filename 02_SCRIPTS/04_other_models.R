
#####################################################
## Section: Explore other models
#####################################################


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


## Section: 
##################################################

library(glmnet)

smod_aa_biohist <- (glm(sp ~ .,
                        family = binomial,  
                        data = xy))


lasso.model <- glmnet(x = xy %>% dplyr::select(-sp),
                y = as.factor(xy$sp),
                family = "binomial",  relax=TRUE,
                alpha = 1)
lasso_pred <- predict(lasso, newx=as.matrix(x[5:23], s = c(0.01, 0.005)))

plot(lasso.model, xvar="lambda")

## Section:
##################################################
library(dismo)



xy <- right_join(x %>% dplyr::select("CUADRICULA", starts_with("Bio")), y)

xy.sf <- inner_join(utm10, xy, by="CUADRICULA")

sp <- xy.sf %>% dplyr::filter(AQUADA==1) %>% st_coordinates()

xy <- xy %>% select(-"CUADRICULA") %>% rename("sp"="AQUADA")

bc <- bioclim(xy %>% dplyr::select(-"sp") %>% as.matrix(), xy$sp)

## RF
##################################################
library(randomForest)


aquada.rf <- randomForest(x = x[5:23],y = as.factor(xy$sp), ntree=101, proximity=TRUE, oob.prox=FALSE)

Fav(aquada.rf)

print(aquada.rf)

importance(aquada.rf)
varImpPlot(aquada.rf)

## Favourability
##################################################

fav_to_spatial_grid <- function(model, sp_df, grid_code_predictor, utm_grid) {
  # build a spatial grid based on favourability values
  fav <- Fav(model)
  df <- data.frame(CUADRICULA=grid_code_predictor[["CUADRICULA"]], occ=sp_df[["sp"]], f=fav)
  aq.sf <- inner_join(utm_grid, df, by="CUADRICULA")
  aq.sf$category <- cut(aq.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  return(aq.sf)
}


aq.sf <- fav_to_spatial_grid(aquada.rf, xy, x, utm10)

aq10cat.sf <- fav_to_spatial_grid_10cat(model.glm, xy, x, utm10)

## Section: Plot favourability
##################################################


p <- plot_fav(aq.sf, FALSE)
p


#produce variable importance plot
varImpPlot(model) 

ext.raster <- rast(x = '01_DATA/INPUT/extremadura_10.tiff')  %>% 
  scale(center=TRUE, scale=TRUE) # present

names(ext.raster)  <- c(paste0("Bio0", 1:9), paste0("Bio", 10:19))

plot(predict(ext.raster, model))

plot(Fav(pred=predict(ext.raster, model), obs=xy$sp))


new <- data.frame(Solar.R=150, Wind=8, Temp=70, Month=5, Day=5)

#use fitted bagged model to predict Ozone value of new observation
predict(newdata=tiff_files_2021_2040_scaled[[1]], model)
