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
library(blorr)
library(gridExtra)
library(broom)
library(flextable)

here::i_am('Aquila_adalberti_SDM.Rproj')

# x <- st_read(dsn = here::here('01_DATA/INPUT/shp/CUTM10_AQUADA.shp')) %>% select('CUADRICULA', 'AQUADA') %>% st_drop_geometry()
# head(x)


## Load data
##################################################

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")


xy <- right_join(x,y, by="CUADRICULA")


# filter variables

xy <- xy %>% dplyr::select(-c("CUADRICULA", "Ori_W", "Ori_S",
                              "Area_ha", "Area_en_Ex", "Percent_Ex",
                              "La", "Lo", "Lo2", "La2", "LaLoVaria",
                              "CazaMe", "CazaMa", 
                              "DenArena", "DenCaliza",  "DenGranito", "DenPizarra",
                              "DenEmbalse", "DenHidro",
                              "S", "SE", "SW", "W", "N", "N", "NW", "NE", "E"))

paste(names(xy)[1:dim(xy)[2]], collapse = " + ")


## Functions
##################################################




fav_step_models <- function(model) {
  # Save intermediate models in a list
  n_models <- dim(model[["keep"]])[2]
  l <- list()
  R2 <- list()
  for(i in seq(1:n_models)) {
    l[[i]] <- Fav(model[["keep"]][["model", i]])
    R2[[i]] <-  with(summary(model[["keep"]][["model", i]]), 1 - deviance/null.deviance)
  }
  
  l_R2 <- list(l, R2)
  
  return(l_R2)
}

plot_fav_step_models <- function(fav_R2_list) {
  # Save favourability ggplot from intermediate models in a list
  n_fav <-length(fav_R2_list[[1]])
  p <- list()
  for(i in seq(1:n_fav)) {
    fav <- fav_R2_list[[1]][[i]]
    aq <- data.frame(CUADRICULA=x$CUADRICULA, sp=xy$AQUADA, f=fav)
    
    aq.sf <- inner_join(utm10, aq, by="CUADRICULA")
    
    aq.sf$category <- cut(aq.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
    
    p[[i]] <- ggplot(data=aq.sf) +
      geom_sf(aes(fill=category)) +
      scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
      # geom_sf(data=st_centroid(aq.sf) %>% dplyr::filter(sp==1)) +
      annotate(geom = "text", x = 316696.7, y = 4240000, 
               label = paste("Paso: ", i, "\n", "R2=", round(fav_R2_list[[2]][[i]],3)),
               color = "black", size = 2) +
      guides(fill="none") +
      theme_void()
    
  }
  return(p)
}


fav_to_spatial_grid <- function(model, sp_df, grid_code_predictor, utm_grid) {
  # build a spatial grid based on favourability values
  fav <- Fav(model.glm)
  df <- data.frame(CUADRICULA=grid_code_predictor[["CUADRICULA"]], occ=sp_df[["AQUADA"]], f=fav)
  aq.sf <- inner_join(utm_grid, df, by="CUADRICULA")
  aq.sf$category <- cut(aq.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  return(aq.sf)
}




##################################################
## Section: Model not scaled
##################################################

## FDR - Correlation
##################################################

xy.fdr <- FDR(data = xy, sp.cols = length(xy), var.cols = 1:(length(xy)-1), q = 0.1)

xy.fdr$exclude

paste(rownames(xy.fdr$select)[1:length(rownames(xy.fdr$select))], collapse = " + ")

## Fit
##################################################

null.model <- glm(AQUADA ~ 1, data=xy, family = binomial)

model.glm <- step(null.model, direction='forward',
                  keep =  function(model, aic) list(model = model, aic = aic),
                  scope = (~Cul_len + TABSMAX1 + Ciervo + PSpr + QUESUR + Reg + 
                             Cul_het + Sup_arti + RAINDAY1 + DenCap18 + QFAGPY + 
                             P_DIAS + QUEILE + Ptot + P_prim + NumPoblaD + Deh +
                             LongCarrD + Psum + DenPobla + EUCSPP + Pwin + Arroz + 
                             SinVeg + DistPobla + Pene + Pjul + LongElectD + CASSAT +
                             Jabali + Pvar + Tmax7 + Paut + RAINDAY7 + Conejo +
                             RAINMAX7 + TABSMIN1 + Tsum))


saveRDS(model.glm, file="04_RESULTS/01_models/glm_model.rds")

## Favourability
##################################################

aq.sf <- fav_to_spatial_grid(model.glm, xy, x, utm10)


## Section: Plot favourability
##################################################

plot_fav <- function(sf_dataframe, include_presence_points=TRUE) {
  p <- ggplot(data=sf_dataframe) +
    geom_sf(aes(fill=category)) +
    scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green")) +
    theme_void()

  if(include_presence_points) {
    p <- p + geom_sf(data=st_centroid(sf_dataframe) %>% dplyr::filter(occ==1))
  }
  return(p)
} 


p <- plot_fav(aq.sf, FALSE)
p

ggsave(filename = '04_RESULTS/02_figures/fav_historical_model.png', 
       plot = p, dpi = "retina")


## Plot intermediate models 
##################################################


fav_list <- fav_step_models(model.glm)
fav_plots <- plot_fav_step_models(fav_list)

fav_plots <-  do.call("grid.arrange", c(fav_plots, ncol=3))

ggsave(filename = '04_RESULTS/02_figures/intermediate_models_fav.png', 
       plot = fav_plots, dpi = "retina", width = 210, height = 297, units = "mm")



## Trimmed model
##################################################

model.glm.trimmed <- modelTrim(model.glm)

saveRDS(model.glm.trimmed, file="04_RESULTS/01_models/glm_model_trimmed.rds")


## Variance partition
##################################################

dput(colnames(model.glm.trimmed$data))


var_groups <- data.frame(vars =c("PSpr", "TABSMAX1", "QUESUR", "RAINDAY1", "Ciervo", 
                                 "Paut", "DenCap18", "Reg", "Conejo", "Pene", "Pvar", "LongElectD"),
                         groups = c("Factores ambientales Abióticos", "Factores ambientales Abióticos",
                                    "Factores ambientales Bióticos",   "Factores ambientales Abióticos", 
                                    "Factores ambientales Bióticos",  "Factores ambientales Abióticos", 
                                    "Factores ambientales Antrópicos", "Factores ambientales Abióticos",
                                    "Factores ambientales Bióticos",  "Factores ambientales Abióticos", 
                                    "Factores ambientales Abióticos", "Factores ambientales Antrópicos"))

png(filename="04_RESULTS/02_figures/variance_plot.png", width = 800, height = 800)
varPart(model = model.glm.trimmed, groups = var_groups,plot.unexpl=FALSE, cex.names	=1)
dev.off()

png(filename="04_RESULTS/02_figures/variance_plot_color.png", width = 800, height = 800)
varPart(model = model.glm.trimmed, groups = var_groups, pred.type = "P", colr = TRUE,plot.unexpl=FALSE, cex.names	=1)
dev.off()


##################################################
## Section: Coefficients, metrics, etc
##################################################  


## Coefficients
##################################################

tidy_coef <- broom::tidy(model.glm.trimmed)

odd.ratio <- exp(coef(model.glm.trimmed))


foo <- cbind(tidy_coef, odd.ratio)


library(kableExtra)
tidy_coef.kable <-  knitr::kable(tidy_coef, "html") %>%
  kable_styling("striped") %>% 
  row_spec(0, background = "#F7CAAC") %>%
  row_spec(seq(1,dim(tidy_coef)[1],2), background = "#FBE4D5")

tidy_coef.flex <- flextable(tidy_coef)

saveRDS(tidy_coef.flex, file = "01_DATA/OUTPUT/foo_flex.rds")

saveRDS(tidy_coef.kable, file = "01_DATA/OUTPUT/foo_kable.rds")

save_kable(tidy_coef.kable, file = "01_DATA/OUTPUT/foo_kable.png")

coeficients_df <- tidy(model.glm.trimmed) # coefficients


glance(model.glm.trimmed) #  summary statistics

head(augment(model.glm)) #  fitted values and residuals

## de Hosmer y Lemeshow 
##################################################


blr <- blr_test_hosmer_lemeshow(model.glm.trimmed)
blr

## Wald 
##################################################

summary(model.glm.trimmed, statistic = "Wald")

w <- waldtest(object = model.glm)
drop1(model.glm, test = "LRT") 


library(aod)

foo <- aod::wald.test(Sigma = vcov(model.glm.trimmed), b = coef(model.glm.trimmed), Terms = 1)

foo$result$chi2


## Confussion matrix and metrics
##################################################

library(mecofun)


crosspred_glm <- mecofun::crossvalSDM(model.glm.trimmed, traindat= xy,
                                      colname_pred=colnames(xy)[1:(dim(xy)[2]-1)], 
                                      colname_species = "AQUADA", kfold= 10)

mecofun::evalSDM(xy$AQUADA, crosspred_glm)

confusionMatrix(obs = xy$AQUADA, pred =  Fav(model.glm), thresh = 0.5)

tm <- threshMeasures(model = model.glm, ylim = c(0, 1), thresh = 0.5)
tm <- threshMeasures(model = model.glm, ylim = c(0, 1), thresh = "preval")

cm <- tm$ConfusionMatrix

cm <- as.data.frame(cm)

rownames(cm) <- c("Favorable", "Desfavorable")
colnames(cm) <- c("Presencia", "Ausencia")
cm
