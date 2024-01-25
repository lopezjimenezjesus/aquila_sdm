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


here::i_am('Aquila_adalberti_SDM.Rproj')

## Load data
##################################################

utm10 <- st_read('01_DATA/INPUT/shp/CUTM10_AQUADA.shp') %>% dplyr::select("CUADRICULA")

x <- read_excel(path = here::here('01_DATA/INPUT/CUTM10extVariables.xlsx'), sheet = "CUTM10extVariables")

y <- read_excel(path = here::here('01_DATA/INPUT/Presencia_Imperial_CUTM10.xlsx'), sheet = "Presencia")

xy <- right_join(x,y, by="CUADRICULA")


paste(names(xy)[1:dim(xy)[2]], collapse = " + ")

metadata <- readr::read_delim(file = "01_DATA/OUTPUT/metadata_predictors.csv")


paisaje <- read_csv(file = "01_DATA/OUTPUT/paisaje_predictors.csv")

xy <- xy %>% dplyr::select(all_of(c(paisaje$predictors_names, "AQUADA")))   %>% 
   dplyr::select(-c("La", "Lo", "Lo2", "La2", "LaLoVaria"))  %>% 
   dplyr::select(-c("CazaMe", "CazaMa")) %>%  # correlated with other predictors
   dplyr::select(-c( "TABSMAX1", "TABSMAX7", "TABSMIN1", "TABSMIN7")) %>% 
   dplyr::select(-c( "RAINMAX1", "RAINMAX7","RAINDAY1","RAINDAY7"))# correlated with other predictors
xy

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
## Section: Model not scaled
##################################################

## Collinearity: VIF
##################################################

xy.nosf <- xy %>% st_drop_geometry()


#https://quantifyinghealth.com/vif-threshold/

mult1 <- multicol(subset(xy, select = 1:(ncol(xy)-1)))  # especificar las columnas con variables!
mult1  # muestra la multicolinealidad


VIF_threshold <- config$configuration$vif_threshold

xy.vif <- subset(mult1, VIF < VIF_threshold)

xy.subset <- xy %>% dplyr::select(rownames(xy.vif), AQUADA) 

## FDR - Correlation
##################################################

xy.fdr <- FDR(data = xy.subset, sp.cols = length(xy.subset), var.cols = 1:(length(xy.subset)-1), q = 0.05)

xy.fdr$exclude

variables_to_model <- paste(rownames(xy.fdr$select)[1:length(rownames(xy.fdr$select))], collapse = " + ")

## Fit model
##################################################


# null.model <- glm(AQUADA ~ 1, data=xy, family = binomial)
# 
# model.glm <- step(null.model, direction='forward',
#                   keep =  function(model, aic) list(model = model, aic = aic),
#                   scope = paste0("~", variables_to_model))
# 
# summary(model.glm)

variables_to_model_vector <- str_trim(as.vector(str_split(string = variables_to_model, pattern = "\\+", simplify = TRUE)))

xy_selection <- xy %>% dplyr::select(all_of(variables_to_model_vector), "AQUADA")

# foo <- glm(AQUADA ~ Ciervo + QUESUR + Cul_het + DenCap18 + QFAGPY + NumPoblaD + LongCarrD + EUCSPP + Arroz + DistPobla + LongElectD + CASSAT + Jabali + Conejo, data=foox, family = binomial)
# 
# 
# fooo <- modelTrim(foo)


model.glm  <- fuzzySim::stepwise(xy_selection, length(xy_selection), 1:(length(xy_selection) -1), family = binomial(link="logit"), 
                           simplif=FALSE, direction="forward", trace=2,preds=TRUE, Favourability=FALSE,Wald=TRUE)



## Trimmed model
##################################################

model.glm.trimmed <- modelTrim(model.glm$model) # not necessary but ok to avoid further modification of script due to changes on modelling strategy (see above)

saveRDS(model.glm.trimmed, file= file.path(paths$model_path, 'glm_model_trimmed.rds'))
saveRDS(model.glm.trimmed$predictions, file=file.path(paths$model_path, "glm_model_predictions.rds"))


# Names of our variables:

pred <- names(coefficients(model.glm$model)[2:length(coefficients(model.glm$model))])


## save metadata

saveRDS(model.glm, file=file.path(paths$model_path, "glm_model.rds"))
saveRDS(model.glm$predictions, file=file.path(paths$model_path, "glm_model_predictions.rds"))


##################################################
## Section: Favourability
##################################################

aq.sf <- fav_to_spatial_grid(model.glm.trimmed, xy, x, utm10)

aq10cat.sf <- fav_to_spatial_grid_10cat(model.glm.trimmed, xy, x, utm10)

## Section: Plot favourability
##################################################


p <- plot_fav(aq.sf, FALSE)
p

ggsave(filename = file.path(paths$figure_path, 'fav_historical_model_c3.png'), 
       plot = p, dpi = "retina")


p10 <- plot_fav_10cat(aq10cat.sf, FALSE)
p10
ggsave(filename = file.path(paths$figure_path, 'fav_historical_model_c10.png'), 
       plot = p10, dpi = "retina")


library(patchwork)

q <- (p10 + p) 

dpi = 96

ggsave(filename =  file.path(paths$figure_path, 'fav_historical_model_c3_10.png'), 
         plot = q, width = 1000 / dpi, height = 500 / dpi,
       dpi = dpi)


## Plot intermediate models 
##################################################

# fav_list <- fav_step_models(model.glm) # use model not trimmed for exploring intermediate steps

# fav_list <- fav_step_models(model.glm)
fav_list <- fav_step_models(model.glm$predictions, obs=xy_selection$AQUADA) # use model not trimmed for exploring intermediate steps


# fav_list <- lapply(X = as.list(model.glm$predictions),FUN = function(x) Fav(obs=xy_selection$AQUADA, pred=x ))


fav_plots <- plot_fav_step_models(fav_list)

fav_plots <-  do.call("grid.arrange", c(fav_plots, ncol=3))

ggsave(filename = file.path(paths$figure_path, 'intermediate_models_fav.png'), 
       plot = fav_plots, dpi = "retina")

ggsave(filename =  file.path(paths$figure_path, 'intermediate_models_fav.svg'), 
       plot = fav_plots, dpi = "retina", width = 210, height = 297, units = "mm")

## Variance partition
##################################################

coef_model <-dput(names(coef(model.glm.trimmed)))[-1]

#coef_model <- dput(names(coef(model.glm)))[-1]

predictors_table <- readr::read_delim(file = "01_DATA/OUTPUT/metadata_predictors.csv")


# get groups from predictor_table

groups_vars <-  subset(predictors_table, predictors_names %in% coef_model)

var_groups <- data.frame(variables=groups_vars$predictors_names, factores=groups_vars$factores)

# # get groups from predictor_table
# var_groups <- data.frame(vars=names(coef(model.glm.trimmed))[-1],
#                          groups=predictors_table[predictors_table$predictors_names %in% dput(names(coef(model.glm.trimmed)))[-1], "factores"])
# 

# with(predictors_table, factores[names(coef(model.glm.trimmed))[-1] %in% predictors_names])

png(filename=file.path(paths$figure_path, "variance_plot.png"), width = 800, height = 800)
varPart(model = model.glm.trimmed, groups = var_groups,plot.unexpl=FALSE, cex.names	=1, pred.type = "P")
dev.off()

png(filename=file.path(paths$figure_path, "variance_plot_color.png"), width = 800, height = 800)
varPart(model = model.glm.trimmed, groups = var_groups, pred.type = "P", colr = TRUE,plot.unexpl=FALSE, cex.names	=1)
dev.off()

# add sign of coefficients

var_groups <- var_groups %>% mutate(
  signo = if_else(do.call(
    what = function(x) {x>0}, list(coef(model.glm.trimmed)[2:length(coef(model.glm.trimmed))])),
                "+", "-")) %>% 
  mutate(var_sign=paste0(variables, "(", signo, ")"))

# format to report

var_groups_w <- pivot_wider(data = var_groups, names_from =  factores, values_from = var_sign) 

B <- paste(unique(na.omit(var_groups_w$`Factores Ambientales Bióticos`)))
A <- paste(unique(na.omit(var_groups_w$`Factores Ambientales Abióticos`)))
An <-  paste(unique(na.omit(var_groups_w$`Factores Antrópicos`)))

BAAn <- data.frame(`Factores ambientales Bióticos`= toString(B),
          `Factores ambientales Abióticos`= toString(A),
          `Factores ambientales Antrópicos` = toString(An))


saveRDS(BAAn, file = file.path(paths$model_path, "BAAn.rds"))
saveRDS(var_groups, file= file.path(paths$model_path, "model_variables_factor.rds"))

#

knitr::kable(BAAn, col.names = c("Factores ambientales Bióticos",
                                 "Factores ambientales Abióticos",
                                 "Factores ambientales Antrópicos")) %>%
  kable_styling(full_width = F)


##################################################
## Section: Coefficients, matrix and metrics  ----
##################################################  

## Coefficients ----
##################################################

# tidy_coef <- round(summaryWald(model.glm.trimmed),3) # this is the proper way...but keep function below..lazy

tidy_coef <- build_coef_table(model.glm.trimmed)

write_rds(tidy_coef, file = file.path(paths$model_path, "tabla_coef.rds"))

## Confussion matrix  ----
##################################################

# just temp... need to calculate model using well know function because of bug on fuzzy package (stepwise)

null.model_2 <- glm(AQUADA ~ 1, data=xy, family = binomial)

model.glm_2 <- step(null.model_2, direction='forward',
                  keep =  function(model, aic) list(model = model, aic = aic),
                  scope = paste0("~", variables_to_model))



model.glm.trimmed_2 <- modelTrim(model.glm_2 )


crosspred_glm <- mecofun::crossvalSDM(model=model.glm.trimmed_2, traindat= as.data.frame(xy),
                                      colname_pred=names(coef(model.glm.trimmed))[2:length(coef(model.glm.trimmed))] , 
                                      colname_species = "AQUADA", kfold= 10)

evalsdm <-  mecofun::evalSDM(xy$AQUADA, crosspred_glm)


#tm <- threshMeasures(model = model.glm.trimmed, ylim = c(0, 1), thresh = 0.5)

# using fav instead model

tm <- threshMeasures(pred = Fav((model.glm.trimmed_2)), obs=xy$AQUADA, ylim = c(0, 1), thresh = 0.5)



c_m <-  as.data.frame(tm$ConfusionMatrix)

rownames(c_m) <- c("Favorable", "Desfavorable")
colnames(c_m) <- c("Presencia", "Ausencia")
c_m

saveRDS(c_m, file =  file.path(paths$model_path, "confusion_matrix.rds"))

## Metrics
##################################################

blr <- blr_test_hosmer_lemeshow(model.glm.trimmed)

model.metrics <- data.frame(AUC=evalsdm[c("AUC")], 
                            UPR=tm$ThreshMeasures[c("UPR"),],
                            OPR=tm$ThreshMeasures[c("OPR"),],
                            HyL=blr$pvalue)
model.metrics

saveRDS(object = model.metrics, file =  file.path(paths$model_path, "model_metrics.rds"))

