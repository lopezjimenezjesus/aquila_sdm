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

here::i_am('Aquila_adalberti_SDM.Rproj')

## Load data
##################################################

metadata_predictors <- readr::read_delim('01_DATA/OUTPUT/metadata_predictors.csv') %>%
  filter(!is.na(C_NestSite)) %>% select(predictors_names)

xy <- read_excel(path = here::here('01_DATA/INPUT/VARIABLES_NestSite_AQUADA.xlsx'), sheet = "Variables")

xy <- xy %>% 
  mutate(presence= if_else(grepl('^Nestsite', ID), 1, 0)) %>% select(-ID)


##################################################
## Section: configure
##################################################

library(yaml)

config = yaml.load_file("config.yml")


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

variables_to_model <- paste(rownames(xy.fdr$select)[1:length(rownames(xy.fdr$select))], collapse = " + ")


## Fit model
##################################################


null.model <- glm(presence ~ 1, data=xy, family = binomial)

model.glm <- step(null.model, direction='forward',
                  keep =  function(model, aic) list(model = model, aic = aic),
                  scope = paste0("~", variables_to_model))

summary(model.glm)

model.glm <- modelTrim(model.glm)

fav <- data.frame(fav=Fav(model.glm))
colnames(fav) <- "fav"



ggplot(data = fav, mapping = aes(fav))+
  geom_histogram(bins = 10, col="white") + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  # scale_x_continuous(labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4","0.4-0.5", "0.5-0-6", "0.6-0.7", "0.7-0.8","0.8-0.9", "0.9-1")) +
  theme_minimal()


fav$interval <- cut(fav$fav,10)

fav$interval_2 <- cut(fav$fav,breaks=c(0,0.2,0.8,1))



xy$fav <- fav$fav
xy$interval <- as.factor(fav$interval)
xy$interval_2 <- as.factor(fav$interval_2)
xy.summarize <- xy %>% filter(presence==1) %>% group_by(presence, interval) %>% summarise(n=n())
xy.summarize_2 <- xy %>% filter(presence==1) %>% group_by(presence, interval_2) %>% summarise(n=n())

# ggplot(data = fav, mapping = aes(x = interval))+
#   geom_bar( fill="#70AD47") + 
#   labs(x="Favorabilidad", y="Nº de localizaciones") +
#   scale_x_discrete(labels = c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4","0.4 -0.5 ", "0.5 - 0-6", "0.6 - 0.7", "0.7 - 0.8","0.8 - 0.9", "0.9 - 1")) +
# 
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 


p <- ggplot(data = xy.summarize , mapping = aes(x = interval, y=n))+
  geom_bar( fill="#70AD47", stat = "identity") + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  scale_x_discrete(labels = c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4","0.4 -0.5 ", "0.5 - 0-6", "0.6 - 0.7", "0.7 - 0.8","0.8 - 0.9", "0.9 - 1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 


ggsave(filename = file.path(paths$figure_path, 'nest_barplot_10.png'), 
       plot = p, dpi = "retina")


q <- ggplot(data = xy.summarize_2 , mapping = aes(x = interval_2, y=n))+
  geom_col(aes(fill=interval_2)) + 
  labs(x="Favorabilidad", y="Nº de localizaciones") +
  scale_x_discrete(labels = c("0 - 0.2", "0.2 - 0.8", "0.8 - 1")) +
  scale_fill_manual( values = c("red", "yellow",  "#70AD47")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") 

ggsave(filename = file.path(paths$figure_path, 'nest_barplot_3.png'), 
       plot = q, dpi = "retina")

