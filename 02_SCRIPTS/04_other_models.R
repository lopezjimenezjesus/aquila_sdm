
#####################################################
## Section: Explore other models
#####################################################



## Section: 
##################################################

library(glmnet)

smod_aa_biohist <- (glm(sp ~ Bio01 + Bio02,
                        family = binomial,  
                        data = xy))


lasso <- glmnet(x = xy %>% dplyr::select("Bio01", "Bio02"),
                y = xy$sp,
                family = "binomial",
                alpha = 1)
lasso_pred <- predict(lasso, xy$sp, type = "response")


## Section:
##################################################
library(dismo)



xy <- right_join(x %>% dplyr::select("CUADRICULA", starts_with("Bio")), y)

xy.sf <- inner_join(utm10, xy, by="CUADRICULA")

sp <- xy.sf %>% dplyr::filter(AQUADA==1) %>% st_coordinates()

xy <- xy %>% select(-"CUADRICULA") %>% rename("sp"="AQUADA")

bc <- bioclim(xy %>% dplyr::select(-"sp") %>% as.matrix(), xy$sp)
