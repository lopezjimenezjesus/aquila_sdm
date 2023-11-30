

## Functions
##################################################


# filter paisajes variables 

filter_predictors_by_type_analysis <- function(type_analisys="A_Paisaje", name_of_columns=c("tipologia","predictors_names")) {
  
  predictors_table[!is.na(predictors_table[type_analisys]),c(name_of_columns)] 
}


# 
# fav_step_models <- function(model) {
#   # Save intermediate glm models in a list
#   n_models <- dim(model[["keep"]])[2]
#   l <- list()
#   R2 <- list()
#   for(i in seq(1:n_models)) {
#     l[[i]] <- Fav(model[["keep"]][["model", i]])
#     R2[[i]] <-  with(summary(model[["keep"]][["model", i]]), 1 - deviance/null.deviance)
#   }
#   
#   l_R2 <- list(l, R2)
#   
#   return(l_R2)
# }

fav_step_models <- function(model) {
  # Save intermediate glm models in a list
  n_models <- dim(model[["keep"]])[2]
  l <- list()
  R2 <- list()
  for(i in seq(1:n_models)) {
    l[[i]] <- Fav(model[["keep"]][["model", i]])
    R2[[i]] <- summary(lm(Fav(model) ~ Fav(model$keep[,i]$model)))$adj.r.squared

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
  fav <- Fav(model)
  df <- data.frame(CUADRICULA=grid_code_predictor[["CUADRICULA"]], occ=sp_df[["AQUADA"]], f=fav)
  aq.sf <- inner_join(utm_grid, df, by="CUADRICULA")
  aq.sf$category <- cut(aq.sf$f, breaks=c(-Inf, 0.2, 0.8, Inf), labels=c("low","middle","high"))
  return(aq.sf)
}

fav_to_spatial_grid_10cat <- function(model, sp_df, grid_code_predictor, utm_grid) {
  # build a spatial grid based on favourability values (ten classes)
  fav <- Fav(model)
  df <- data.frame(CUADRICULA=grid_code_predictor[["CUADRICULA"]], occ=sp_df[["AQUADA"]], f=fav)
  aq.sf <- inner_join(utm_grid, df, by="CUADRICULA")
  aq.sf$category <- cut(aq.sf$f, breaks=c(0.0, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  return(aq.sf)
}


plot_fav <- function(sf_dataframe, include_presence_points=TRUE) {
  p <- ggplot(data=sf_dataframe) +
    geom_sf(aes(fill=category)) +
    labs(fill="Favorabilidad") + 
    scale_fill_manual(values = c("low" = "red", "middle" = "yellow", "high" = "green"),
                      labels = c("0.00 - 0.20", "0.20 - 0.80", "0.80 - 1.00")) +
    theme_void() +
    theme(legend.position = c(.9, .15))  +
    theme(legend.key.size = unit(0.75, 'cm'))
  
  if(include_presence_points) {
    p <- p + geom_sf(data=st_centroid(sf_dataframe) %>% dplyr::filter(occ==1))
  }
  return(p)
}


plot_fav_10cat <- function(sf_dataframe, include_presence_points=TRUE) {
  p <- ggplot(data=sf_dataframe) +
    geom_sf(aes(fill=f)) +
    labs(fill="Favorabilidad") + 
    scale_fill_steps(low ="#DCF5E9", high = "#226633", 
                     breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                     labels=c("0.00 - 0.10", "0.10 - 0.20", "0.20 - 0.30", 
                              "0.30 - 0.40",  "0.40 - 0.50", "0.50 - 0.60",
                              "0.60 - 0.70", "0.70 - 0.80", "0.80 - 0.90","0.90 - 1.00"),
                     limits = c(0,1)) +
    theme_void() +
    theme(legend.position = c(.9, .145)) +
    theme(legend.key.size = unit(0.7, 'cm'))
  
  
  if(include_presence_points) {
    p <- p + geom_sf(data=st_centroid(sf_dataframe) %>% dplyr::filter(occ==1))
  }
  return(p)
}


load_list_of_tif <- function (path, period) {
  x <- list.files(path = here::here(paste0(path_to_tif, period)), 
                  pattern = "*.tif", full.names = TRUE)
  return(x)
}



# create_folder <- function(path1="04_RESULTS/01_models", path2= "04_RESULTS/02_figures", folder_name="") {
#   path_for_models <- file.path((path1), folder_name)
#   path_for_figures <- file.path((path2), folder_name)
#   
#   if(dir.exists(path_for_models) || dir.exists(path_for_figures)) {
#     print("Folder exists, using it")
#     return(list(model_path=path_for_models, figure_path=path_for_figures))
#   } else {
#     dir.create(file.path((path1), folder_name))
#     dir.create(file.path((path2), folder_name))
#   }
#   
#   return(list(model_path=path_for_models, figure_path=path_for_figures))
#   
# }

create_folder <- function(path="04_RESULTS", folder_name="") {
  
  path_for_models <- file.path(path, folder_name, "01_Models")
  path_for_figures <- file.path(path,folder_name, "02_Figures")
  path_for_report <- file.path(path, folder_name, "03_Report")
  
  if(dir.exists(path_for_models) || dir.exists(path_for_figures) ||  dir.exists(path_for_report)) {
    print("Folder exists, using it")
    return(list(model_path=path_for_models, figure_path=path_for_figures, report_path=path_for_report))
  } else {
    dir.create(path_for_models, recursive =TRUE)
    dir.create(path_for_figures, recursive =TRUE)
    dir.create(path_for_report, recursive =TRUE)
  }
  return(list(model_path=path_for_models, figure_path=path_for_figures))
}

get_wald <- function(model) {
  wald_coef <- sapply(1:length(coef(model)), function(x){
    temp <- aod::wald.test(b=coef(model), Sigma=vcov(model), Terms = x)
    temp$result$chi2[1]
  })
  return(wald_coef)
  
}

build_coef_table <- function(model) {
  tidy_coef <- broom::tidy(model) %>%
    mutate(odd.ratio= exp(coef(model)), wald_coef=get_wald(model))
  
  colnames(tidy_coef) <- c("variables", "β", "ET", "statistic", "Sig.", "Exp(B)", "Wald")
  
  tidy_coef <- tidy_coef[-1,] %>% dplyr::select(c("variables", "β", "ET","Wald", "Sig.", "Exp(B)"))
  return(tidy_coef)
}




