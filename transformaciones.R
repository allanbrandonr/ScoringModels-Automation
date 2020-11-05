################################################################################
################################################################################

############ preliminares ######################################################

# opciones
options(scipen = 999)

# librerias
library(magrittr)
library(dplyr)
library(ltm)

# cargar las funciones
# source("tratamiento_bases.R", encoding = "UTF-8")


################################################################################
################################################################################

############ funciones para filtros y transformaciones #########################


#### filtro de variables explicativas ##########################################

f_filtro <- function(bd_filtro, 
                     noms_xs = "x", nom_y = "y", nom_b = "base",
                     vesps = list("nohay" = -9999999990),
                     tipo_xs = list("x" = "numerica"), tipo_y = "numerica",
                     solo_train = F, prop_vesp_max = 1, cor_min = 0, cor_max = 1, cv_min = 0){
  
  cat(paste0("----------------------------------------",
             "\n", "----------------------------------------",
             "\n", 
             "\n", "Filtros iniciales",
             "\n"))
  
  cat(paste0("\n", "----------------------------------------",
             "\n", 
             "\n", "Filtro de variables explicativas",
             "\n",
             "\n"))
  
  # base
  bd_filtro %<>% filter(!get(nom_y) %in% vesps[[nom_y]]) %>%
                 data.frame(stringsAsFactors = F)
  
  # base train
  bd_filtro_eval <- bd_filtro %>% 
                    select(nom_b, nom_y, noms_xs) %>% 
                    data.frame(stringsAsFactors = F)
  
  if (solo_train){ bd_filtro_eval %<>% filter(get(nom_b) == "train") }
  
  # variables para filtrar
  list_xs_elim <- lapply(noms_xs, function(nom_x){
    
    cat(paste0(" > ", nom_x,
               "\n"))
    
    bd_filtro_eval_nesp <- bd_filtro_eval %>% filter(!get(nom_x) %in% vesps[[nom_x]])
    
    p_vesp <- NA
    med <- NA
    desv <- NA
    cv <- NA
    cor <- NA
    eliminar <- 0
    
    if (nrow(bd_filtro_eval_nesp) == 0){ 
      
      eliminar <- 1
      p_vesp <- 1
      
    } else { 
      
      # porcentaje de valores especiales
      p_vesp <- 1 - nrow(bd_filtro_eval_nesp)/nrow(bd_filtro_eval)
      p_vesp %<>% round(5)
      
      if (p_vesp > prop_vesp_max){ eliminar <- 1 }
      
      # numericas
      if (tipo_xs[[nom_x]] == "numerica"){
    
        # media
        med <- mean(bd_filtro_eval_nesp[[nom_x]]) %>% round(5)
        
        # desviacion
        desv <- sd(bd_filtro_eval_nesp[[nom_x]]) %>% round(5)
        
        if (desv == 0){ eliminar <- 1 }
        
        # cv
        cv <- ifelse(med == 0, yes = Inf, no = desv/med) %>% round(5)
        
        if (abs(cv) < cv_min){ eliminar <- 1 }
        
        # correlacion
        if (tipo_y == "numerica"){
          
          cor <- ifelse(desv == 0 | is.na(desv), 
                        yes = 0,
                        no = cor(bd_filtro_eval_nesp[[nom_x]], bd_filtro_eval_nesp[[nom_y]])) %>% 
                 round(5)
          
        } else if (tipo_y == "dicotomica"){
          
          cor <- ifelse(desv == 0 | is.na(desv), 
                        yes = 0,
                        no = biserial.cor(bd_filtro_eval_nesp[[nom_x]], bd_filtro_eval_nesp[[nom_y]], level = 2)) %>% 
                 round(5)
          
        }
        
        if (abs(cor) > cor_max | abs(cor) < cor_min){ eliminar <- 1 }  
        
      }
      
    }
    
    # juntar
    x_elim <- data.frame("variable" = nom_x, "p_vesp" = p_vesp, "med" = med, 
                         "desv" = desv, "cv" = cv, "cor" = cor, "eliminar" = eliminar,
                         stringsAsFactors = F)
    
    return(x_elim)
    
  })
  
  # juntar en una base
  xs_elim <- do.call(rbind, list_xs_elim, quote = F)
  
  # base filtrada
  bd_filtro <- bd_filtro[, names(bd_filtro)[!names(bd_filtro) %in% xs_elim$variable[xs_elim$eliminar == 1]]]
  
  filtros <- list("bd" = bd_filtro, 
                  "vars_elim" = xs_elim,
                  "vals_esp" = vesps,
                  "tipos_var" = tipo_xs)

  cat(paste0("\n", "----------------------------------------",
             "\n", "----------------------------------------",
             "\n",
             "\n"))
  
  return(filtros)
  
}


#### transformaciones de variables explicativas  ###############################

f_transformaciones <- function(bd_transf, 
                               noms_xs = "x", nom_y = "y", nom_b = "base",
                               vesps = list("nohay" = -999999999),
                               tipo_xs = list("x" = "numerica"), tipo_y = "numerica",
                               transf_estab_xs = F, solo_train = F,
                               num_lags = 1, cor_min = 0, n_top_cor = Inf){

  
  cat(paste0("----------------------------------------",
             "\n", "----------------------------------------",
             "\n", 
             "\n", "Calculo y eleccion de transformaciones",
             "\n"))
  
  
  # transformaciones de variables explicativas ---------------------------------
  
  if (transf_estab_xs){
    
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Transformaciones de variables explicativas",
               "\n",
               "\n"))
    
    # calculo de transformaciones estabilizadoras por variable explicativas
    list_xs_trans <- lapply(noms_xs, function(nom_x){
      
      # nom_x <- noms_xs[20]
      
      cat(paste0(" > ", nom_x,
                 "\n"))
      
      orig_x <- bd_transf[[nom_x]]
      
      transfs_x <- data.frame(orig_x)
      names(transfs_x) <- "orig_x"
      
      # logaritmo
      if (sum(orig_x[!orig_x %in% vesps[[nom_x]]] <= 0) == 0){
        
        log_x <- ifelse(!orig_x %in% vesps[[nom_x]],
                        yes = log(orig_x),
                        no = orig_x)
        
        # agregar
        transfs_x <- data.frame(transfs_x, log_x)
        names(transfs_x) <- c("orig_x", paste0(nom_x, "_log"))
        
      }
      
      transfs_x %<>% select(-orig_x)
      
      return(transfs_x)
      
    })
    
    # juntar en una base
    bd_xs_trans <- do.call(cbind, list_xs_trans, quote = F)
    noms_xs %<>% c(names(bd_xs_trans))
    
    # actualizar valores especiales
    bd_xs_trans %<>% f_calificar_esp_bd()
    vesps %<>% c(bd_xs_trans$vals_esp)
    bd_xs_trans <- bd_xs_trans$bd
    
    # agregar a base original
    bd_transf %<>% cbind(bd_xs_trans)
    
  }
  
  
  # lags de variables explicativas ---------------------------------------------
  
  cat(paste0("\n", "----------------------------------------",
             "\n", 
             "\n", "Lags de variables explicativas",
             "\n",
             "\n"))
  
  # calculo de mejores lags 
  list_lags_xs_opt <- lapply(noms_xs, function(nom_x){
    
    # nom_x <- noms_xs[1]
    
    i_p <- which(noms_xs == nom_x)[1]
    
    cat(paste0(" > ", nom_x, " (", round((i_p - 1)/length(noms_xs)*100, 2) ,"%)",
               "\n"))
    
    # lags y filtros por correlacion minima deseada
    list_lags_x <- lapply(num_lags, function(l){
      
      # l <- 1
      
      # lags
      orig_x <- bd_transf[[nom_x]]
      lag_x <- lag(orig_x, n = l)  
      lag_dif_x <- lag_x/orig_x - 1
      
      # agregar
      lags_x <- data.frame(bd_transf[[nom_b]], bd_transf[[nom_y]], lag_x, lag_dif_x, stringsAsFactors = F)
      names(lags_x) <- c(nom_b, nom_y, paste0(nom_x, "_lag_", l), paste0(nom_x, "_lag_dif_", l))
      
      if (solo_train){ lags_x_eval <- lags_x %>% filter(!get(nom_y) %in% vesps[[nom_y]] & get(nom_b) == "train") 
      } else { lags_x_eval <- lags_x }
      
      desv_lag_x <- sd(lags_x_eval[[paste0(nom_x, "_lag_", l)]], na.rm = T)
      desv_lag_dif_x <- sd(lags_x_eval[[paste0(nom_x, "_lag_dif_", l)]], na.rm = T)
      
      # correlacion de lag
      if (desv_lag_x == 0 | is.na(desv_lag_x)){
        
        cors_lag_x <- 0
        
      } else if (tipo_y == "numerica"){
        
        cors_lag_x <- abs(cor(lags_x_eval[[paste0(nom_x, "_lag_", l)]], 
                              lags_x_eval[[nom_y]]))  
        
      } else if (tipo_y == "dicotomica"){
        
        cors_lag_x <- abs(biserial.cor(lags_x_eval[[paste0(nom_x, "_lag_", l)]], 
                                       lags_x_eval[[nom_y]], level = 2))  
        
      }
      
      # correlacion de dif
      if (desv_lag_dif_x == 0 | is.na(desv_lag_dif_x)){
        
        cors_lag_dif_x <- 0
        
      } else if (tipo_y == "numerica"){
        
        cors_lag_dif_x <- abs(cor(lags_x_eval[[paste0(nom_x, "_lag_dif_", l)]], 
                                  lags_x_eval[[nom_y]]))  
        
      } else if (tipo_y == "dicotomica"){
        
        cors_lag_dif_x <- abs(biserial.cor(lags_x_eval[[paste0(nom_x, "_lag_dif_", l)]], 
                                           lags_x_eval[[nom_y]], level = 2))  
        
      }
      
      # checar correlacion minima deseada
      cors_lags_x <- rbind(data.frame("var" = paste0(nom_x, "_lag_", l), "cor" = cors_lag_x, stringsAsFactors = F), 
                           data.frame("var" = paste0(nom_x, "_lag_dif_", l), "cor" = cors_lag_dif_x, stringsAsFactors = F)) %>%
                     filter(cor > cor_min)

      # mantener variables 
      lags_x %<>% select(cors_lags_x$var)
      
      lags_x <- list("bd_lags_x" = lags_x, "bd_cors_lags_x" = cors_lags_x)
      
      return(lags_x)
      
    })
    
    # juntar en una base
    bd_lags_x_opt <- data.frame("nom_x" = bd_transf[[nom_x]])
    bd_cors_lags_x_opt <- data.frame()
    
    for (i in 1:length(list_lags_x)){
      
      bd_lags_x_opt %<>% cbind(list_lags_x[[i]][["bd_lags_x"]])
      bd_cors_lags_x_opt %<>% rbind(list_lags_x[[i]][["bd_cors_lags_x"]])
      
    }
    
    # elejir lags con mejores correlaciones 
    bd_cors_lags_x_opt %<>% arrange(desc(cor)) %>% filter(row_number() <= n_top_cor)
    bd_lags_x_opt %<>% select(bd_cors_lags_x_opt$var)
    
    lags_x_opt <- list("bd_lags_x" = bd_lags_x_opt, "bd_cors_lags_x" = bd_cors_lags_x_opt)
    
    return(lags_x_opt)
    
  })
  
  # juntar en una base
  bd_lags_xs_opt <- data.frame("nom_y" = bd_transf[[nom_y]])
  bd_cors_lags_xs_opt <- data.frame()
  
  for (i in 1:length(list_lags_xs_opt)){
    
    bd_lags_xs_opt %<>% cbind(list_lags_xs_opt[[i]][["bd_lags_x"]])
    bd_cors_lags_xs_opt %<>% rbind(list_lags_xs_opt[[i]][["bd_cors_lags_x"]])
    
  }
  
  # correlaciones
  bd_cors_lags_xs_opt %<>% arrange(desc(cor)) %>% mutate(cor = round(cor, 5))
  
  bd_lags_xs_opt %<>% select(-nom_y)
  
  # actualizar valores especiales
  bd_lags_xs_opt %<>% f_calificar_esp_bd()
  vesps %<>% c(bd_lags_xs_opt$vals_esp)
  bd_lags_xs_opt <- bd_lags_xs_opt$bd
  
  # agregar a base original
  bd_transf %<>% cbind(bd_lags_xs_opt)
  
  transfs_info <- f_calificar_esp_bd(bd_transf)
  
  transfs <- list("bd" = bd_transf, 
                  "correlaciones" = bd_cors_lags_xs_opt,
                  "vals_esp" = transfs_info$vals_esp,
                  "tipos_var" = transfs_info$tipos_var)
  
  cat(paste0("\n", "----------------------------------------",
             "\n", "----------------------------------------",
             "\n",
             "\n"))
  
  return(transfs)
  
}
