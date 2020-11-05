################################################################################
################################################################################

############ preliminares ######################################################

# opciones
options(scipen = 999)

# librerias
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(weights)
library(coin)
library(pROC)
library(ltm)


# Posibles valores
# tipo_x: "numerica", "nominal", "ordinal"
# tipo_y: "numerica", "dicotomica"
# nom_tests: "welch", "anova", "wald", "wilcoxon", "mood"
# nom_medidas_biv: "correlacion", "correlacion_x", "ecm", "dam", "gini", "auc"
# nom_info_calif: "bucket", "ym", "woe" y cualquier columna de info_biv
# nom_graficas_biv: "bivariado", "dispersion", "series", "densidades_clasif", "boxplot_clasif"


################################################################################
################################################################################

############ funciones para estimacion de buckets ##############################


#### estimacion de buckets #####################################################

f_analisis_buckets <- function(bd_x, tipo_analisis = "train", info_biv = NULL,
                               nom_x = "x", nom_y = "y", nom_p = "peso", nom_r = "ref", nom_id = "id",
                               vesp_x = -9999999990, 
                               tipo_x = "numerica", tipo_y = "numerica", 
                               nom_tests = c("welch", "anova"), nom_medidas_biv = "correlacion", 
                               nom_info_calif = c("bucket", "ym"), nom_graficas_biv = c("bivariado", "dispersion"),
                               pobl_prop_min = .05, signi = .1){
  
  # nombres
  if (nom_r == nom_p){
    
    bd_x[[paste0(nom_r, "_r")]] <- bd_x[[nom_r]]
    nom_r <- paste0(nom_r, "_r")
    
  }
  
  # mantener variables
  bd_x %<>% select(nom_x, nom_p, nom_y, nom_r, nom_id) %>% data.frame(stringsAsFactors = F)
  
  
  # analisis en test -----------------------------------------------------------------------
  
  if (tipo_analisis == "test" & !is.null(info_biv)){
    
    tendencia <- info_biv$tendencia[info_biv$tipo_valor == "no_especial"][1]
    
    if(tipo_x == "numerica"){
      
      info_biv %<>% mutate(desde = as.numeric(desde), 
                           hasta = as.numeric(hasta))
      
    }
    
    # calificar base
    bd_biv <- f_califica_buckets(bd_x, info_biv, nom_x, nom_id, tipo_x, unique(c(nom_info_calif, "ym", "bucket", "tipo_valor"))) %>%
              data.frame(stringsAsFactors = F)
    
    # re calcular info de tabla de bivariados para corregir errores de redondeo
    info_biv <- left_join(info_biv %>% select(desde, hasta, bucket),
                          f_info_biv(bd_biv, nom_x, nom_y, nom_p, nom_r, nom_id, paste0(nom_x, "_bucket"), paste0(nom_x, "_tipo_valor"), tipo_x, tipo_y, tendencia, nom_tests),
                          by = c("bucket" = "bucket")) %>%
                data.frame(stringsAsFactors = F)
    
    # medidas
    medidas_biv <- f_medidas_biv(bd_biv, nom_x, nom_y, paste0(nom_x, "_ym"), vesp_x, tipo_x, tipo_y, nom_medidas_biv)
    
    # ordenar variables
    info_biv <- cbind(info_biv %>% select(variable, tendencia, tipo_variable, tipo_valor, desde, hasta),
                      info_biv %>% select(-variable, -tendencia, -tipo_variable, -tipo_valor, -desde, -hasta))
    
    
    # analisis en train -----------------------------------------------------------------------
    
  } else if (tipo_analisis == "train"){
  
    
    # trameado valores especiales -----------------------------------------------------------------------
    
    # base de valores especiales
    bd_esp <- bd_x %>% 
              filter(get(nom_x) %in% vesp_x) %>% 
              rename(bucket = nom_x) %>% 
              data.frame(stringsAsFactors = F)
    
    info_biv_esp <- data.frame()
    
    if (nrow(bd_esp) > 0){
      
      # tabla de bivariados
      info_biv_esp <- bd_esp %>% 
                      mutate(vyp = get(nom_y)*get(nom_p)) %>%
                      group_by(bucket) %>%
                      summarise(ym = sum(vyp)/sum(get(nom_p))) %>%
                      data.frame(stringsAsFactors = F) %>%
                      mutate(desde = bucket, hasta = bucket, tipo_valor = "especial") %>%
                      data.frame(stringsAsFactors = F)
      
      # ordenar variables
      info_biv_esp %<>% select(desde, hasta, ym, bucket, tipo_valor)
      
    }
    
    
    # trameado valores no especiales -----------------------------------------------------------------------
    
    # numero de observaciones totales y de observaciones minimas para garantizar porcentaje
    n_tot <- sum(!bd_x[[nom_x]] %in% vesp_x)
    n_pct <- max(c(ceiling(n_tot*pobl_prop_min), 1))
    
    
    # numericas u ordinarias -----------------------------------------------------------------------
  
    if (tipo_x %in% c("numerica", "ordinal")){
      
      
      # tendencia decreciente -----------------------------------------------------------------------
      
      # poner en mismo orden
      bd_nesp_dec <- bd_x %>% 
                     filter(!get(nom_x) %in% vesp_x) %>% 
                     mutate(vyp = get(nom_y)*get(nom_p)) %>% 
                     arrange(get(nom_x), vyp) %>% 
                     select(-vyp) %>% 
                     data.frame(stringsAsFactors = F)
      
      # vectores de observaciones pendientes de tramear
      vp_post_dec <- bd_nesp_dec[[nom_p]]
      vx_post_dec <- bd_nesp_dec[[nom_x]]
      vy_post_dec <- bd_nesp_dec[[nom_y]]
      ym_post_dec <- cumsum(vy_post_dec*vp_post_dec)/cumsum(vp_post_dec)
      ym_post_dec[is.na(ym_post_dec)] <- -Inf
      ym_prev_dec <- Inf
      
      # numero de ultima observacion ya trameada y numero de observaciones pendientes de tramear
      acum_cut_dec <- 0 
      n_post_dec <- n_tot
      
      # construccion de buckets mientras haya observaciones pendientes de tramear
      info_biv_dec <- data.frame()
      
      while (acum_cut_dec < n_tot) {
        
        # si numero de observaciones pendientes es menor o igual al minimo para garantizar porcentaje
        if (n_post_dec <= n_pct) {
          
          ym_post_dec_optim <- ym_post_dec[n_post_dec]
          cut_iter_dec <- n_post_dec
          
          # limites y ym
          desde_iter <- ifelse(tipo_x == "numerica", yes = -Inf, no = "-Inf")
          hasta_iter <- vx_post_dec[cut_iter_dec]
          ym_iter <- sum(vy_post_dec[1:cut_iter_dec]*vp_post_dec[1:cut_iter_dec])/sum(vp_post_dec[1:cut_iter_dec])
          
          # acumular info
          info_biv_dec <- rbind(info_biv_dec,
                                data.frame("desde" = desde_iter, "hasta" = hasta_iter, "ym" = ym_iter, stringsAsFactors = F))
          
          # si numero de observaciones pendientes es mayor al minimo para garantizar porcentaje
        } else {
          
          # calcular maximo
          tend_mal_dec <- TRUE
          while (tend_mal_dec){
            
            # definir numero de observacion donde de limite superior de bucket actual
            ym_post_dec_optim <- max(ym_post_dec[n_pct:n_post_dec])
            cut_iter_dec <- max(which(ym_post_dec == ym_post_dec_optim))
            
            # no separar observaciones con mismo valor
            while ((cut_iter_dec < n_post_dec) & (vx_post_dec[cut_iter_dec] == vx_post_dec[cut_iter_dec + 1])) {
              
              ym_post_dec[1:(max(which(vx_post_dec == vx_post_dec[cut_iter_dec])) - 1)] <- -Inf
              ym_post_dec_optim <- max(ym_post_dec[n_pct:n_post_dec])
              cut_iter_dec <- max(which(ym_post_dec == ym_post_dec_optim))
              
            }
            
            # si bucket no tiene tendencia adecuada por errores de redondeo
            if (sum(vy_post_dec[1:cut_iter_dec]*vp_post_dec[1:cut_iter_dec])/sum(vp_post_dec[1:cut_iter_dec]) > ym_prev_dec){
              
              ym_post_dec[cut_iter_dec] <- -Inf
              
              # terminar si bucket tiene tendencia adecuada considerando errores de redondeo
            } else {
              
              tend_mal_dec <- FALSE
              
            }
            
          }
          
          # limites y ym
          desde_iter <- ifelse(tipo_x == "numerica", yes = -Inf, no = "-Inf")
          hasta_iter <- vx_post_dec[cut_iter_dec]
          ym_iter <- sum(vy_post_dec[1:cut_iter_dec]*vp_post_dec[1:cut_iter_dec])/sum(vp_post_dec[1:cut_iter_dec])
          
          # acumular info
          info_biv_dec <- rbind(info_biv_dec,
                                data.frame("desde" = desde_iter, "hasta" = hasta_iter, "ym" = ym_iter, stringsAsFactors = F))
                                 
          # media de y del bucket
          ym_prev_dec <- info_biv_dec$ym[nrow(info_biv_dec)]
          
          # actualizar 
          vp_post_dec <- vp_post_dec[(cut_iter_dec + 1):n_post_dec]
          vx_post_dec <- vx_post_dec[(cut_iter_dec + 1):n_post_dec]
          vy_post_dec <- vy_post_dec[(cut_iter_dec + 1):n_post_dec]
          ym_post_dec <- cumsum(vy_post_dec*vp_post_dec)/cumsum(vp_post_dec)
          ym_post_dec[is.na(ym_post_dec)] <- -Inf
          n_post_dec <- length(vy_post_dec)
          
        }
        
        # actualizar numero de observaciones pendientes
        acum_cut_dec <- acum_cut_dec + cut_iter_dec
        
      }
      
      # variables de limites de intervalos
      n_buckets_dec <- nrow(info_biv_dec)
      info_biv_dec$desde[-1] <- info_biv_dec$hasta[-n_buckets_dec]
      info_biv_dec$bucket <- ifelse(tipo_x == "numerica",
                                    yes = paste0("(", round(info_biv_dec$desde, 5), ", ", round(info_biv_dec$hasta, 5), "]"),
                                    no = paste0("(", info_biv_dec$desde, ", ", info_biv_dec$hasta, "]"))
      info_biv_dec$hasta[n_buckets_dec] <- ifelse(tipo_x == "numerica", yes = Inf, no = "Inf")
      info_biv_dec$bucket[n_buckets_dec] <- ifelse(tipo_x == "numerica",
                                                   yes = paste0("(", round(info_biv_dec$desde[n_buckets_dec], 5), ", ", round(info_biv_dec$hasta[n_buckets_dec], 5), ")"),
                                                   no = paste0("(", info_biv_dec$desde[n_buckets_dec], ", ", info_biv_dec$hasta[n_buckets_dec], ")"))
      
      # agregar valores especiales
      info_biv_dec %<>% mutate(tipo_valor = "no_especial") %>% rbind(info_biv_esp)
      
      if (tipo_y == "dicotomica"){ info_biv_dec %<>% mutate(woe = log((1 - ym)/ym))}
      
      # calificar base
      bd_dec <- f_califica_buckets(bd_calif = bd_x, info_biv_calif = info_biv_dec, nom_x, nom_id, tipo_x, nom_calif = unique(c(nom_info_calif, "ym", "bucket", "tipo_valor"))) %>%
                data.frame(stringsAsFactors = F)
      
      # re calcular info de tabla de bivariados para corregir errores de redondeo
      info_biv_dec <- left_join(info_biv_dec %>% select(desde, hasta, bucket),
                                f_info_biv(bd_buckets = bd_dec, nom_x, nom_y, nom_p, nom_r, nom_id, nom_bu = paste0(nom_x, "_bucket"), nom_tv = paste0(nom_x, "_tipo_valor"), tipo_x, tipo_y, tendencia = "decreciente", nom_tests),
                                by = c("bucket" = "bucket")) %>%
                      data.frame(stringsAsFactors = F)
      
      # medidas
      medidas_dec <- f_medidas_biv(bd_medidas = bd_dec, nom_x, nom_y, nom_ym = paste0(nom_x, "_ym"), vesp_x, tipo_x, tipo_y, nom_medidas_biv)
      
      
      # tendencia creciente -----------------------------------------------------------------------
      
      # poner en orden opuesto
      bd_nesp_cre <- bd_x %>% 
                     filter(!get(nom_x) %in% vesp_x) %>% 
                     mutate(vyp = get(nom_y)*get(nom_p)) %>% 
                     arrange(get(nom_x), desc(vyp)) %>% 
                     select(-vyp) %>% 
                     data.frame(stringsAsFactors = F)
      
      # vectores de observaciones pendientes de tramear
      vp_post_cre <- bd_nesp_cre[[nom_p]]
      vx_post_cre <- bd_nesp_cre[[nom_x]]
      vy_post_cre <- bd_nesp_cre[[nom_y]]
      ym_post_cre <- cumsum(vy_post_cre*vp_post_cre)/cumsum(vp_post_cre)
      ym_post_cre[is.na(ym_post_cre)] <- Inf
      ym_prev_cre <- -Inf
      
      # numero de ultima observacion ya trameada y numero de observaciones pendientes de tramear
      acum_cut_cre <- 0 
      n_post_cre <- n_tot
      
      # construccion de buckets mientras haya observaciones pendientes de tramear
      info_biv_cre <- data.frame()
      
      while (acum_cut_cre < n_tot) {
        
        # si numero de observaciones pendientes es menor o igual al minimo para garantizar porcentaje
        if (n_post_cre <= n_pct) {
          
          ym_post_cre_optim <- ym_post_cre[n_post_cre]
          cut_iter_cre <- n_post_cre
          
          # limites y ym
          desde_iter <- ifelse(tipo_x == "numerica", yes = -Inf, no = "-Inf")
          hasta_iter <- vx_post_cre[cut_iter_cre]
          ym_iter <- sum(vy_post_cre[1:cut_iter_cre]*vp_post_cre[1:cut_iter_cre])/sum(vp_post_cre[1:cut_iter_cre])
          
          # acumular info
          info_biv_cre <- rbind(info_biv_cre,
                                data.frame("desde" = desde_iter, "hasta" = hasta_iter, "ym" = ym_iter, stringsAsFactors = F))
                                 
          # si numero de observaciones pendientes es mayor al minimo para garantizar porcentaje
        } else {
          
          # calcular minimo
          tend_mal_cre <- TRUE
          while (tend_mal_cre){
            
            # definir numero de observacion donde de limite superior de bucket actual
            ym_post_cre_optim <- min(ym_post_cre[n_pct:n_post_cre])
            cut_iter_cre <- max(which(ym_post_cre == ym_post_cre_optim))
            
            # no separar observaciones con mismo valor
            while ((cut_iter_cre < n_post_cre) & (vx_post_cre[cut_iter_cre] == vx_post_cre[cut_iter_cre + 1])) {
              
              ym_post_cre[1:(max(which(vx_post_cre == vx_post_cre[cut_iter_cre])) - 1)] <- Inf
              ym_post_cre_optim <- min(ym_post_cre[n_pct:n_post_cre])
              cut_iter_cre <- max(which(ym_post_cre == ym_post_cre_optim))
              
            }
    
            # si bucket no tiene tendencia adecuada por errores de redondeo
            if (sum(vy_post_cre[1:cut_iter_cre]*vp_post_cre[1:cut_iter_cre])/sum(vp_post_cre[1:cut_iter_cre]) < ym_prev_cre){
              
              ym_post_cre[cut_iter_cre] <- Inf
              
              # terminar si bucket tiene tendencia adecuada considerando errores de redondeo
            } else {
              
              tend_mal_cre <- FALSE
              
            }
            
          }
          
          # limites y ym
          desde_iter <- ifelse(tipo_x == "numerica", yes = -Inf, no = "-Inf")
          hasta_iter <- vx_post_cre[cut_iter_cre]
          ym_iter <- sum(vy_post_cre[1:cut_iter_cre]*vp_post_cre[1:cut_iter_cre])/sum(vp_post_cre[1:cut_iter_cre])
          
          # acumular info
          info_biv_cre <- rbind(info_biv_cre,
                                data.frame("desde" = desde_iter, "hasta" = hasta_iter, "ym" = ym_iter, stringsAsFactors = F))
          
          # media de y del bucket
          ym_prev_cre <- info_biv_cre$ym[nrow(info_biv_cre)]
          
          # actualizar 
          vp_post_cre <- vp_post_cre[(cut_iter_cre + 1):n_post_cre]
          vx_post_cre <- vx_post_cre[(cut_iter_cre + 1):n_post_cre]
          vy_post_cre <- vy_post_cre[(cut_iter_cre + 1):n_post_cre]
          ym_post_cre <- cumsum(vy_post_cre*vp_post_cre)/cumsum(vp_post_cre)
          ym_post_cre[is.na(ym_post_cre)] <- Inf
          n_post_cre <- length(vy_post_cre)
          
        }
        
        # actualizar numero de observaciones pendientes
        acum_cut_cre <- acum_cut_cre + cut_iter_cre
        
      }
      
      # variables de limites de intervalos
      n_buckets_cre <- nrow(info_biv_cre)
      info_biv_cre$desde[-1] <- info_biv_cre$hasta[-n_buckets_cre]
      info_biv_cre$bucket <- ifelse(tipo_x == "numerica",
                                    yes = paste0("(", round(info_biv_cre$desde, 5), ", ", round(info_biv_cre$hasta, 5), "]"),
                                    no = paste0("(", info_biv_cre$desde, ", ", info_biv_cre$hasta, "]"))
      info_biv_cre$hasta[n_buckets_cre] <- ifelse(tipo_x == "numerica", yes = Inf, no = "Inf")
      info_biv_cre$bucket[n_buckets_cre] <- ifelse(tipo_x == "numerica",
                                                   yes = paste0("(", round(info_biv_cre$desde[n_buckets_cre], 5), ", ", round(info_biv_cre$hasta[n_buckets_cre], 5), ")"),
                                                   no = paste0("(", info_biv_cre$desde[n_buckets_cre], ", ", info_biv_cre$hasta[n_buckets_cre], ")"))
      
      # agregar valores especiales
      info_biv_cre %<>% mutate(tipo_valor = "no_especial") %>% rbind(info_biv_esp)
      
      if (tipo_y == "dicotomica"){ info_biv_cre %<>% mutate(woe = log((1 - ym)/ym))}
      
      # calificar base
      bd_cre <- f_califica_buckets(bd_x, info_biv_cre, nom_x, nom_id, tipo_x, unique(c(nom_info_calif, "ym", "bucket", "tipo_valor"))) %>%
                data.frame(stringsAsFactors = F)
      
      # re calcular info de tabla de bivariados para corregir errores de redondeo
      info_biv_cre <- left_join(info_biv_cre %>% select(desde, hasta, bucket),
                                f_info_biv(bd_cre, nom_x, nom_y, nom_p, nom_r, nom_id, paste0(nom_x, "_bucket"), paste0(nom_x, "_tipo_valor"), tipo_x, tipo_y, "creciente", nom_tests),
                                by = c("bucket" = "bucket")) %>%
                      data.frame(stringsAsFactors = F)
      
      # medidas
      medidas_cre <- f_medidas_biv(bd_cre, nom_x, nom_y, paste0(nom_x, "_ym"), vesp_x, tipo_x, tipo_y, nom_medidas_biv)
  
      
      # tendencia optima -----------------------------------------------------------------------
      
      # mantener tendencia con mejor medida
      if (medidas_cre[[nom_medidas_biv[1]]] >= medidas_dec[[nom_medidas_biv[1]]]){
        
        tendencia <- "creciente"
        info_biv <- info_biv_cre
        bd_biv <- bd_cre
        medidas_biv <- medidas_cre
        
      } else {
        
        tendencia <- "decreciente"
        info_biv <- info_biv_dec
        bd_biv <- bd_dec
        medidas_biv <- medidas_dec
        
      }
      
      
      # nominales -----------------------------------------------------------------------
      
    } else if (tipo_x == "nominal"){
      
      tendencia <- "decreciente"
      
      # base para hacer trameado
      bd_biv_nesp <- bd_x %>% 
                     filter(!get(nom_x) %in% vesp_x) %>% 
                     rename(bucket = nom_x) %>% 
                     data.frame(stringsAsFactors = F)
      
      # tabla de bivariados
      info_biv <- bd_biv_nesp %>% 
                  mutate(vyp = get(nom_y)*get(nom_p)) %>%
                  group_by(bucket) %>%
                  summarise(ym = sum(vyp)/sum(get(nom_p))) %>%
                  data.frame(stringsAsFactors = F) %>%
                  mutate(desde = bucket, hasta = bucket, tipo_valor = "no_especial") %>%
                  select(desde, hasta, ym, bucket, tipo_valor) %>%
                  rbind(info_biv_esp)
      
      if (tipo_y == "dicotomica"){ info_biv %<>% mutate(woe = log((1 - ym)/ym))}
      
      # calificar base
      bd_biv <- f_califica_buckets(bd_x, info_biv, nom_x, nom_id, tipo_x, unique(c(nom_info_calif, "ym", "bucket", "tipo_valor"))) %>%
                data.frame(stringsAsFactors = F)
      
      # re calcular info de tabla de bivariados para corregir errores de redondeo
      info_biv <- left_join(info_biv %>% select(desde, hasta, bucket),
                            f_info_biv(bd_biv, nom_x, nom_y, nom_p, nom_r, nom_id, paste0(nom_x, "_bucket"), paste0(nom_x, "_tipo_valor"), tipo_x, tipo_y, "decreciente", nom_tests),
                            by = c("bucket" = "bucket")) %>%
                  data.frame(stringsAsFactors = F)
      
      # medidas
      medidas_biv <- f_medidas_biv(bd_biv, nom_x, nom_y, paste0(nom_x, "_ym"), vesp_x, tipo_x, tipo_y, nom_medidas_biv)
      
    }
    
    # ordenar variables
    info_biv <- cbind(info_biv %>% select(variable, tendencia, tipo_variable, tipo_valor, desde, hasta),
                      info_biv %>% select(-variable, -tendencia, -tipo_variable, -tipo_valor, -desde, -hasta))
    
  
    # pegado de buckets valores no especiales ----------------------------------------------------------------------
    
    comparar_buckets <- T
    
    while (comparar_buckets){
      
      # valores por buckets
      info_biv_bucket <- info_biv %>% 
                         filter(tipo_valor == "no_especial") %>% 
                         select(bucket, ym, pobl, pesos, paste0("pval_", nom_tests)) %>% 
                         data.frame(stringsAsFactors = F) %>%
                         group_by(bucket) %>%
                         filter(row_number() == 1) %>%
                         data.frame(stringsAsFactors = F)
      
      # revisar condiciones para pegado de buckets
      idx <- ((tipo_y == "dicotomica" & sum(info_biv_bucket$ym %in% c(0, 1)) > 0) | 
                (info_biv_bucket$pobl <= n_pct) |
                (apply(info_biv_bucket[paste0("pval_", nom_tests)] > signi, 1, sum) > 0))
      
      # pegar buckets
      if (sum(idx) > 0 & nrow(info_biv_bucket) > 1){
        
        # valores por limites y buckets
        info_biv %<>% filter(tipo_valor == "no_especial") %>% 
                      select(desde, hasta, bucket, ym, pobl, pesos, paste0("pval_", nom_tests))
        
        # indice de buckets que deben pegarse
        k <- min(which(idx)[1], nrow(info_biv_bucket) - 1)
        buckets_pegar <- c(k, k + 1)
        
        # y media
        ym_p <- sum(info_biv_bucket$ym[buckets_pegar]*info_biv_bucket$pesos[buckets_pegar])/sum(info_biv_bucket$pesos[buckets_pegar])
        
        
        # numericas u ordinarias  -----------------------------------------------------------------------
        
        if (tipo_x %in% c("numerica", "ordinal")){
          
          # valores de buckets nuevos
          desde_p <- info_biv$desde[buckets_pegar[1]]
          hasta_p <- info_biv$hasta[buckets_pegar[2]]
          bucket_p <- ifelse(tipo_x == "numerica",
                             yes = paste0("(", round(desde_p, 5), ", ", round(hasta_p, 5), ifelse(hasta_p == Inf, yes = ")", "]")),
                             no = paste0("(", desde_p, ", ", hasta_p, ifelse(hasta_p == "Inf", yes = ")", "]")))
          
          # tabla con buckets pegados
          info_biv$hasta[buckets_pegar[1]] <- hasta_p
          info_biv$ym[buckets_pegar[1]] <- ym_p
          info_biv$bucket[buckets_pegar[1]] <- bucket_p
          
          info_biv <- info_biv[-buckets_pegar[2], ] 
          
          
          # nominales  -----------------------------------------------------------------------
          
        } else if (tipo_x == "nominal"){
          
          # valores de buckets nuevos
          bucket_p <- info_biv_bucket$ym[buckets_pegar[1]]
          
          # tabla con buckets pegados
          info_biv$ym[info_biv$bucket %in% info_biv_bucket$bucket[buckets_pegar]] <- ym_p
          info_biv$bucket[info_biv$bucket %in% info_biv_bucket$bucket[buckets_pegar]] <- bucket_p
          
        }
        
        # agregar valores especiales
        info_biv %<>% select(desde, hasta, bucket, ym) %>% mutate(tipo_valor = "no_especial") %>% rbind(info_biv_esp)
        
        if (tipo_y == "dicotomica"){ info_biv %<>% mutate(woe = log((1 - ym)/ym))}
        
        # calificar base
        bd_biv <- f_califica_buckets(bd_x, info_biv, nom_x, nom_id, tipo_x, unique(c(nom_info_calif, "ym", "bucket", "tipo_valor"))) %>%
                  data.frame(stringsAsFactors = F)
        
        # re calcular info de tabla de bivariados para corregir errores de redondeo
        info_biv <- left_join(info_biv %>% select(desde, hasta, bucket),
                              f_info_biv(bd_biv, nom_x, nom_y, nom_p, nom_r, nom_id, paste0(nom_x, "_bucket"), paste0(nom_x, "_tipo_valor"), tipo_x, tipo_y, tendencia, nom_tests),
                              by = c("bucket" = "bucket")) %>%
                    data.frame(stringsAsFactors = F)
        
        # medidas
        medidas_biv <- f_medidas_biv(bd_biv, nom_x, nom_y, paste0(nom_x, "_ym"), vesp_x, tipo_x, tipo_y, nom_medidas_biv)
        
        # ordenar variables
        info_biv <- cbind(info_biv %>% select(variable, tendencia, tipo_variable, tipo_valor, desde, hasta),
                          info_biv %>% select(-variable, -tendencia, -tipo_variable, -tipo_valor, -desde, -hasta))
        
      } else {
        
        comparar_buckets <- F
        
      }
      
    }
  
  }

  
  # graficas  -----------------------------------------------------------------------
  
  graficas_biv <- f_graficas_biv(bd_graf = bd_biv, info_biv_graf = info_biv, nom_x, nom_y, nom_r, tipo_x, tipo_y, nom_graficas_biv)
  
  # mantener varivalores calificados
  bd_biv %<>% select(nom_id, paste0(nom_x, "_", nom_info_calif))
  
  # buckets, base calificada y graficas
  buckets_estim <- list("info_biv" = info_biv,
                        "bd_calif" = bd_biv,
                        "medidas_biv" = medidas_biv,
                        "graficas_biv" = graficas_biv)
    
  return(buckets_estim)
  
}


################################################################################
################################################################################

############ funciones para tablas de bivariados ###############################


#### calculo de tabla de bivariados ############################################

f_info_biv <- function(bd_buckets, 
                       nom_x = "x", nom_y = "y", nom_p = "peso", nom_r = "ref", nom_id = "id", nom_bu = "bucket", nom_tv = "tipo_valor",
                       tipo_x = "numerica", tipo_y = "numerica", 
                       tendencia = "decreciente", nom_tests = c("welch", "anova")){

  # base para hacer trameado
  bd_buckets %<>% select(nom_bu, nom_tv, nom_p, nom_y, nom_r, nom_id) %>% data.frame(stringsAsFactors = F)
  
  # tabla de bivariados
  info_biv <- bd_buckets %>% mutate(vyp = get(nom_y)*get(nom_p)) %>%
                             group_by(get(nom_tv), get(nom_bu)) %>%
                             summarise(y = sum(get(nom_y)),
                                       ym = sum(vyp)/sum(get(nom_p)),
                                       pobl = n(),
                                       ref = sum(get(nom_r)),
                                       pesos = sum(get(nom_p))) %>%
                             data.frame(stringsAsFactors = F) %>%
                             mutate(y_prop = y/sum(y),
                                    pobl_prop = pobl/sum(pobl),
                                    ref_prop = ref/sum(ref),
                                    pesos_prop = pesos/sum(pesos)) %>%
                             data.frame(stringsAsFactors = F)
  
  if (tipo_y == "dicotomica"){ info_biv %<>% mutate(woe = log((1 - ym)/ym))}
  
  names(info_biv)[names(info_biv) == "get.nom_tv."] <- "tipo_valor"
  names(info_biv)[names(info_biv) == "get.nom_bu."] <- "bucket"
  
  # tablas de bivariados de valores especiales y no especiales
  info_biv_esp <- info_biv %>% filter(tipo_valor == "especial")
  info_biv_nesp <- info_biv %>% filter(tipo_valor == "no_especial")
  
  # ordenar por tendencia
  if (tendencia == "decreciente"){ info_biv_nesp %<>% arrange(desc(ym)) %>% data.frame(stringsAsFactors = F)
  } else if (tendencia == "creciente"){ info_biv_nesp %<>% arrange(ym) %>% data.frame(stringsAsFactors = F) }

  # tests estadisticos para buckets de valores no especiales
  list_pvals_tests <- lapply(1:nrow(info_biv_nesp), function(i){
    
    # i <- 1
    
    if (i < nrow(info_biv_nesp)){
    
      # filtrar buckets y ordenar por ym para test
      info_biv_iter <- info_biv_nesp[c(i, i + 1), ]
        
      # filtrar buckets
      bd_buckets_iter <- bd_buckets %>% filter(get(nom_bu) %in% info_biv_iter$bucket)
      
      # tests estadisticos para tramos consecutivos
      pval_test_iter <- f_tests_biv(info_biv_test = info_biv_iter, bd_buckets_test = bd_buckets_iter, nom_y, nom_p, nom_bu, tipo_y, nom_tests)
      
    } else {
      
      pval_test_iter <- data.frame(matrix(rep(-99, length(nom_tests)), nrow = 1))
      names(pval_test_iter) <- paste0("pval_", nom_tests)
      
    }
    
    return(pval_test_iter)
    
  })
  
  biv_pvals_tests <- do.call(rbind, list_pvals_tests)
  
  # pegar valores p
  info_biv_nesp %<>% cbind(biv_pvals_tests)
  
  if (nrow(info_biv_esp) > 0){ info_biv_esp[paste0("pval_", nom_tests)] <- -99 }
  
  # test de diferencia de yms
  info_biv <- data.frame("variable" = nom_x, "tendencia" = tendencia, "tipo_variable" = tipo_x, rbind(info_biv_nesp, info_biv_esp), stringsAsFactors = F)

  return(info_biv)

}


#### calculo de tests de tabla de bivariados ###################################

f_tests_biv <- function(info_biv_test, bd_buckets_test, 
                        nom_y = "y", nom_p = "peso", nom_bu = "bucket",
                        tipo_y = "numerica", 
                        nom_tests = c("welch", "anova")){

  # ordenar ym para pruebas
  info_biv_test %<>% arrange(ym)
  
  # variable categoria para test
  bd_buckets_test %<>% mutate(bucket_fact = factor(get(nom_bu), levels = info_biv_test$bucket),
                              peso_test = get(nom_p))
  
  
  # test Welch diferencia yms ponderadas  -----------------------------------------------------------------------

  if (sum(nom_tests == "welch")){

    # si no hay observaciones para estimar varianza
    if (sum(info_biv_test$pobl == 1) > 0){

      info_biv_test$pval_welch <- 1

      # si hay observaciones para estimar varianza
    } else {

      # test de welch para buckets consecutivos (aternativa mu1 < mu2)
      welch_test <- wtd.t.test(x = bd_buckets_test[[nom_y]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[1]],
                               y = bd_buckets_test[[nom_y]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[2]],
                               weight = bd_buckets_test[[nom_p]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[1]],
                               weighty = bd_buckets_test[[nom_p]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[2]],
                               samedata = F, alternative = "less")

      # valor p de test
      info_biv_test$pval_welch <- welch_test$coefficients["p.value"] %>%
                                   round(5)

    }

  }


  # test ANOVA significanca regresion lineal -----------------------------------------------------------------------

  if (sum(nom_tests == "anova") > 0){
    
    info_biv_test$pval_anova <- -99
    
    if (tipo_y == "numerica"){
      
      # glm de buckets consecutivos
      lm_test <- lm(as.formula(paste0(nom_y, "~bucket_fact")), 
                    data = bd_buckets_test, weights = peso_test) %>%
                 summary()
      
      # valor p de test (aternativa mu1 < mu2)
      info_biv_test$pval_anova <- lm_test$coefficients[-1, "t value"] %>%
                                  pt(df = nrow(bd_buckets_test) - 2, lower.tail = F) %>%
                                  round(5)

    }
    
  }
  
  
  # test Wald significanca regresion logistica -----------------------------------------------------------------------
  
  if (sum(nom_tests == "wald") > 0){
    
    info_biv_test$pval_wald <- -99
    
    if (tipo_y == "dicotomica"){
      
      # glm de buckets consecutivos
      glm_test <- glm(as.formula(paste0(nom_y, "~bucket_fact")), 
                      data = bd_buckets_test, family = binomial, weights = peso_test) %>%
                  summary()
      
      # valor p de test (aternativa mu1 < mu2)
      info_biv_test$pval_wald <- glm_test$coefficients[-1, "z value"] %>%
                                 pnorm(lower.tail = F) %>%
                                 round(5)
    
    } 
    
  }
  
  
  
  # test Mann Whitney Wilcoxon diferencia dsitribuciones -----------------------------------------------------------------------
  
  if (sum(nom_tests == "wilcoxon") > 0){
    
    info_biv_test$pval_wilcoxon <- -99
    
    if (tipo_y == "numerica"){
      
      # glm de buckets consecutivos
      wilcoxon_test <- wilcox.test(x = bd_buckets_test[[nom_y]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[1]],
                                   y = bd_buckets_test[[nom_y]][bd_buckets_test[[nom_bu]] == info_biv_test$bucket[2]], 
                                   paired = F, alternative = "less")  
      
      # valor p de test (aternativa localizacion 1 < localizacion 2)
      info_biv_test$pval_wilcoxon <- wilcoxon_test$p.value %>%
                                     round(5)
      
    } 
    
  }
  
  
  # test Mood diferencia medianas -----------------------------------------------------------------------
  
  if (sum(nom_tests == "mood") > 0){
    
    info_biv_test$pval_mood <- -99
    
    if (tipo_y == "numerica"){
      
      # glm de buckets consecutivos
      mood_test <- median_test(as.formula(paste0(nom_y, "~bucket_fact")), 
                               data = bd_buckets_test, alternative = "less")  
      
      # valor p de test (aternativa mediana 1 < mediana 2)
      info_biv_test$pval_mood <- pvalue(mood_test) %>%
                                 round(5)
      
    } 
    
  }
  
  pvals <- info_biv_test %>% 
           data.frame(stringsAsFactors = F) %>%
           filter(row_number() == 1) %>% 
           select(paste0("pval_", nom_tests))
                   
  return(pvals)

}


################################################################################
################################################################################

############ funciones para calificacion de bd #################################


#### calificacion de bd ########################################################

f_califica_buckets <- function(bd_calif, info_biv_calif, 
                               nom_x = "x", nom_id = "id",
                               tipo_x = "numerica", 
                               nom_calif = c("bucket", "ym")){

  # tablas de bivariados de valores no especiales
  info_biv_calif_nesp <- info_biv_calif %>%
                         filter(tipo_valor == "no_especial") %>%
                         mutate(hasta = ifelse(hasta == "Inf" & tipo_x == "ordinal",
                                               yes = paste(rep("z", 200), collapse = ""),
                                               no = hasta)) %>%
                         select(hasta, nom_calif) %>%
                         data.frame(stringsAsFactors = F)
  
  names(info_biv_calif_nesp)[-1] <- paste0(nom_x, "_", nom_calif)

  # tablas de bivariados de valores especiales
  info_biv_calif_esp <- info_biv_calif %>%
                        filter(tipo_valor == "especial")  %>%
                        select(hasta, nom_calif) %>%
                        data.frame(stringsAsFactors = F)
  
  names(info_biv_calif_esp)[-1] <- paste0(nom_x, "_", nom_calif)

  # base para calificar
  bd_calif <- bd_calif[!names(bd_calif) %in% paste0(nom_x, "_", nom_calif)] %>% 
              data.frame(stringsAsFactors = F) %>% 
              mutate(id_order = row_number())
  
  # base para calificar de valores especiales y no especiales
  bd_calif_nesp <- bd_calif %>% filter(!get(nom_x) %in% info_biv_calif_esp$hasta)
  bd_calif_esp <- bd_calif %>% filter(get(nom_x) %in% info_biv_calif_esp$hasta)


  # calificar valores no especiales -----------------------------------------------------------------------


  # numericas u ordinarias -----------------------------------------------------------------------

  if (tipo_x %in% c("numerica", "ordinal")){

    # calificar base por pedazos
    if (nrow(bd_calif_nesp) > 100000){
      idx_iter <- seq(1, nrow(bd_calif_nesp), by = 100000)
      idx_iter[length(idx_iter)] <- nrow(bd_calif_nesp)
    } else {
      idx_iter <- c(1, nrow(bd_calif_nesp))
    }

    info_biv_calif_nesp %<>% mutate(aux = "aux")

    list_bd_calif_final <- lapply(2:length(idx_iter), function(i){

      # pegar bucket
      bd_nesp_iter <- bd_calif_nesp[idx_iter[i - 1]:idx_iter[i], ] %>%
                      mutate(aux = "aux") %>%
                      left_join(info_biv_calif_nesp, by = c("aux" = "aux")) %>%
                      filter(get(nom_x) <= hasta) %>%
                      group_by(get(nom_id)) %>%
                      filter(hasta == min(hasta))  %>%
                      data.frame(stringsAsFactors = F) %>%
                      select(-hasta, -aux, -"get.nom_id.") %>%
                      data.frame(stringsAsFactors = F)
      
      return(bd_nesp_iter)

    })

    bd_calif_final <- do.call(rbind, list_bd_calif_final, quote = F)


    # nominales -----------------------------------------------------------------------

  } else if (tipo_x == "nominal"){

    vals_join <- "hasta"
    names(vals_join) <- nom_x
      
    bd_calif_final <- left_join(bd_calif_nesp, info_biv_calif_nesp, by = vals_join) %>%
                      data.frame(stringsAsFactors = F)

  }


  # calificar valores especiales -----------------------------------------------------------------------

  if (nrow(bd_calif_esp) > 0){

    vals_join <- "hasta"
    names(vals_join) <- nom_x
    
    bd_calif_esp %<>% left_join(info_biv_calif_esp, by = vals_join) %>%
                      data.frame(stringsAsFactors = F)

    bd_calif_final %<>% rbind(bd_calif_esp)

  }
  
  # ordenar base
  bd_calif_final %<>% arrange(id_order) %>% select(-id_order)

  return(bd_calif_final)

}


################################################################################
################################################################################

############ funciones para calculo de estadisticas ############################


#### calculo de medidas de relacion con y ######################################

f_medidas_biv <- function(bd_medidas, 
                          nom_x = "x", nom_y = "y", nom_ym = "ym",
                          vesp_x = -9999999990,
                          tipo_x = "numerica", tipo_y = "numerica", 
                          nom_medidas_biv = "correlacion"){

  # base para calificar
  bd_medidas %<>% select(nom_x, nom_y, nom_ym) %>% data.frame(stringsAsFactors = F)

  # base para calificar de valores especiales y no especiales
  bd_medidas_nesp <- bd_medidas %>% filter(!get(nom_x) %in% vesp_x)

  medidas <- data.frame("variable" = nom_x, "tipo_variable" = tipo_x)
    
  
  # correlacion escala de x calificada con nom_ym -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "correlacion") > 0){
    
    if (tipo_y == "numerica"){
      
      medidas[["correlacion_nesp"]] <- ifelse(sd(bd_medidas_nesp[[nom_ym]]) == 0, yes = 0, no = abs(cor(bd_medidas_nesp[[nom_ym]], bd_medidas_nesp[[nom_y]])))
      medidas[["correlacion"]] <- ifelse(sd(bd_medidas[[nom_ym]]) == 0, yes = 0, no = abs(cor(bd_medidas[[nom_ym]], bd_medidas[[nom_y]])))
    
    } else if (tipo_y == "dicotomica"){
      
      medidas[["correlacion_nesp"]] <- ifelse(sd(bd_medidas_nesp[[nom_ym]]) == 0, yes = 0, no = abs(biserial.cor(bd_medidas_nesp[[nom_ym]], bd_medidas_nesp[[nom_y]], level = 2)))
      medidas[["correlacion"]] <- ifelse(sd(bd_medidas[[nom_ym]]) == 0, yes = 0, no = abs(biserial.cor(bd_medidas[[nom_ym]], bd_medidas[[nom_y]], level = 2)))
      
    }
    
  }
  
  
  # correlacion escala de x original -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "correlacion_x") > 0){
    
    medidas[["correlacion_x_nesp"]] <- -9999999990
    medidas[["correlacion_x"]] <- -9999999990
    
    if (tipo_x == "numerica"){
      
      if (tipo_y == "numerica"){
        
        medidas[["correlacion_x_nesp"]] <- ifelse(sd(bd_medidas_nesp[[nom_x]]) == 0, yes = 0, no = abs(cor(bd_medidas_nesp[[nom_x]], bd_medidas_nesp[[nom_y]])))
        medidas[["correlacion_x"]] <- ifelse(sd(bd_medidas[[nom_x]]) == 0, yes = 0, no = abs(cor(bd_medidas[[nom_x]], bd_medidas[[nom_y]])))
      
      } else if (tipo_y == "dicotomica"){
        
        medidas[["correlacion_x_nesp"]] <- ifelse(sd(bd_medidas_nesp[[nom_x]]) == 0, yes = 0, no = abs(biserial.cor(bd_medidas_nesp[[nom_x]], bd_medidas_nesp[[nom_y]], level = 2)))
        medidas[["correlacion_x"]] <- ifelse(sd(bd_medidas[[nom_x]]) == 0, yes = 0, no = abs(biserial.cor(bd_medidas[[nom_x]], bd_medidas[[nom_y]], level = 2)))
        
      }
      
    }
    
  }
  
  
  # ecm -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "ecm") > 0){
    
    medidas[["ecm_nesp"]] <- -9999999990
    medidas[["ecm"]] <- -9999999990
    
    if (tipo_y == "numerica"){
      
      medidas[["ecm_nesp"]] <- mean((bd_medidas_nesp[[nom_ym]] - bd_medidas_nesp[[nom_y]])^2)
      medidas[["ecm"]] <- mean((bd_medidas[[nom_ym]] - bd_medidas[[nom_y]])^2)
    
    }
    
  }
  
  
  # dam -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "dam") > 0){
    
    medidas[["dam_nesp"]] <- -9999999990
    medidas[["dam"]] <- -9999999990
    
    if (tipo_y == "numerica"){
      
      medidas[["dam_nesp"]] <- mean(abs(bd_medidas_nesp[[nom_ym]] - bd_medidas_nesp[[nom_y]]))
      medidas[["dam"]] <- mean(abs(bd_medidas[[nom_ym]] - bd_medidas[[nom_y]]))
    
    }
    
  }
  
  
  # gini -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "gini") > 0){
    
    medidas[["gini_nesp"]] <- -9999999990
    medidas[["gini"]] <- -9999999990
    
    if (tipo_y == "dicotomica"){
      
      medidas[["gini_nesp"]] <- auc(roc(bd_medidas_nesp[[nom_y]], bd_medidas_nesp[[nom_ym]]))*2 - 1
      medidas[["gini"]] <- auc(roc(bd_medidas[[nom_y]], bd_medidas[[nom_ym]]))*2 - 1
    
    }
    
  }
  
  
  # auc -----------------------------------------------------------------------
  
  if (sum(nom_medidas_biv == "auc") > 0){
    
    medidas[["auc_nesp"]] <- -9999999990
    medidas[["auc"]] <- -9999999990
    
    if (tipo_y == "dicotomica"){
      
      medidas[["auc_nesp"]] <- auc(roc(bd_medidas_nesp[[nom_y]], bd_medidas_nesp[[nom_ym]]))
      medidas[["auc"]] <- auc(roc(bd_medidas[[nom_y]], bd_medidas[[nom_ym]]))
    
    }
    
  }
  
  return(medidas)

}


################################################################################
################################################################################

############ funciones para graficas ###########################################


#### graficas de relacion con y ################################################

f_graficas_biv <- function(bd_graf, info_biv_graf, 
                           nom_x = "x", nom_y = "y", nom_r = "ref",
                           tipo_x = "numerica", tipo_y = "numerica", 
                           nom_graficas_biv = c("bivariado", "dispersion")){
  
  info_biv_graf <- info_biv_graf[!duplicated(info_biv_graf$bucket),]
  
  # tablas de bivariados de valores especiales y no especiales
  info_biv_graf_nesp <- info_biv_graf %>% filter(tipo_valor == "no_especial") 
  info_biv_graf_esp <- info_biv_graf %>% filter(tipo_valor == "especial")
  
  info_biv_graf <- rbind(info_biv_graf_esp, info_biv_graf_nesp)
  
  # base para calificar de valores especiales y no especiales
  bd_graf_nesp <- bd_graf %>% filter(!get(nom_x) %in% info_biv_graf_esp$hasta)
  bd_graf_esp <- bd_graf %>% filter(get(nom_x) %in% info_biv_graf_esp$hasta)
  
  graficas <- list()
  
  # grafica bivariados -----------------------------------------------------------------------
  
  if (sum(nom_graficas_biv == "bivariado")){
    
    # variable para gráfico
    info_biv_graf$plot_group <- sort(paste0(letters[1:min(26, nrow(info_biv_graf))], 1:nrow(info_biv_graf)))
    
    # grafico de ym y distribucion de observaciones
    rango_ym <- max(info_biv_graf$ym) - min(info_biv_graf$ym)
    min_ym <- max(min(info_biv_graf$ym), 0)
    
    graficas$bivariados_pobl <- ggplot(info_biv_graf) + 
                                geom_bar(aes(x = plot_group, y = pobl_prop), stat = "identity", fill = "#006EC1", group = 1) +
                                geom_point(aes(x = plot_group, y = (ym - min_ym)/rango_ym), size = 1.5, color = "#191970", group = 1) +
                                geom_text(aes(x = plot_group, y = (ym - min_ym)/rango_ym, label = round(ym, 2)), hjust = .5, vjust = -1, size = 2.5, group = 1) +
                                geom_line(data = info_biv_graf %>% filter(tipo_valor == "no_especial"),
                                          aes(x = plot_group, y = (ym - min_ym)/rango_ym), group = 1, size = 1, color = "#191970") +
                                theme_minimal() +
                                scale_y_continuous(limits = c(0, 1.2), sec.axis = sec_axis(~.*rango_ym + min_ym, name = "y media")) +
                                scale_x_discrete(breaks = info_biv_graf$plot_group,
                                                 labels = info_biv_graf$bucket) +
                                labs(title = "Bivariado", subtitle = nom_y, x = nom_x, y = "distribucion de observaciones") +
                                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    # grafico de ym y distribucion de variable de referencia
    graficas$bivariados_ref <- ggplot(info_biv_graf) + 
                               geom_bar(aes(x = plot_group, y = ref_prop), stat = "identity", fill = "#006EC1", group = 1) +
                               geom_point(aes(x = plot_group, y = (ym - min_ym)/rango_ym), size = 1.5, color = "#191970", group = 1) +
                               geom_text(aes(x = plot_group, y = (ym - min_ym)/rango_ym, label = round(ym, 2)), hjust = .5, vjust = -1, size = 2.5, group = 1) +
                               geom_line(data = info_biv_graf %>% filter(tipo_valor == "no_especial"),
                                         aes(x = plot_group, y = (ym - min_ym)/rango_ym), group = 1, size = 1, color = "#191970") +
                               theme_minimal() +
                               scale_y_continuous(limits = c(0, 1.2), sec.axis = sec_axis(~.*rango_ym + min_ym, name = "y media")) +
                               scale_x_discrete(breaks = info_biv_graf$plot_group,
                                                labels = info_biv_graf$bucket) +
                               labs(title = "Bivariado", subtitle = nom_y, x = nom_x, y = paste0("distribucion de ", nom_r)) +
                               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica dispersion -----------------------------------------------------------------------
  
  if (sum(nom_graficas_biv == "dispersion")){
    
    # grafico de dispersion
    graficas$dispersion <- ggplot(bd_graf_nesp) +
                           geom_point(aes(x = get(nom_x), y = get(nom_y)), alpha = .2, color = "#191970", group = 1) +
                           theme_minimal() +
                           labs(title = "Dispersion", subtitle = "valores no especiales", x = "x", y = "y") +
                           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de series -----------------------------------------------------------------------
  
  if (sum(nom_graficas_biv == "series")){
    
    graficas$series <- ggplot(data = bd_graf_nesp %>% mutate(num_obs = row_number())) +
                       geom_point(aes(x = num_obs, y = get(nom_x), color = "x"), alpha = .8, group = 1) +
                       geom_point(aes(x = num_obs, y = get(nom_y), color = "y"), alpha = .8, group = 1) +
                       theme_minimal() +
                       scale_color_manual(name = "tipo", values = c("x" = "#191970", "y" = "#006EC1")) +
                       labs(title = "Series", subtitle = "valores no especiales", x = "numero de observacion", y = "") +
                       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de densidades densidades por valor observado -----------------------------------------------------------------------
  
  if (sum(nom_graficas_biv == "densidades_clasif") & tipo_y == "dicotomica"){
    
    graficas$densidades_clasif <- ggplot() +
                                  geom_density(data = bd_graf_nesp %>% filter(get(nom_y) == 1), aes(get(nom_x), fill = "positivos"), alpha = .8) +
                                  geom_density(data = bd_graf_nesp %>% filter(get(nom_y) == 0), aes(get(nom_x), fill = "negativos"), alpha = .8) +
                                  theme_minimal() +
                                  scale_fill_manual(name = "tipo", values = c("positivos" = "#191970", "negativos" = "#006EC1")) +
                                  labs(title = "Densidades", subtitle = "valores no especiales", x = "y", y = "densidad") +
                                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
                                
  }
  
  
  # grafica de boxplot por valor observado -----------------------------------------------------------------------
  
  if (sum(nom_graficas_biv == "boxplot_clasif")  & tipo_y == "dicotomica"){
    
    graficas$boxplot_clasif <- ggplot(data = bd_graf_nesp) +
                               geom_boxplot(aes(x = factor(get(nom_y)), y = get(nom_x)), color = "#191970", group = 1) +
                               theme_minimal() +
                               labs(title = "Boxplot", subtitle = "valores no especiales", x = "y", y = "distribucion") +
                               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
                            
  }
  
  
  return(graficas)       
  
}


################################################################################
################################################################################

############ funciones para analisis bivariados ################################


#### analisis de bivariados ####################################################

f_analisis_biv <- function(bd_analisis, nom_analisis = c("train", "test"), info_biv_eval = NULL,
                           noms_xs = "x", nom_y = "y", nom_p = "peso", nom_r = "ref", nom_id = "id", nom_b = "base",
                           vesps = list("nohay" = -9999999990),
                           tipo_xs = list("x" = "numerica"), tipo_y = "numerica", 
                           nom_tests = c("welch", "anova"), nom_medidas_biv = "correlacion", 
                           nom_info_calif = c("bucket", "ym"), nom_graficas_biv = c("bivariado", "dispersion"),
                           pobl_prop_min = .05, signi = .1){

  cat(paste0("----------------------------------------",
             "\n", "----------------------------------------",
             "\n", 
             "\n", "Analisis de bivariados",
             "\n"))
  
  # nombres
  if (nom_r == nom_p){
    
    bd_analisis[[paste0(nom_r, "_r")]] <- bd_analisis[[nom_r]]
    nom_r <- paste0(nom_r, "_r")
    
  }
  
  # filtrar base
  bd_analisis %<>% select(nom_y, nom_p, nom_r, nom_id, nom_b, noms_xs) %>% 
                   filter(!get(nom_y) %in% vesps[[nom_y]]) %>%
                   data.frame(stringsAsFactors = F)
  
  # bases train, test
  bd_analisis_train <- bd_analisis %>% filter(get(nom_b) == "train")
  bd_analisis_test <- bd_analisis %>% filter(get(nom_b) == "test")

  analisis <- list()
  bd_calif_f <- data.frame()
  
  # analisis en train -----------------------------------------------------------------------
  
  if (nrow(bd_analisis_train) > 0 & sum(nom_analisis == "train") > 0){
  
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Analisis en train",
               "\n",
               "\n"))

    # esitmacion de buckets
    list_bucekts_estim <- lapply(noms_xs, function(nom_var){
      
      # nom_var <- noms_xs[1]
      
      i_p <- which(noms_xs == nom_var)[1]
      
      cat(paste0(" > ", nom_var, " (", round((i_p - 1)/length(noms_xs)*100, 2) ,"%)",
                 "\n"))
      
      buckets_estim_var <- f_analisis_buckets(bd_x = bd_analisis_train %<>% select(nom_y, nom_p, nom_r, nom_id, nom_var), 
                                              tipo_analisis = "train", 
                                              info_biv = NULL,
                                              nom_x = nom_var, 
                                              nom_y = nom_y, 
                                              nom_p = nom_p, 
                                              nom_r = nom_r, 
                                              nom_id = nom_id,
                                              vesp_x = vesps[[nom_var]], 
                                              tipo_x = tipo_xs[[nom_var]], 
                                              tipo_y = tipo_y, 
                                              nom_tests = nom_tests, 
                                              nom_medidas_biv = nom_medidas_biv, 
                                              nom_info_calif = nom_info_calif, 
                                              nom_graficas_biv = nom_graficas_biv,
                                              pobl_prop_min = pobl_prop_min, 
                                              signi = signi)

      return(buckets_estim_var)
      
    })
  
    # juntar en una base
    info_biv_train <- data.frame()
    medidas_biv_train <- data.frame()
    graficas_biv_train <- list()
    
    for (i in 1:length(list_bucekts_estim)){
      
      bd_analisis_train %<>% left_join(list_bucekts_estim[[i]][["bd_calif"]], by = nom_id)
      
      info_biv_train %<>% rbind(list_bucekts_estim[[i]][["info_biv"]])
      medidas_biv_train %<>% rbind(list_bucekts_estim[[i]][["medidas_biv"]])
      graficas_biv_train %<>% c(list_bucekts_estim[[i]][["graficas_biv"]])
      
    }
    
    # info para calificar test
    info_biv_eval <- info_biv_train
    
    # acumular
    bd_calif_f %<>% rbind(bd_analisis_train)
    analisis[["info_biv"]][["train"]] <- info_biv_train
    analisis[["medidas_biv"]][["train"]] <- medidas_biv_train
    analisis[["graficas_biv"]][["train"]] <- graficas_biv_train
    
  }
  
  
  # analisis en test -----------------------------------------------------------------------
  
  if (nrow(bd_analisis_test) > 0 & sum(nom_analisis == "test") > 0 & !is.null(info_biv_eval)){
    
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Analisis en test",
               "\n",
               "\n"))
    
    # esitmacion de buckets
    list_bucekts_eval <- lapply(noms_xs, function(nom_var){
      
      # nom_var <- noms_xs[4]
      
      i_p <- which(noms_xs == nom_var)[1]
      
      cat(paste0(" > ", nom_var, " (", round((i_p - 1)/length(noms_xs)*100, 2) ,"%)",
                 "\n"))
      
      buckets_eval_var <- f_analisis_buckets(bd_x = bd_analisis_test %<>% select(nom_y, nom_p, nom_r, nom_id, nom_var), 
                                             tipo_analisis = "test", 
                                             info_biv = info_biv_eval %>% filter(variable == nom_var),
                                             nom_x = nom_var, 
                                             nom_y = nom_y, 
                                             nom_p = nom_p, 
                                             nom_r = nom_r, 
                                             nom_id = nom_id,
                                             vesp_x = vesps[[nom_var]], 
                                             tipo_x = tipo_xs[[nom_var]], 
                                             tipo_y = tipo_y, 
                                             nom_tests = nom_tests, 
                                             nom_medidas_biv = nom_medidas_biv, 
                                             nom_info_calif = nom_info_calif, 
                                             nom_graficas_biv = nom_graficas_biv)
      
      return(buckets_eval_var)
      
    })
    
    # juntar en una base
    info_biv_test <- data.frame()
    medidas_biv_test <- data.frame()
    graficas_biv_test <- list()
    
    for (i in 1:length(list_bucekts_eval)){
      
      bd_analisis_test %<>% left_join(list_bucekts_eval[[i]][["bd_calif"]], by = nom_id)
      
      info_biv_test %<>% rbind(list_bucekts_eval[[i]][["info_biv"]])
      medidas_biv_test %<>% rbind(list_bucekts_eval[[i]][["medidas_biv"]])
      graficas_biv_test %<>% c(list_bucekts_eval[[i]][["graficas_biv"]])
      
    }
    
    # acumular
    bd_calif_f %<>% rbind(bd_analisis_test)
    analisis[["info_biv"]][["test"]] <- info_biv_test
    analisis[["medidas_biv"]][["test"]] <- medidas_biv_test
    analisis[["graficas_biv"]][["test"]] <- graficas_biv_test
    
  }
  
  
  # salida
  analisis[["bd_calif"]] <- bd_calif_f
  analisis[["vals_esp"]] <- vesps
  analisis[["tipos_var"]] <- tipo_xs

  return(analisis)

}
