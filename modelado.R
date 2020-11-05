################################################################################
################################################################################

############ preliminares ######################################################

# opciones
options(scipen = 999)

# librerias
library(magrittr)
library(dplyr)
library(reshape2)
library(pROC)
library(MLmetrics)
library(ltm)


# Posibles valores
# tipo_y: "numerica", "dicotomica"
# nom_medidas_fit: "correlacion", "recm", "rdam", "gini", "auc", "f1", "acc", "ppv", "npv", "tpr", "tnr"
# nom_graficas_fit: "dispersion", "series", "densidades", "densidades_clasif", "boxplot_clasif", "mapa_calor"


################################################################################
################################################################################

############ Funciones para seleccion de variables #############################


#### Funcion para evaluacion o seleccion de variables ##########################

f_analisis_modelo <- function(bd_biv_mod, meds_biv, tipo_analisis = "train", modelo = NULL,
                              noms_xs = "x", nom_xs_elec = noms_xs, nom_y = "y", nom_p = "peso", nom_id = "id", 
                              tipo_y = "numerica", nom_medida_biv = "correlacion", 
                              nom_medidas_fit = c("recm", "correlacion", "rdam"),
                              nom_graficas_fit = c("dispersion", "series", "densidades"),
                              p_nas_max = 1, medida_biv_min = 0, cor_xs_max = .5, 
                              signi = .05, peso_max = .3, peso_min = .05, n_vars_min = 5,
                              cortes_clasif = seq(0, 1, by = .05)){
  
  # nombres de variables explicatibas para ajuste de modelo
  if (tipo_y == "numerica"){ noms_xs_to_fit <- paste0(noms_xs, "_ym") 
  } else if (tipo_y == "dicotomica"){ noms_xs_to_fit <- paste0(noms_xs, "_woe") }
  
  # base
  bd_biv_mod %<>% select(nom_y, nom_p, nom_id, unique(noms_xs, nom_xs_elec), noms_xs_to_fit) %>% 
                  data.frame(stringsAsFactors = F)
  
  # eliminar woes negativos de valores especiales 
  if (tipo_y == "dicotomica"){
    
    for (nom_x_fit in noms_xs_to_fit){
      
      bd_biv_mod %<>% filter(!get(nom_x_fit) %in% c(-Inf, Inf))
      
    }
    
  }
  
  
  # medidas de variables -----------------------------------------------------------------------
  
  # medidas de variables
  list_medidas_xs <- lapply(noms_xs, function(nom_x){
    
    # nom_x <- "ultimo_iot"
    
    nom_x_fit <- ifelse(tipo_y == "numerica", yes = paste0(nom_x, "_ym"), no = paste0(nom_x, "_woe"))
    
    eliminar <- 0
    
    # porcentaje de nas
    p_nas <- sum(bd_biv_mod[[nom_x]] == -9999999990)/nrow(bd_biv_mod) %>% round(5)
    
    if (p_nas > p_nas_max){ eliminar <- 1 }
      
    # desviacion
    desv <- sd(bd_biv_mod[[nom_x_fit]]) %>% round(5)
        
    if (desv == 0){ eliminar <- 1 }
    
    # medida de relacion con y
    medida_biv <- meds_biv[[nom_medida_biv]][meds_biv$variable == nom_x] %>% round(5)
        
    if (abs(medida_biv) < medida_biv_min){ eliminar <- 1 }  
      
    # juntar
    x_eval <- data.frame("variable" = nom_x, "p_nas" = p_nas, "desv" = desv,
                         stringsAsFactors = F)
    
    x_eval[[nom_medida_biv]] <- medida_biv
    x_eval[["eliminar"]] <- eliminar
    
    return(x_eval)
    
  })
  
  medidas_xs <- do.call(rbind, list_medidas_xs)
  
  # ordenar por medida de relacion con y
  medidas_xs %<>% arrange(desc(get(nom_medida_biv)), p_nas)
  
  
  # medidas entre variables -----------------------------------------------------------------------
  
  # correlaciones entre variables
  cor_xs <- bd_biv_mod %>% select(medidas_xs$variable[medidas_xs$desv > 0]) %>% cor() %>% data.frame()
  
  cor_xs$variable_1 <- row.names(cor_xs)
  row.names(cor_xs) <- NULL
  
  cor_xs %<>% gather(variable_2, cor, -variable_1) %>% mutate(cor = abs(cor))
  
  
  # analisis en test -----------------------------------------------------------------------
  
  if (tipo_analisis == "test" & !is.null(modelo)){
    
    cat(paste0(" > ", "evaluacion de modelo",
               "\n"))
    
    modelo_eval <- f_evaluacion_modelo(bd_biv_mod, tipo_analisis, modelo,
                                       nom_xs_elec, nom_y, nom_p, nom_id, 
                                       tipo_y, nom_medidas_fit,
                                       cortes_clasif)
      
   
    # analisis en train -----------------------------------------------------------------------
    
  } else if (tipo_analisis == "train"){
    
    
    # filtros -----------------------------------------------------------------------
    
    cat(paste0(" > ", "filtro de variables",
               "\n"))
    
    # filtro por medidas de variables
    medidas_xs %<>% filter(eliminar == 0) 
    
    # filtro por correlacion entre variables
    nom_xs_elec <- medidas_xs$variable[1]
    for (nom_x_iter in medidas_xs$variable[-1]){
      
      cor_xs_iter <- cor_xs %>% filter(variable_1 %in% nom_xs_elec & variable_2 == nom_x_iter)
      
      if(sum(abs(cor_xs_iter$cor) > cor_xs_max) == 0){
        
        nom_xs_elec %<>% c(nom_x_iter)
        
      }
      
    }
    
    
    # seleccion -----------------------------------------------------------------------
    
    cat(paste0(" > ", "eleccion de modelo",
               "\n"))
    
    seguir_eleccion <- T
    
    while(seguir_eleccion){
      
      modelo_eval <- f_evaluacion_modelo(bd_biv_eval = bd_biv_mod, tipo_analisis, NULL,
                                         nom_xs_eval = nom_xs_elec, nom_y, nom_p, nom_id, 
                                         tipo_y, nom_medidas_fit,
                                         cortes_clasif)
      
      if (tipo_y == "numerica"){ nom_xs_elec <- gsub("_ym", "", row.names(modelo_eval$summary_modelo))
      } else if (tipo_y == "dicotomica"){ nom_xs_elec <- gsub("_woe", "", row.names(modelo_eval$summary_modelo)) }
      
      medidas_vars_iter <- modelo_eval$medidas_vars %>% data.frame()
      
      # si no se supera el minimo de variables
      if (length(nom_xs_elec) <= n_vars_min){
        
        seguir_eleccion <- F
      
        # checar signiicancia
      } else if (sum(medidas_vars_iter$valor_p > signi) > 0){
        
        medidas_vars_iter %<>% filter(valor_p > signi) %>% arrange(desc(valor_p))
        nom_xs_elec <- nom_xs_elec[nom_xs_elec != medidas_vars_iter$variable[1]]
        
        # checar peso maximo
      } else if (sum(medidas_vars_iter$peso > peso_max) > 0){
        
        medidas_vars_iter %<>% filter(peso > peso_max) %>% arrange(desc(peso))
        nom_xs_elec <- nom_xs_elec[nom_xs_elec != medidas_vars_iter$variable[1]]
        
        # checar peso minimo
      } else if (sum(medidas_vars_iter$peso < peso_min) > 0){
        
        medidas_vars_iter %<>% filter(peso < peso_min) %>% arrange(peso)
        nom_xs_elec <- nom_xs_elec[nom_xs_elec != medidas_vars_iter$variable[1]]
        
      } else {
        
        seguir_eleccion <- F
        
      }
      
    }
    
    # valores de variables elegidas
    medidas_xs %<>% filter(variable %in% nom_xs_elec)
    cor_xs %<>% filter(variable_1 %in% nom_xs_elec & variable_2 %in% nom_xs_elec)
    
  } 
  
  
  # graficas de ajuste
  graficas_fit <- f_graficas_biv(modelo_eval$bd_fit, cor_xs,
                                 nom_y, paste0(nom_y, "_fit"),
                                 tipo_y, nom_graficas_fit)
  
  # agregar correlaciones
  analisis_mod <- modelo_eval
  
  analisis_mod$medidas_vars %<>% left_join(medidas_xs, by = "variable")
  analisis_mod$correlaciones <- cor_xs
  
  bd_biv_mod %<>% left_join(analisis_mod$bd_fit %>% select(-nom_y), by = nom_id)
  analisis_mod$bd_fit <- bd_biv_mod
  analisis_mod$graficas_fit <- graficas_fit
  
  return(analisis_mod)
  
}


################################################################################
################################################################################

############ Funciones para evaluacion de modelos ##############################


#### Funcion para evaluacion de modelo de regresion ############################

f_evaluacion_modelo <- function(bd_biv_eval, tipo_analisis = "train", modelo_fit = NULL, 
                                nom_xs_eval, nom_y = "y", nom_p = "peso", nom_id = "id",
                                tipo_y = "numerica", 
                                nom_medidas_fit = c("recm", "correlacion", "rdam"),
                                cortes_clasif = seq(0, 1, by = .05)) {
  
  # nombres de variables explicatibas para ajuste de modelo
  if (tipo_y == "numerica"){ nom_xs_eval_fit <- paste0(nom_xs_eval, "_ym") 
  } else if (tipo_y == "dicotomica"){ nom_xs_eval_fit <- paste0(nom_xs_eval, "_woe") }
  
  # variable categoria para test
  bd_biv_eval %<>% mutate(peso_test = get(nom_p))
  
  # modelo
  if (tipo_analisis == "train"){
    
    if (tipo_y == "numerica"){
      
      modelo_fit <- lm(as.formula(paste(nom_y, paste(nom_xs_eval_fit, collapse = "+"), sep = "~")), 
                   data = bd_biv_eval, weights = peso_test)  
      
    } else if (tipo_y == "dicotomica"){
      
      modelo_fit <- glm(as.formula(paste(nom_y, paste(nom_xs_eval_fit, collapse = "+"), sep = "~")), 
                    data = bd_biv_eval, family = binomial, weights = peso_test)  
      
    }
    
    
  } 
  

  # meidas de variables del modelo -----------------------------------------------------------------------
  
  # resumen de modelo
  summary_modelo <- summary(modelo_fit)$coefficients
  summary_modelo <- summary(modelo_fit)$coefficients[-1, ]
  nom_xs_eval_fit <- row.names(summary_modelo)
  
  # valores p de variables
  if (tipo_y == "numerica"){ pvals_modelo <- summary_modelo[, "Pr(>|t|)"] %>% round(5) 
  } else if (tipo_y == "dicotomica"){ pvals_modelo <- summary_modelo[, "Pr(>|z|)"] %>% round(5) }
  
  # varianzas de variables
  varianzas_modelo <- diag(var(bd_biv_eval[, nom_xs_eval_fit]))
  
  # betas y pesos de variables
  betas_modelo <- summary_modelo[, "Estimate"]
  pesos_modelo <- sqrt(varianzas_modelo*betas_modelo^2) 
  betas_modelo %<>% round(5)
  pesos_modelo <- pesos_modelo/sum(pesos_modelo) %>% round(5)
  
  # vifs de variables
  vifs_modelo <- numeric()
  for(j in 1:length(nom_xs_eval_fit)){
    
    f_vifs_iter <- as.formula(paste(nom_xs_eval_fit[j], '~', paste(nom_xs_eval_fit[-j], collapse = ' + ')))
    model_vifs_iter <- lm(f_vifs_iter, data = bd_biv_eval)
    
    vifs_modelo %<>% c(1/(1 - summary(model_vifs_iter)$r.squared))
    
  }
  
  vifs_modelo %<>% round(5)
  
  # nombres de variables explicatibas para ajuste de modelo
  if (tipo_y == "numerica"){ nom_xs_eval <- gsub("_ym", "", nom_xs_eval_fit) 
  } else if (tipo_y == "dicotomica"){ nom_xs_eval <- gsub("_woe", "", nom_xs_eval_fit) }
  
  
  # juntar 
  medidas_vars <- data.frame("variable" = nom_xs_eval, "beta" = betas_modelo, "valor_p" = pvals_modelo, 
                             "peso" = pesos_modelo, "vif" = vifs_modelo, row.names = NULL, stringsAsFactors = F)
  
  
  # meidas de variables del modelo -----------------------------------------------------------------------
  
  nom_y_fit <- paste0(nom_y, "_fit")
  
  # estimacion
  bd_biv_eval[[nom_y_fit]] <- predict(modelo_fit, newdata = bd_biv_eval, type = "response", weights = peso_test) 
  
  bd_biv_eval %<>% select(nom_id, nom_y, nom_y_fit)
  
  # medidas de ajuste
  medidas_mod <- f_medidas_fit(bd_medidas = bd_biv_eval, nom_y, nom_y_fit, tipo_y, nom_medidas_fit, cortes_clasif)
  
  # juntar
  evaluacion <- list("modelo" = modelo_fit, 
                     "summary_modelo" = summary_modelo, 
                     "medidas_modelo" = medidas_mod, 
                     "medidas_vars" = medidas_vars, 
                     "bd_fit" = bd_biv_eval)
  
  return(evaluacion)
  
}


################################################################################
################################################################################

############ funciones para calculo de estadisticas ############################


#### calculo de medidas de ajuste ##############################################

f_medidas_fit <- function(bd_medidas, 
                          nom_y = "y", nom_y_fit = "y_fit",
                          tipo_y = "numerica", 
                          nom_medidas_fit = "correlacion", 
                          cortes_clasif = seq(0, 1, by = .05)){
  
  # base para calificar
  bd_medidas %<>% select(nom_y, nom_y_fit) %>% data.frame(stringsAsFactors = F)
  
  # bases para acumular
  medidas <- data.frame("tipo_y" = tipo_y)
  
  medidas_clasif <- NULL
  
  if (sum(nom_medidas_fit %in% c("f1", "acc", "ppv", "npv", "tpr", "tnr")) > 0){
    
    list_cortes <- lapply(cortes_clasif, function(corte){
      
      casif_real <- bd_medidas[[nom_y]]
      clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
      
      tp <- sum(casif_real == 1 & clasif_corte == 1)
      tn <- sum(casif_real == 0 & clasif_corte == 0)
      fp <- sum(casif_real == 0 & clasif_corte == 1)
      fn <- sum(casif_real == 1 & clasif_corte == 0)
      
      data_opt <- data.frame("corte" = corte, 
                             "tn" = tn, "fp" = fp, "n" = tn + fp,
                             "tp" = tp, "fn" = fn, "p" = tp + fn)
      
      return(data_opt)
      
    })
    
    medidas_clasif <- do.call(rbind, list_cortes)
    
  }
  
  
  # correlacion -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "correlacion") > 0){
    
    if (tipo_y == "numerica"){
      
      medidas[["correlacion"]] <- ifelse(sd(bd_medidas[[nom_y_fit]]) == 0, yes = 0, no = abs(cor(bd_medidas[[nom_y_fit]], bd_medidas[[nom_y]])))
    
    } else if (tipo_y == "dicotomica"){
      
      medidas[["correlacion"]] <- ifelse(sd(bd_medidas[[nom_y_fit]]) == 0, yes = 0, no = abs(biserial.cor(bd_medidas[[nom_y_fit]], bd_medidas[[nom_y]], level = 2)))
      
    }
    
  }
  
  
  # recm -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "recm") > 0){
    
    medidas[["recm"]] <- -9999999990
    
    if (tipo_y == "numerica"){
      
      medidas[["recm"]] <- sqrt(mean((bd_medidas[[nom_y_fit]] - bd_medidas[[nom_y]])^2))
      
    }
    
  }
  
  
  # rdam -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "rdam") > 0){
    
    medidas[["rdam"]] <- -9999999990
    
    if (tipo_y == "numerica"){
      
      medidas[["rdam"]] <- sqrt(mean(abs(bd_medidas[[nom_y_fit]] - bd_medidas[[nom_y]])))
      
    }
    
  }
  
  
  # gini -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "gini") > 0){
    
    medidas[["gini"]] <- -9999999990
    
    if (tipo_y == "dicotomica"){
      
      medidas[["gini"]] <- auc(roc(bd_medidas[[nom_y]], bd_medidas[[nom_y_fit]]))*2 - 1
      
    }
    
  }
  
  
  # auc -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "auc") > 0){
    
    medidas[["auc"]] <- -9999999990
    
    if (tipo_y == "dicotomica"){
      
      medidas[["auc"]] <- auc(roc(bd_medidas[[nom_y]], bd_medidas[[nom_y_fit]]))
      
    }
    
  }
  
  
  # score f1 con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "f1") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
         
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        if (length(unique(clasif_corte)) == 1){ f1 <- NA
        } else { f1 <- F1_Score(casif_real, clasif_corte) }
        
        data_opt <- data.frame("corte" = corte, "f1" = f1)
        
        return(data_opt)
        
      })
      
      f1_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(f1_opt, by = "corte")
      
    }
    
  }
  
  
  # acc con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "acc") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
        
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        tp <- sum(casif_real == 1 & clasif_corte == 1)
        tn <- sum(casif_real == 0 & clasif_corte == 0)
        fp <- sum(casif_real == 0 & clasif_corte == 1)
        fn <- sum(casif_real == 1 & clasif_corte == 0)
        
        acc <- (tp + tn)/(tp + tn + fp + fn)
        
        data_opt <- data.frame("corte" = corte, "acc" = acc)
        
        return(data_opt)
        
      })
      
      acc_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(acc_opt, by = "corte")
      
    }
    
  }
  
  
  # ppv con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "ppv") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
        
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        tp <- sum(casif_real == 1 & clasif_corte == 1)
        tn <- sum(casif_real == 0 & clasif_corte == 0)
        fp <- sum(casif_real == 0 & clasif_corte == 1)
        fn <- sum(casif_real == 1 & clasif_corte == 0)
        
        ppv <- tp/(tp + fp)
        
        data_opt <- data.frame("corte" = corte, "ppv" = ppv)
        
        return(data_opt)
        
      })
      
      ppv_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(ppv_opt, by = "corte")
      
    }
    
  }
  
  
  # npv con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "npv") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
        
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        tp <- sum(casif_real == 1 & clasif_corte == 1)
        tn <- sum(casif_real == 0 & clasif_corte == 0)
        fp <- sum(casif_real == 0 & clasif_corte == 1)
        fn <- sum(casif_real == 1 & clasif_corte == 0)
        
        npv <- tn/(tn + fn)
        
        data_opt <- data.frame("corte" = corte, "npv" = npv)
        
        return(data_opt)
        
      })
      
      npv_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(npv_opt, by = "corte")
      
    }
    
  }
  
  
  # tpr con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "tpr") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
        
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        tp <- sum(casif_real == 1 & clasif_corte == 1)
        tn <- sum(casif_real == 0 & clasif_corte == 0)
        fp <- sum(casif_real == 0 & clasif_corte == 1)
        fn <- sum(casif_real == 1 & clasif_corte == 0)
        
        tpr <- tp/(tp + fn)
        
        data_opt <- data.frame("corte" = corte, "tpr" = tpr)
        
        return(data_opt)
        
      })
      
      tpr_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(tpr_opt, by = "corte")
      
    }
    
  }
  
  
  # tnr con mejor punto de corte -----------------------------------------------------------------------
  
  if (sum(nom_medidas_fit == "tnr") > 0){
    
    if (tipo_y == "dicotomica"){
      
      list_cortes <- lapply(cortes_clasif, function(corte){
        
        casif_real <- bd_medidas[[nom_y]]
        clasif_corte <- as.numeric(bd_medidas[[nom_y_fit]] > corte)
        
        tp <- sum(casif_real == 1 & clasif_corte == 1)
        tn <- sum(casif_real == 0 & clasif_corte == 0)
        fp <- sum(casif_real == 0 & clasif_corte == 1)
        fn <- sum(casif_real == 1 & clasif_corte == 0)
        
        tnr <- tn/(tn + fp)
        
        data_opt <- data.frame("corte" = corte, "tnr" = tnr)
        
        return(data_opt)
        
      })
      
      tnr_opt <- do.call(rbind, list_cortes)
      
      medidas_clasif %<>% left_join(tnr_opt, by = "corte")
      
    }
    
  }
  
  
  meds <- list("medidas_fit" = medidas %>% select(-tipo_y), 
               "medidas_clasif" = medidas_clasif)
  
  return(meds)
  
}


################################################################################
################################################################################

############ funciones para graficas ###########################################


#### graficas de ajuste ########################################################

f_graficas_biv <- function(bd_graf, bd_cor_graf = NULL,
                           nom_y = "y", nom_y_fit = "y_fit",
                           tipo_y = "numerica", 
                           nom_graficas_fit = c("dispersion", "series", "densidades")){
 
  graficas <- list()
  
  
  # grafica de dispersion -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "dispersion")){
    
    graficas$dispersion <- ggplot(data = bd_graf) +
                           geom_point(aes(x = get(nom_y), y = get(nom_y_fit)), alpha = .2, color = "#191970", group = 1) +
                           theme_minimal() +
                           labs(title = "Dispersion", subtitle = "Observado vs ajustado", x = "observado", y = "ajustado") +
                           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de series -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "series")){
    
    graficas$series <- ggplot(data = bd_graf %>% mutate(num_obs = row_number())) +
                       geom_point(aes(x = num_obs, y = get(nom_y), color = "observado"), alpha = .8, group = 1) +
                       geom_point(aes(x = num_obs, y = get(nom_y_fit), color = "ajustado"), alpha = .8, group = 1) +
                       theme_minimal() +
                       scale_color_manual(name = "tipo", values = c("observado" = "#191970", "ajustado" = "#006EC1")) +
                       labs(title = "Series", subtitle = "Observado vs ajustado", x = "numero de observacion", y = "y") +
                       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de densidades -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "densidades")){
    
    graficas$densidades <- ggplot(data = bd_graf) +
                           geom_density(aes(get(nom_y), fill = "observado"), alpha = .8) +
                           geom_density(aes(get(nom_y_fit), fill = "ajustado"), alpha = .8) +
                           theme_minimal() +
                           scale_fill_manual(name = "tipo", values = c("observado" = "#191970", "ajustado" = "#006EC1")) +
                           labs(title = "Densidades", subtitle = "Observado vs ajustado", x = "y", y = "densidad") +
                           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de densidades densidades por valor observado -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "densidades_clasif") & tipo_y == "dicotomica"){
    
    graficas$densidades_clasif <- ggplot() +
                                  geom_density(data = bd_graf %>% filter(get(nom_y) == 1), aes(get(nom_y_fit), fill = "positivos"), alpha = .8) +
                                  geom_density(data = bd_graf %>% filter(get(nom_y) == 0), aes(get(nom_y_fit), fill = "negativos"), alpha = .8) +
                                  theme_minimal() +
                                  scale_fill_manual(name = "tipo", values = c("positivos" = "#191970", "negativos" = "#006EC1")) +
                                  labs(title = "Densidades", subtitle = "Ajustado vs clasificacion observada", x = "y", y = "densidad") +
                                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # grafica de boxplot por valor observado -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "boxplot_clasif")  & tipo_y == "dicotomica"){
    
    graficas$boxplot_clasif <- ggplot(data = bd_graf) +
                               geom_boxplot(aes(x = factor(get(nom_y)), y = get(nom_y_fit)), color = "#191970", group = 1) +
                               theme_minimal() +
                               labs(title = "Boxplot", subtitle = "Ajustado vs clasificacion observada", x = "y", y = "distribucion") +
                               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
  }
  
  
  # mapa de calor de correlaciones -----------------------------------------------------------------------
  
  if (sum(nom_graficas_fit == "mapa_calor") & !is.null(bd_cor_graf)){
    
    graficas$mapa_calor <- ggplot(data = bd_cor_graf, aes(variable_2, variable_1, fill = cor)) +
                           geom_tile(color = "white") +
                           theme_minimal() +
                           scale_fill_gradient2(low = "yellow", high = "red", mid = "white", 
                                                 midpoint = 0, limit = c(-1,1), space = "Lab", 
                                                 name = "correlacion") +
                           labs(title = "Mapa de calor", subtitle = "Correlaciones", x = "", y = "") +
                           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                           coord_fixed()
    
  }
  
  return(graficas)
   
}


################################################################################
################################################################################

############ Funciones para modelado ###########################################


#### Eleccion y evaluacion de modelos ##########################################

f_analisis_modelado <- function(bd_modelado, medidas_biv, nom_analisis = c("train", "test"), modelo = NULL,
                                noms_xs = "x", nom_xs_elec = noms_xs, nom_y = "y", nom_p = "peso", nom_id = "id", nom_b = "base",
                                tipo_y = "numerica", nom_medida_biv = "correlacion", 
                                nom_medidas_fit = c("recm", "correlacion", "rdam"),
                                nom_graficas_fit = c("dispersion", "series", "densidades"),
                                p_nas_max = 1, medida_biv_min = 0, cor_xs_max = .5, 
                                signi = .05, peso_max = .3, peso_min = .05, n_vars_min = 5,
                                cortes_clasif = seq(0, 1, by = .05)){
  
  
  cat(paste0("----------------------------------------",
             "\n", "----------------------------------------",
             "\n", 
             "\n", "Analisis de modelado",
             "\n"))
  
  
  # nombres de variables explicatibas para ajuste de modelo
  if (tipo_y == "numerica"){ noms_xs_to_fit <- paste0(noms_xs, "_ym") 
  } else if (tipo_y == "dicotomica"){ noms_xs_to_fit <- paste0(noms_xs, "_woe") }
  
  # base
  bd_modelado %<>% select(nom_y, nom_p, nom_id, nom_b, noms_xs, noms_xs_to_fit) %>% 
                   data.frame(stringsAsFactors = F)
  
  # bases train, test y oot
  bd_modelado_train <- bd_modelado %>% filter(get(nom_b) == "train")
  bd_modelado_test <- bd_modelado %>% filter(get(nom_b) == "test")
  bd_modelado_oot <- bd_modelado %>% filter(get(nom_b) == "oot")
  
  medidas_biv_train <- medidas_biv$train
  medidas_biv_test <- medidas_biv$test
  medidas_biv_oot <- medidas_biv$oot
  
  analisis <- list()
  bd_fit_f <- data.frame()
  
  # analisis en train -----------------------------------------------------------------------
  
  if (nrow(bd_modelado_train) > 0 & sum(nom_analisis == "train") > 0){
    
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Analisis en train",
               "\n",
               "\n"))
    
    # esitmacion de modelo
    analisis_train <- f_analisis_modelo(bd_biv_mod = bd_modelado_train, 
                                        meds_biv = medidas_biv_train, 
                                        tipo_analisis = "train", 
                                        modelo = NULL,
                                        noms_xs, 
                                        nom_xs_elec,
                                        nom_y, 
                                        nom_p, 
                                        nom_id, 
                                        tipo_y, 
                                        nom_medida_biv, 
                                        nom_medidas_fit,
                                        nom_graficas_fit,
                                        p_nas_max, 
                                        medida_biv_min, 
                                        cor_xs_max, 
                                        signi, 
                                        peso_max, 
                                        peso_min, 
                                        n_vars_min,
                                        cortes_clasif)
    
    # info para calificar test
    modelo <- analisis_train$modelo
    nom_xs_elec <- analisis_train$medidas_vars$variable
    
    # acumular
    bd_fit_f %<>% rbind(analisis_train$bd_fit)
    
    analisis[["modelo"]] <- analisis_train$modelo
    analisis[["summary_modelo"]] <- analisis_train$summary_modelo
    analisis[["medidas_modelo"]][["train"]] <- analisis_train$medidas_modelo
    analisis[["medidas_vars"]][["train"]] <- analisis_train$medidas_vars
    analisis[["correlaciones"]][["train"]] <- analisis_train$correlaciones
    analisis[["graficas_fit"]][["train"]] <- analisis_train$graficas_fit
    
  }
  
  
  # analisis en test -----------------------------------------------------------------------
  
  if (nrow(bd_modelado_test) > 0 & sum(nom_analisis == "test") > 0 & !is.null(modelo)){
    
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Analisis en test",
               "\n",
               "\n"))
    
    # esitmacion de modelo
    analisis_test <- f_analisis_modelo(bd_biv_mod = bd_modelado_test, 
                                       meds_biv = medidas_biv_test, 
                                       tipo_analisis = "test", 
                                       modelo,
                                       noms_xs, 
                                       nom_xs_elec,
                                       nom_y, 
                                       nom_p, 
                                       nom_id, 
                                       tipo_y, 
                                       nom_medida_biv, 
                                       nom_medidas_fit,
                                       nom_graficas_fit,
                                       p_nas_max, 
                                       medida_biv_min, 
                                       cor_xs_max, 
                                       signi, 
                                       peso_max, 
                                       peso_min, 
                                       n_vars_min,
                                       cortes_clasif)
      
    # acumular
    bd_fit_f %<>% rbind(analisis_test$bd_fit)
    
    analisis[["modelo"]] <- analisis_test$modelo
    analisis[["summary_modelo"]] <- analisis_test$summary_modelo
    analisis[["medidas_modelo"]][["test"]] <- analisis_test$medidas_modelo
    analisis[["medidas_vars"]][["test"]] <- analisis_test$medidas_vars
    analisis[["correlaciones"]][["test"]] <- analisis_test$correlaciones
    analisis[["graficas_fit"]][["test"]] <- analisis_test$graficas_fit
    
  }
  
  
  # analisis en oot -----------------------------------------------------------------------
  
  if (nrow(bd_modelado_oot) > 0 & sum(nom_analisis == "oot") > 0 & !is.null(modelo)){
    
    cat(paste0("\n", "----------------------------------------",
               "\n", 
               "\n", "Analisis en oot",
               "\n",
               "\n"))
    
  }
  
  cat(paste0("\n", "----------------------------------------",
             "\n", "----------------------------------------",
             "\n",
             "\n"))
  
  # salida
  analisis[["bd_fit"]] <- bd_fit_f
  
  return(analisis)
  
}