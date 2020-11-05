source("01 Auxiliares/01 Codigos/trameado.R")

data <- data.frame("id" = 1:5000,
                   "num1" = c(rnorm(1000, 18000, 700),
                              rnorm(1000, 15000, 700),
                              rnorm(1000, 10000, 700),
                              rnorm(1000, 5000, 700),
                              rnorm(1000, 1000, 700)),
                   "num2" = c(rpois(1000, 1000),
                              rpois(1000, 3000),
                              rpois(1000, 6000),
                              rpois(1000, 9000),
                              rpois(1000, 12000)),
                   "ordi1" = c(rep("a", 750),
                               rep("d", 250),
                               rep("e", 250),
                               rep("f", 250),
                               rep("g", 1250),
                               rep("l", 250),
                               rep("m", 250),
                               rep("n", 250),
                               rep("o", 250),
                               rep("p", 250),
                               rep("q", 250),
                               rep("r", 250),
                               rep("s", 250),
                               rep("t", 250)),
                   "nomi1" = c(sample(c("z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10"), size = 1500, replace = T),
                               sample(c("z6", "z7", "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15"), size = 2000, replace = T),
                               sample(c("z11", "z12", "z13", "z14", "z15", "z16", "z17", "z18", "z19", "z20"), size = 1500, replace = T)),
                   "y" = c(round(runif(1000, 0, .55)),
                           round(runif(1000, 0, .8)),
                           round(runif(1000, 0, 1)),
                           round(runif(1000, .2, 1)),
                           round(runif(1000, .45, 1))),
                   "peso" = round(runif(5000, 0, 10000)),
                   "ref" = round(runif(5000, 1000, 100000)),
                   "base" = c(rep("train", 4000), rep("test", 1000)))


data$nomi1 <- as.character(data$nomi1)
data$ordi1 <- as.character(data$ordi1)
data$nomi1[sample(1:nrow(data), 500, replace = F)] <- "-999999"
data$ordi1[sample(1:nrow(data), 500, replace = F)] <- "-999999"
data$num1[sample(1:nrow(data), 100, replace = F)] <- -999999

list_esp <- list("nomi1" = "-999999",
                 "ordi1" = "-999999",
                 "num1" = -999999)

tipo_xs <- list("num1" = "numerica", 
                "num2" = "numerica", 
                "ordi1" = "ordinal", 
                "nomi1" = "nominal")

analisis_prueba <- f_analisis_biv(bd_analisis = data, 
                                  nom_analisis = c("train", "test"),
                                  info_biv_eval = NULL, 
                                  noms_xs = c("num1", "num2", "ordi1", "nomi1"), 
                                  nom_y = "y", 
                                  nom_p = "peso", 
                                  nom_r = "ref", 
                                  nom_id = "id", 
                                  nom_b = "base",
                                  vesp_xs = list_esp,
                                  pobl_prop_min = .05, 
                                  signi = .1,
                                  tipo_xs = tipo_xs, 
                                  tipo_y = "dicotomica", 
                                  nom_tests = c("welch", "wald"), 
                                  nom_medidas_biv = "gini", 
                                  nom_info_calif = c("bucket", "ym", "woe"), 
                                  nom_graficas_biv = "bivariado")

