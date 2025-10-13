#--------------------------------------------------
# DESARROLLO DE LA PROBLEMÁTICA: CASOS DE DENGUE  #
#--------------------------------------------------

required_pkgs <- c("sf","spdep","dplyr","readxl","stringi","tmap","ggplot2","MASS","matrixStats", "stringr", "tmaptools", "RColorBrewer", "ggpubr",
                   "gridExtra", "tidyr", "spatialreg", "Metrics", "FNN", "viridis")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"]) ]
if(length(new_pkgs)) install.packages(new_pkgs)

library(sf)
library(spdep)
library(dplyr)
library(readxl)
library(stringi)
library(tmap)
library(tmaptools)
library(RColorBrewer)
library(ggpubr)
library(gridExtra) 
library(ggplot2)
library(MASS)
library(matrixStats)
library(stringr)
library(tidyr)
library(spatialreg)
library(Metrics)
library(FNN)
library(viridis)

setwd("~/project_2/datasets/datasets/Situación_3") 

datos <- st_read("datos_dengue.shp", quiet=TRUE)

#EDA CLÁSICO:

# A. Histograma + densidad de la tasa de dengue
g1 <- ggplot(datos %>% st_drop_geometry(), aes(x = tas_dng)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#56B4E9", alpha = 0.7, color = "black") +
  geom_density(color = "#D55E00", size = 1, na.rm=TRUE) +
  labs(title = "Distribución de la tasa de dengue por barrio", x = "Casos por 10,000 hab.", y = "Densidad") +
  theme_minimal()

# B. Boxplot de la tasa de dengue
g2 <- ggplot(datos %>% st_drop_geometry(), aes(y = tas_dng)) +
  geom_boxplot(fill = "#D55E00", color = "black", na.rm = TRUE) +
  labs(title = "Boxplot tasa de dengue", y = "Casos por 10,000 hab.") +
  theme_minimal()

# C. Dispersión entre tasa de dengue y sumideros con larvas
g3 <- ggplot(datos %>% st_drop_geometry(), aes(x = NSUMIDE, y = tas_dng)) +
  geom_point(alpha = 0.7, color = "#009E73") +
  geom_smooth(method = "lm", color = "red", se = FALSE, na.rm=TRUE) +
  labs(title = "Tasa de dengue vs Sumideros con larvas", x = "Sumideros con larvas", y = "Casos por 10,000 hab.") +
  theme_minimal()

# D. Mapa coroplético de la tasa de dengue
g4 <- ggplot(datos) +
  geom_sf(aes(fill = tas_dng), color = NA) +
  scale_fill_viridis_c(option = "inferno", na.value = "grey90") +
  labs(title = "Mapa: Tasa dengue por barrio", fill = "Casos/10,000 hab") +
  theme_minimal()

grid.arrange(g1, g2, g3, g4, ncol=2)

#Pasamos al análisis espacial de los casos de dengue por barrio, con la asociación espacial a nivel global y local.
#Para esto, debemos de tener las matrices de pesos por vecindad y las haremos de tres formas: Queen, Rook y KNN.
#Establecemos cuántos k usar con el método del codo

coords <- st_coordinates(st_centroid(datos))

kmax <- 15

distancias_k <- knn.dist(coords, k = kmax)

distancias_medias <- apply(distancias_k, 2, mean, na.rm=TRUE)

plot(1:kmax, distancias_medias, type='b', pch=19,
     xlab="Número de vecinos (k)",
     ylab="Distancia media al k-ésimo vecino",
     main="Método del codo para KNN")

#Se identificó el codo entre K = 3 o 4, se escoge k= 4. Ahora se establecen las matrices:

# Queen
nb_queen <- poly2nb(datos, queen=TRUE)
listw_queen <- nb2listw(nb_queen, style="W", zero.policy=TRUE)

# Rook
nb_rook <- poly2nb(datos, queen=FALSE)
listw_rook <- nb2listw(nb_rook, style="W", zero.policy=TRUE)

# KNN
k_optimo <- 4
nb_knn <- knn2nb(knearneigh(coords, k=k_optimo))
listw_knn <- nb2listw(nb_knn, style="W", zero.policy=TRUE)

#Ahora, con la tasa de dengue, sacamos el índice global de Moran y Geary para cada una de las matrices y comparamos

moran_queen <- moran.test(datos$tas_dng, listw_queen, zero.policy=TRUE)
moran_rook  <- moran.test(datos$tas_dng, listw_rook,  zero.policy=TRUE)
moran_knn   <- moran.test(datos$tas_dng, listw_knn,   zero.policy=TRUE)

geary_queen <- geary.test(datos$tas_dng, listw_queen, zero.policy=TRUE)
geary_rook  <- geary.test(datos$tas_dng, listw_rook,  zero.policy=TRUE)
geary_knn   <- geary.test(datos$tas_dng, listw_knn,   zero.policy=TRUE)

# Tabla 1. Índice Global de Moran:

tabla_moran <- data.frame(
  Spatial_Weight_Matrix = c("Queen Contiguity (vertices and edges)",
                            "Rook Contiguity (edges only)",
                            paste0("K-Nearest Neighbors (k=", k_optimo, ")")),
  Moran_I = c(moran_queen$estimate[1], moran_rook$estimate[1], moran_knn$estimate[1]),
  p_value = c(moran_queen$p.value, moran_rook$p.value, moran_knn$p.value)
)
print(tabla_moran)

# Tabla 2. La C de Geary:

tabla_geary <- data.frame(
  Spatial_Weight_Matrix = c("Queen Contiguity (vertices and edges)",
                            "Rook Contiguity (edges only)",
                            paste0("K-Nearest Neighbors (k=", k_optimo, ")")),
  Geary_C = c(geary_queen$estimate[1], geary_rook$estimate[1], geary_knn$estimate[1]),
  p_value = c(geary_queen$p.value, geary_rook$p.value, geary_knn$p.value)
)
print(tabla_geary)

#Calculamos el índice local de Moran de la tasa de dengue

lisa_queen <- localmoran(datos$tas_dng, listw_queen, zero.policy=TRUE)
lisa_rook <- localmoran(datos$tas_dng, listw_rook, zero.policy=TRUE)
lisa_knn <- localmoran(datos$tas_dng, listw_knn, zero.policy=TRUE)

#Los ponemos en el objeto espacial para hacer el mapa

datos$lisa_q_I <- lisa_queen[,1]
datos$lisa_q_p <- lisa_queen[,5]
datos$lisa_r_I <- lisa_rook[,1]
datos$lisa_r_p <- lisa_rook[,5]
datos$lisa_k_I <- lisa_knn[,1]
datos$lisa_k_p <- lisa_knn[,5]

g_lisa_queen <- ggplot(datos) +
  geom_sf(aes(fill = lisa_q_I), color=NA) +
  scale_fill_viridis(option="D") +
  labs(title = "LISA (Queen)", fill = "I local") +
  theme_minimal()

g_lisa_rook <- ggplot(datos) +
  geom_sf(aes(fill = lisa_r_I), color=NA) +
  scale_fill_viridis(option="C") +
  labs(title = "LISA (Rook)", fill = "I local") +
  theme_minimal()

g_lisa_knn <- ggplot(datos) +
  geom_sf(aes(fill = lisa_k_I), color=NA) +
  scale_fill_viridis(option="B") +
  labs(title = "LISA (KNN, k=4)", fill = "I local") +
  theme_minimal()

grid.arrange(g_lisa_queen, g_lisa_rook, g_lisa_knn, ncol=3)

#Graficamos también con la significancia o p-value

datos$lisa_q_sig <- datos$lisa_q_p < 0.05
datos$lisa_r_sig <- datos$lisa_r_p < 0.05
datos$lisa_k_sig <- datos$lisa_k_p < 0.05

g_sig_q <- ggplot(datos) +
  geom_sf(aes(fill = lisa_q_sig), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (Queen)", fill="Significativo") +
  theme_minimal()

g_sig_r <- ggplot(datos) +
  geom_sf(aes(fill = lisa_r_sig), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (Rook)", fill="Significativo") +
  theme_minimal()

g_sig_k <- ggplot(datos) +
  geom_sf(aes(fill = lisa_k_sig), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (KNN, k=4)", fill="Significativo") +
  theme_minimal()

grid.arrange(g_sig_q, g_sig_r, g_sig_k, ncol=3)

#Haremoslo mismo con el número de sumideros,primero calcular el índice de Moran y C de Geary global

moran_queen_sum <- moran.test(datos$NSUMIDE, listw_queen, zero.policy=TRUE)
moran_rook_sum  <- moran.test(datos$NSUMIDE, listw_rook,  zero.policy=TRUE)
moran_knn_sum   <- moran.test(datos$NSUMIDE, listw_knn,   zero.policy=TRUE)

geary_queen_sum <- geary.test(datos$NSUMIDE, listw_queen, zero.policy=TRUE)
geary_rook_sum  <- geary.test(datos$NSUMIDE, listw_rook,  zero.policy=TRUE)
geary_knn_sum   <- geary.test(datos$NSUMIDE, listw_knn,   zero.policy=TRUE)

# Tabla 3. Índice Global de Moran para sumideros

tabla_moran_sum <- data.frame(
  Spatial_Weight_Matrix = c("Queen Contiguity (vertices and edges)",
                            "Rook Contiguity (edges only)",
                            paste0("K-Nearest Neighbors (k=", k_optimo, ")")),
  Moran_I = c(moran_queen_sum$estimate[1], moran_rook_sum$estimate[1], moran_knn_sum$estimate[1]),
  p_value = c(moran_queen_sum$p.value, moran_rook_sum$p.value, moran_knn_sum$p.value)
)
print(tabla_moran_sum)

# Tabla 4. La C de Geary para sumideros

tabla_geary_sum <- data.frame(
  Spatial_Weight_Matrix = c("Queen Contiguity (vertices and edges)",
                            "Rook Contiguity (edges only)",
                            paste0("K-Nearest Neighbors (k=", k_optimo, ")")),
  Geary_C = c(geary_queen_sum$estimate[1], geary_rook_sum$estimate[1], geary_knn_sum$estimate[1]),
  p_value = c(geary_queen_sum$p.value, geary_rook_sum$p.value, geary_knn_sum$p.value)
)
print(tabla_geary_sum)

#Pasamos al índice local de Moran LISA 

lisa_queen_sum <- localmoran(datos$NSUMIDE, listw_queen, zero.policy=TRUE)
lisa_rook_sum  <- localmoran(datos$NSUMIDE, listw_rook,  zero.policy=TRUE)
lisa_knn_sum   <- localmoran(datos$NSUMIDE, listw_knn,   zero.policy=TRUE)

#Y se guardan en el objeto espacial

datos$lisa_q_I_sum <- lisa_queen_sum[,1]
datos$lisa_q_p_sum <- lisa_queen_sum[,5]
datos$lisa_r_I_sum <- lisa_rook_sum[,1]
datos$lisa_r_p_sum <- lisa_rook_sum[,5]
datos$lisa_k_I_sum <- lisa_knn_sum[,1]
datos$lisa_k_p_sum <- lisa_knn_sum[,5]

#Mapa del índice I local

g_lisa_queen_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_q_I_sum), color=NA) +
  scale_fill_viridis(option="D") +
  labs(title = "LISA (Queen) - Sumideros", fill = "I local") +
  theme_minimal()

g_lisa_rook_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_r_I_sum), color=NA) +
  scale_fill_viridis(option="C") +
  labs(title = "LISA (Rook) - Sumideros", fill = "I local") +
  theme_minimal()

g_lisa_knn_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_k_I_sum), color=NA) +
  scale_fill_viridis(option="B") +
  labs(title = "LISA (KNN, k=4) - Sumideros", fill = "I local") +
  theme_minimal()

gridExtra::grid.arrange(g_lisa_queen_sum, g_lisa_rook_sum, g_lisa_knn_sum, ncol=3)

#Mapas por significancia, donde se identifican clusters significativos

datos$lisa_q_sig_sum <- datos$lisa_q_p_sum < 0.05
datos$lisa_r_sig_sum <- datos$lisa_r_p_sum < 0.05
datos$lisa_k_sig_sum <- datos$lisa_k_p_sum < 0.05

g_sig_q_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_q_sig_sum), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (Queen) - Sumideros", fill="Significativo") +
  theme_minimal()

g_sig_r_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_r_sig_sum), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (Rook) - Sumideros", fill="Significativo") +
  theme_minimal()

g_sig_k_sum <- ggplot(datos) +
  geom_sf(aes(fill = lisa_k_sig_sum), color=NA) +
  scale_fill_manual(values=c("grey90", "red")) +
  labs(title = "Barrios con LISA significativo (KNN, k=4) - Sumideros", fill="Significativo") +
  theme_minimal()

gridExtra::grid.arrange(g_sig_q_sum, g_sig_r_sum, g_sig_k_sum, ncol=3)

#MODELADO ESPACIAL USANDO SAR, SAC Y SEM, junto con las matrices de vecindad:

# Modelos SAR
sar_queen <- lagsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_queen, zero.policy=TRUE)
sar_rook  <- lagsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_rook, zero.policy=TRUE)
sar_knn   <- lagsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_knn, zero.policy=TRUE)

# Modelos SEM
sem_queen <- errorsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_queen, zero.policy=TRUE)
sem_rook  <- errorsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_rook, zero.policy=TRUE)
sem_knn   <- errorsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_knn, zero.policy=TRUE)

# Modelos SAC
sac_queen <- sacsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_queen, zero.policy=TRUE)
sac_rook  <- sacsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_rook, zero.policy=TRUE)
sac_knn   <- sacsarlm(tas_dng ~ NSUMIDE, data=datos, listw=listw_knn, zero.policy=TRUE)

#Hacemos una función para sacar las métricas Pseudo R2, AIC, Loglikelihood, y COeficiente Espacial Wy

#Primero calculamos el pseudo R2 usando el Loglikelihood del modelo y de un modelo nulo (solo con el intercepto)

null_queen <- lagsarlm(tas_dng ~ 1, data = datos, listw = listw_queen, zero.policy = TRUE)
null_rook  <- lagsarlm(tas_dng ~ 1, data = datos, listw = listw_rook,  zero.policy = TRUE)
null_knn   <- lagsarlm(tas_dng ~ 1, data = datos, listw = listw_knn,   zero.policy = TRUE)

#Organizamos los modelos en listas

mod_sar <- list(Queen = sar_queen, Rook = sar_rook, KNN = sar_knn)
mod_sem <- list(Queen = sem_queen, Rook = sem_rook, KNN = sem_knn)
mod_sac <- list(Queen = sac_queen, Rook = sac_rook, KNN = sac_knn)
mod_null <- list(Queen = null_queen, Rook = null_rook, KNN = null_knn)

#Hacemos la función para extraer métricas automáticamente:

get_all_metrics <- function(model, model_null, type="SAR") {
  s <- summary(model)
  logLik_model <- as.numeric(logLik(model))
  logLik_null  <- as.numeric(logLik(model_null))
  pseudoR2     <- 1 - (logLik_model/logLik_null)
  
  coef_esp <- if(type=="SAR") s$rho
  else if(type=="SEM") s$lambda
  else paste0("rho: ", round(s$rho, 4), ", lambda: ", round(s$lambda, 4))
  
  data.frame(
    AIC          = round(AIC(model), 2),
    logLik       = round(logLik_model, 2),
    Pseudo_R2    = round(pseudoR2, 4),
    Spatial_Coef = coef_esp
  )
}

#Juntar todos los modelos en la tabla

resultados <- do.call(rbind, lapply(
  c("Queen", "Rook", "KNN"), function(mat) {
    rbind(
      data.frame(Model="SAR", Vecindad=mat, get_all_metrics(mod_sar[[mat]], mod_null[[mat]], "SAR")),
      data.frame(Model="SEM", Vecindad=mat, get_all_metrics(mod_sem[[mat]], mod_null[[mat]], "SEM")),
      data.frame(Model="SAC", Vecindad=mat, get_all_metrics(mod_sac[[mat]], mod_null[[mat]], "SAC"))
    )
  }
))

rownames(resultados) <- NULL
print(resultados)


#Extraemos los residuos del modelo escogido SEM-KNN y los evaluamos con el índice de Moran

res_sem_knn <- residuals(sem_knn)

moran_resid_sem_knn <- moran.test(res_sem_knn, listw_knn, zero.policy=TRUE)
print(moran_resid_sem_knn)

#Mapa de los residuos:

datos$resid_sem_knn <- res_sem_knn

ggplot(datos) +
  geom_sf(aes(fill = resid_sem_knn), color=NA) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  labs(title="Residuos del modelo SEM-KNN", fill="Residuo") +
  theme_minimal()

#Finalmente sacamos las predicciones de la tasa de dengue:

datos$tas_dng_pred <- as.numeric(predict(sem_knn))

ggplot(datos) +
  geom_sf(aes(fill = tas_dng_pred), color=NA) +
  scale_fill_viridis(option="C") +
  labs(title = "Predicción de tasa de dengue por SEM-KNN", fill = "Casos/10,000 hab. (pred)") +
  theme_minimal()

#Y calculamos métricas de error de la predicción

residuos <- datos$tas_dng - datos$tas_dng_pred

RMSE <- sqrt(mean(residuos^2, na.rm = TRUE))
MAE <- mean(abs(residuos), na.rm = TRUE)
R2 <- 1 - sum(residuos^2, na.rm = TRUE) / sum((datos$tas_dng - mean(datos$tas_dng, na.rm = TRUE))^2, na.rm = TRUE)

cat("RMSE:", RMSE, "\nMAE:", MAE, "\nR^2:", R2, "\n")
