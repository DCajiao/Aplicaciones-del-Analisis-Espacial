### Actividad - Patrones Puntuales
### David Alejandro Cajiao Lazt
##################################

# install.packages(c("readxl","sf","sp","spatstat","spatstat.model","spdep","grDevices","grid","stars"))
install.packages("dplyr")

library(readxl)
library(sf)
library(sp)
library(spdep)
library(spatstat)          # meta-paquete; trae spatstat.geom/core
library(spatstat.model)
library(grDevices); library(grid)
library(stars)
library(readr)

# ----- Rutas -----
csv_path  <- "C:/Users/david/Desktop/Proyecto 2/Aplicaciones-del-Analisis-Espacial/casos/Situación_2/data/Accidentes_2009_2010.csv"
map_dir    <- "C:/Users/david/Desktop/Proyecto 2/Aplicaciones-del-Analisis-Espacial/casos/Situación_2/data/Mapas"

# =========================================================
# 1) CARGA DE DATOS
# =========================================================
raw <- read_csv(csv_path)

# Columnas de interés (según indicación)
df <- data.frame(
  x    = as.numeric(raw[["coordenada X (metros)"]]),
  y    = as.numeric(raw[["coordenada Y (metros)"]]),
  sexo = as.factor(raw[["SEXO"]]),
  tipo = as.factor(raw[["CONDICCION"]])
)

# Filas válidas
df <- subset(df, is.finite(x) & is.finite(y))
df$id <- seq_len(nrow(df))   # para rastrear filas

# =========================================================
# 2) MAPAS Y VENTANA DE ESTUDIO 
# =========================================================
shps <- list.files(map_dir, pattern = "\\.shp$", full.names = TRUE)
stopifnot(length(shps) > 0)

cand_borde <- grep("bord|borde|limite|límite|boundary",
                   basename(shps), ignore.case = TRUE, value = TRUE)
borde_path <- if (length(cand_borde) > 0) file.path(map_dir, cand_borde[1]) else shps[1]

borde_sf <- st_read(borde_path, quiet = TRUE) |> st_make_valid()
stopifnot(!is.na(st_crs(borde_sf)))  # debe traer CRS definido (en tu caso EPSG:32618)

# Unificar geometría del borde a un solo polígono
borde_sf <- st_union(borde_sf) |> st_as_sf()
borde_crs <- st_crs(borde_sf)

# =========================================================
# 3) CRS DE LOS PUNTOS (detección entre 3114/3115/3116) y transformación al borde
# =========================================================
mk_pts_src <- function(epsg_src) {
  st_as_sf(df, coords = c("x","y"), crs = epsg_src) |>
    st_transform(borde_crs)
}
cands_epsg <- c(3114, 3115, 3116)  # Bogotá West / Central / East
pts_list <- lapply(cands_epsg, mk_pts_src)

inside_count <- sapply(pts_list, function(p) {
  sum(st_within(p, st_geometry(borde_sf), sparse = FALSE))
})
names(inside_count) <- paste0("EPSG:", cands_epsg)
print(inside_count)

best_ix  <- which.max(inside_count)
best_eps <- cands_epsg[best_ix]
cat(">>> Mejor EPSG para los puntos:", best_eps, 
    "| puntos dentro (previo a filtrado estricto):", inside_count[best_ix], "\n")
stopifnot(inside_count[best_ix] > 0)

pts_sf <- pts_list[[best_ix]]  # puntos ya en CRS del borde, con columnas sexo/tipo/id

# =========================================================
# 4) FILTRADO AL INTERIOR + MANEJO DE DUPLICADOS
# =========================================================
# (a) Tolerancia de borde (buffer pequeño para absorber redondeos de vértices)
buffer_m <- 10  # puedes ajustar (0, 5, 10)
borde_sf_pad <- st_buffer(borde_sf, buffer_m)

# (b) Quedarnos solo con puntos dentro del borde "acolchado"
pts_sf_in <- st_filter(pts_sf, borde_sf_pad, .pred = st_within)

cat("Puntos totales:", nrow(pts_sf), 
    "| dentro (con buffer", buffer_m, "m):", nrow(pts_sf_in), "\n")

# (c) Duplicados: si hay colisiones exactamente en las mismas coords, jitter leve (0.5 m)
xy_in <- st_coordinates(pts_sf_in)
dup_mask <- duplicated(data.frame(x = xy_in[,1], y = xy_in[,2]))
n_dup <- sum(dup_mask)
cat("Puntos duplicados (coincidencia exacta):", n_dup, "\n")

if (n_dup > 0) {
  pts_sf_in <- st_jitter(pts_sf_in, amount = 0.5)  # jitter de 0.5 m
}

# =========================================================
# 5) CONVERSIÓN A spatstat (as.owin desde sf) y construcción del ppp
# =========================================================
win_owin <- as.owin(borde_sf)            # ventana con el borde real (sin buffer)
xy_ok    <- st_coordinates(pts_sf_in)

acc_ppp <- ppp(
  x = xy_ok[,1],
  y = xy_ok[,2],
  window = win_owin,
  marks  = data.frame(sexo = pts_sf_in$sexo, tipo = pts_sf_in$tipo)
)

print(acc_ppp)
summary(marks(acc_ppp))

# =========================================================
# 6) EDA ESPACIAL GENERAL
# =========================================================

# 6.1) Mapa general
par(mar = c(0,0,0,0))
plot(win_owin, main = "", col = "grey90", border = 1)
plot(acc_ppp, add = TRUE, cex = 0.5, pch = 19)

# 6.2) Conteo por cuadrantes + prueba CSR
par(mfrow = c(1,2), mar = c(2,2,1,1))
q_cnt <- quadratcount(acc_ppp, nx = 4, ny = 4)
plot(q_cnt, main = "Frecuencias por cuadrante")
plot(acc_ppp, add = TRUE, pch = 3, cex = 0.4)
qt_res <- quadrat.test(acc_ppp, nx = 4, ny = 4)
plot(qt_res, main = "Residuos de Pearson")

cat("\n--- Resultado quadrat.test ---\n")
print(qt_res)

# 6.3) Densidad kernel
par(mfrow = c(1,1), mar = c(2,2,1,1))
den <- density(acc_ppp, sigma = bw.diggle)
image(den, main = "Densidad kernel (bw.diggle)")
contour(den, add = TRUE)
plot(win_owin, add = TRUE, border = 1)

# 6.4) Función K
k_acc <- Kest(acc_ppp)
plot(k_acc, main = "Kest: Accidentes fatales (general)")

#####
# =========================================================
# 1) Resumen descriptivo de variables
# =========================================================

# Conteo por sexo
table_sexo <- table(marks(acc_ppp)$sexo)
prop_sexo  <- prop.table(table_sexo)

cat("\n--- Distribución por sexo ---\n")
print(table_sexo)
print(round(100*prop_sexo, 1))

# Conteo por tipo de automotor
table_tipo <- table(marks(acc_ppp)$tipo)
prop_tipo  <- prop.table(table_tipo)

cat("\n--- Distribución por tipo de automotor ---\n")
print(table_tipo)
print(round(100*prop_tipo, 1))

# Tabla cruzada Sexo x Tipo
cat("\n--- Tabla cruzada Sexo x Tipo ---\n")
print(table(marks(acc_ppp)$sexo, marks(acc_ppp)$tipo))

# =========================================================
# 2) EDA por SEXO
# =========================================================

# Dividir por sexo
ppp_by_sexo <- split(acc_ppp, f = marks(acc_ppp)$sexo)

par(mfrow = c(1,2), mar = c(0,0,1,0))
for (s in names(ppp_by_sexo)) {
  plot(win_owin, main = paste("Accidentes - Sexo:", s), col = "grey90")
  plot(ppp_by_sexo[[s]], add = TRUE, pch = 19, cex = 0.5)
}

# Densidades por sexo
par(mfrow = c(1,2), mar = c(2,2,1,1))
for (s in names(ppp_by_sexo)) {
  den_s <- density(ppp_by_sexo[[s]], sigma = bw.diggle)
  image(den_s, main = paste("Densidad kernel - Sexo:", s))
  contour(den_s, add = TRUE)
  plot(win_owin, add = TRUE)
}

# =========================================================
# 1) Función K por SEXO
# =========================================================

ppp_by_sexo <- split(acc_ppp, f = marks(acc_ppp)$sexo)

par(mfrow = c(1,1), mar = c(4,4,2,1))
plot(Kest(ppp_by_sexo$F), main = "Función K por sexo", col = "red")
plot(Kest(ppp_by_sexo$M), add = TRUE, col = "blue")
legend("topleft", legend = c("Femenino","Masculino"),
       col = c("red","blue"), lty = 1)

# =========================================================
# 2) Test / mapa de riesgo relativo por SEXO (M vs F)
# =========================================================

# Construir un ppp con UNA sola marca (sexo)
ppp_sexo <- acc_ppp
marks(ppp_sexo) <- factor(marks(acc_ppp)$sexo, levels = c("F","M"))

# Bandwidth recomendado para relrisk
sig_sexo <- bw.relrisk(ppp_sexo)

# Riesgo relativo espacial: Hombres vs Mujeres
sex_rr <- relrisk(ppp_sexo, casecontrol = TRUE,
                  case = "M", control = "F",
                  sigma = sig_sexo)

par(mfrow = c(1,1), mar = c(2,2,2,4))
plot(sex_rr, main = "Riesgo relativo espacial (Hombres vs Mujeres)")
plot(win_owin, add = TRUE)

# =========================================================
# 3) Función K por TIPO DE AUTOMOTOR (principales)
# =========================================================

tipos_top <- c("PEATON","COND. MOTO","CICLISTA")
ppp_by_tipo <- split(acc_ppp, f = marks(acc_ppp)$tipo)

par(mfrow = c(1,1), mar = c(4,4,2,1))
cols <- c("darkgreen","purple","orange")
i <- 1
for (t in tipos_top) {
  K_t <- Kest(ppp_by_tipo[[t]])
  if (i == 1) {
    plot(K_t, main = "Función K por tipo de automotor", col = cols[i])
  } else {
    plot(K_t, add = TRUE, col = cols[i])
  }
  i <- i + 1
}
legend("topleft", legend = tipos_top, col = cols, lty = 1)

# =========================================================
# 4) Riesgo relativo espacial por TIPO (PEATÓN vs OTROS)
# =========================================================

# Construir un ppp con UNA sola marca binaria (PEATON / OTROS)
ppp_tipo2 <- acc_ppp
marks(ppp_tipo2) <- factor(ifelse(marks(acc_ppp)$tipo == "PEATON", "PEATON", "OTROS"),
                           levels = c("OTROS","PEATON"))

sig_tipo <- bw.relrisk(ppp_tipo2)

risk_peaton <- relrisk(ppp_tipo2, casecontrol = TRUE,
                       case = "PEATON", control = "OTROS",
                       sigma = sig_tipo)

par(mar = c(2,2,2,4))
plot(risk_peaton, main = "Riesgo relativo espacial: PEATÓN vs OTROS")
plot(win_owin, add = TRUE)

