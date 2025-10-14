
## SITUACIÓN 2: Homicidios por Accidentes de tránsito en Cali
#################################################

# =========================================================
# Librerías
# =========================================================
# install.packages(c("readxl","sf","sp","spatstat","spatstat.model","spdep","grDevices","grid","stars"))
# install.packages("dplyr")

library(readxl)
library(sf)
library(sp)
library(spdep)
library(spatstat)
library(spatstat.model)
library(grDevices); library(grid)
library(stars)
library(readr)
library(dplyr)

# =========================================================
# 1) CARGA DE DATOS DESDE GITHUB
# =========================================================
csv_url <- "https://raw.githubusercontent.com/DCajiao/Aplicaciones-del-Analisis-Espacial/refs/heads/main/casos/Situación_2/data/output/Accidentes_2009_2010.csv"
map_dir <- "C:/Users/david/Desktop/Proyecto 2/Aplicaciones-del-Analisis-Espacial/casos/Situación_2/data/Mapas"

raw <- read_csv(csv_url, show_col_types = FALSE)

# Columnas de interés + AÑO y TIPO_AUTOMOTOR explícitos 
# (Estas columnas ya están en el CSV, fueropn procesadas desde el notebook de Python)
df <- data.frame(
  x_m  = as.numeric(raw[["coordenada X (metros)"]]),
  y_m  = as.numeric(raw[["coordenada Y (metros)"]]),
  SEXO = as.factor(raw[["SEXO"]]),
  EDAD_AGRUPADA = as.factor(raw[["EDAD AGRUPADA"]]),
  TIPO_AUTOMOTOR = as.factor(raw[["TIPO_AUTOMOTOR"]]),
  ANIO = as.integer(raw[["AÑO_DATA"]]),
  BARRIO = raw[["BARRIO"]],
  COM = raw[["COM"]]
)

# Limpieza mínima: quitar filas sin coordenadas válidas ni año
# Sólo para estar seguros, de igual manera el csv ya viene limpio
df <- df |>
  filter(is.finite(x_m), is.finite(y_m), !is.na(ANIO))

# =========================================================
# 2) CAPAS GEOGRÁFICAS Y CRS (EPSG:3115)
# =========================================================
# Lee el borde municipal y barrios (forzando EPSG:3115)
# Lectura del borde municipal y barrios (FORZANDO EPSG:3115 para compatibilidad)
borde <- st_read(file.path(map_dir, "BordeComunasMetros.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(3115)

barrios <- st_read(file.path(map_dir, "Barrios.shp"), quiet = TRUE) |>
  st_make_valid() |>
  st_transform(3115)

# Puntos como sf en 3115
pts_sf <- st_as_sf(df, coords = c("x_m","y_m"), crs = 3115)

# Asegurar que los puntos están dentro del borde
pts_sf <- pts_sf[st_within(pts_sf, st_union(borde), sparse = FALSE), ]

# =========================================================
# 3) CONVERTIR A OBJETOS SPATSTAT (owin + ppp)
# =========================================================

# --- Limpieza robusta del polígono ---
borde_clean <- borde |>
  st_zm(drop = TRUE, what = "ZM") |>   # quita dimensiones extra
  st_collection_extract("POLYGON") |>  # asegura geometrías planas
  st_make_valid() |>                   # repara polígonos inválidos
  st_transform(3115)                   # fuerza CRS EPSG:3115

# --- Unir a un solo polígono ---
borde_union <- st_union(borde_clean)

# --- Convertir a ventana owin ---
borde_win <- spatstat.geom::as.owin(borde_union)
plot(borde_win, main = "Ventana de estudio (Cali)")

# =========================================================
# 4) CONSTRUIR PATRÓN DE PUNTOS (PPP)
# =========================================================

# Extraer coordenadas de los puntos
coords <- st_coordinates(pts_sf)

# DataFrame de marcas
marks_df <- pts_sf |>
  st_drop_geometry() |>
  dplyr::select(SEXO, EDAD_AGRUPADA, TIPO_AUTOMOTOR, ANIO)

# Crear objeto ppp
X_all <- spatstat.geom::ppp(
  x = coords[,1],
  y = coords[,2],
  window = borde_win,
  marks = marks_df
)

borde_sp <- as(borde_union, "Spatial")


# Print y plot para probar
# print(X_all)
# plot(X_all, main = "Eventos de homicidios viales (2009–2010)")

# =========================================================
# 4)  FUNCIONES PARA  ESDA
# =========================================================

pal_fun <- function(n) grDevices::hcl.colors(n, palette = "Dark2")

plot_points_by <- function(X, by, main, pch_pt=16, cex_pt=0.55) {
  stopifnot(by %in% colnames(marks(X)))
  v   <- marks(X)[[by]]; if (!is.factor(v)) v <- factor(v)
  lev <- levels(v)
  pal <- pal_fun(length(lev))
  cols <- pal[as.integer(v)]
  
  if (missing(main)) main <- sprintf("Eventos por %s", by)
  plot(borde_sp, main=main)
  points(X$x, X$y, pch=pch_pt, cex=cex_pt, col=cols)   # <<-- COLORES
  legend("topright", legend=lev, pch=16, pt.cex=0.9, col=pal, bty="n", title=by)
}

plot_quadrants <- function(X, nx=6, ny=6, main="Conteo por cuadrantes") {
  plot(borde_sp, main=main)
  points(X$x, X$y, pch=".", cex=0.9, col="#00000080")
  plot(quadratcount(X, nx=nx, ny=ny), add=TRUE, col="blue")
}

run_quadrat_test <- function(X, nx=6, ny=6) {
  if (npoints(X) >= 10) print(quadrat.test(X, nx=nx, ny=ny))
  else message(sprintf("quadrat.test omitido (n=%d < 10)", npoints(X)))
}

plot_kernel <- function(X, main="Densidad Kernel", sigma=NULL) {
  if (is.null(sigma)) sigma <- bw.diggle(X)
  den <- density(X, sigma=sigma)
  image(den, main=main)
  plot(borde_sp, add=TRUE, border=1)
  contour(den, add=TRUE, col=1)
}

subset_year <- function(X, year) {
  idx <- with(marks(X), ANIO == year)
  X[idx]
}

# Esta carpeta espara outputs del ESDA
dir.create("outputs", showWarnings = FALSE)

esda_estructura_anual <- function(year, nx=6, ny=6) {
  stopifnot(year %in% unique(marks(X_all)$ANIO))
  Xy <- subset_year(X_all, year)
  sig <- bw.diggle(Xy)
  
  message(sprintf("== ESDA estructural | Año %d | n=%d ==", year, npoints(Xy)))
  
  png(file.path("outputs", sprintf("puntos_SEXO_%d.png", year)), 1400, 1000, res=140)
  plot_points_by(Xy, by="SEXO", main=sprintf("Eventos por SEXO | %d", year))
  dev.off()
  
  png(file.path("outputs", sprintf("puntos_EDAD_%d.png", year)), 1400, 1000, res=140)
  plot_points_by(Xy, by="EDAD_AGRUPADA", main=sprintf("Eventos por EDAD | %d", year))
  dev.off()
  
  png(file.path("outputs", sprintf("puntos_TIPO_%d.png", year)), 1400, 1000, res=140)
  plot_points_by(Xy, by="TIPO_AUTOMOTOR", main=sprintf("Eventos por TIPO_AUTOMOTOR | %d", year))
  dev.off()
  
  png(file.path("outputs", sprintf("cuadrantes_%d.png", year)), 1400, 1000, res=140)
  plot_quadrants(Xy, nx=nx, ny=ny, main=sprintf("Cuadrantes | %d", year))
  dev.off()
  run_quadrat_test(Xy, nx=nx, ny=ny)
  
  png(file.path("outputs", sprintf("kernel_%d.png", year)), 1400, 1000, res=140)
  plot_kernel(Xy, main=sprintf("Densidad Kernel | %d", year), sigma=sig)
  dev.off()
  
  invisible(Xy)
}

# Uso de las funciones para el esda estructural por año
X_2009 <- esda_estructura_anual(2009, nx=6, ny=6)
X_2010 <- esda_estructura_anual(2010, nx=6, ny=6)


# ==========================================================
# 5) COMPARACIÓN DE MODELOS POR AÑO + MÉTRICAS DE ERROR
# ==========================================================
# Esta carpeta espara outputs del modelado
dir.create("outputs_modelo", showWarnings = FALSE)

# ---------- Intensidad robusta (Poisson/Gibbs) ----------
pal_magma <- colorRampPalette(c("#0d0887","#6a00a8","#b12a90",
                                "#e16462","#fca636","#f0f921"))


compute_intensity_im <- function(fit, Xy, dimyx = 256,
                                 nsim = 100,         # más simulaciones = más liso
                                 sigma_mult = 1.8,   # más suave que bw.diggle
                                 seed = 123) {
  # 1) Intento intensidad total
  im1 <- try(predict(fit, type = "intensity", dimyx = dimyx), silent = TRUE)
  if (!inherits(im1, "try-error") && any(is.finite(as.vector(im1$v)))) return(im1)
  
  # 2) Intento CIF (útil en Gibbs)
  im2 <- try(predict(fit, type = "cif", dimyx = dimyx), silent = TRUE)
  if (!inherits(im2, "try-error") && any(is.finite(as.vector(im2$v)))) return(im2)
  
  # 3) Fallback: simulación + promedio de densidades
  set.seed(seed)
  sims <- simulate(fit, nsim = nsim, progress = FALSE)
  if (!is.list(sims)) sims <- list(sims)
  
  sig0 <- bw.diggle(Xy)
  sig  <- as.numeric(sig0) * sigma_mult
  ims  <- lapply(sims, function(si) density(si, sigma = sig, dimyx = dimyx, at = "pixels"))
  num  <- Reduce(`+`, ims)
  eval.im(num / length(ims))
}

# ---------- utilidades de métrica (robusto + simulación si hace falta) ----------
metrics_from_model <- function(fit, Xy, nx = 6, ny = 6, dimyx = 256, nsim = 60, seed = 123) {
  Q <- quadrats(Window(Xy), nx = nx, ny = ny)
  obs <- as.integer(quadratcount(Xy, tess = Q))
  
  # Esperados vía intensidad; si falla, simula y promedia conteos
  # (Esto lo planteamos de esa manera para que funcione con Gibbs)
  expct <- try({
    lam <- predict(fit, type = "intensity", dimyx = dimyx)
    Ti  <- tiles(Q)
    sapply(Ti, function(Wi) integral.im(lam, Wi))
  }, silent = TRUE)
  
  if (inherits(expct, "try-error") || any(!is.finite(expct))) {
    set.seed(seed)
    sims <- simulate(fit, nsim = nsim, progress = FALSE)
    if (!is.list(sims)) sims <- list(sims)
    sim_counts <- sapply(sims, function(si) as.integer(quadratcount(si, tess = Q)))
    expct <- rowMeans(sim_counts, na.rm = TRUE)
  }
  
  keep <- is.finite(expct)
  obs  <- obs[keep]
  expc <- expct[keep]
  
  RMSE <- sqrt(mean((obs - expc)^2))
  MAE  <- mean(abs(obs - expc))
  MAPE <- mean(ifelse(obs > 0, abs(obs - expc) / obs, NA), na.rm = TRUE) * 100
  
  list(RMSE = RMSE, MAE = MAE, MAPE = MAPE, obs = obs, exp = expc)
}

# ---------- mapa de intensidad (1 panel limpio) ----------
plot_intensity_png <- function(fit, Xy, year, label, file,
                               dimyx = 256, nsim = 100, sigma_mult = 1.8) {
  png(file.path("outputs_modelo", file), width = 1600, height = 1200, res = 160)
  par(mfrow = c(1,1), mar = c(4,4,3,7))
  inten <- compute_intensity_im(fit, Xy, dimyx = dimyx, nsim = nsim, sigma_mult = sigma_mult)
  inten_sqrt <- eval.im(sqrt(pmax(inten, 0)))
  plot(inten_sqrt, col = pal_magma(64),
       main = sprintf("Intensidad (modelo: %s) | %d", label, year),
       ribargs = list(las = 1, cex.axis = 0.9))
  plot(borde_win, add = TRUE, border = 1, lwd = 0.8)
  dev.off()
}

# ---------- mapa de residuos (1 panel limpio) ----------
plot_residuals_png <- function(fit, year, file) {
  png(file.path("outputs_modelo", file), width = 1600, height = 1200, res = 160)
  par(mfrow = c(1,1), mar = c(4,4,3,7))
  res_sm <- Smooth(residuals(fit, type = "pearson"))
  pal_div <- colorRampPalette(c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"))
  plot(res_sm, col = pal_div(64),
       main = sprintf("Residuos de Pearson (suavizados) | %d", year),
       ribargs = list(las = 1, cex.axis = 0.9))
  plot(borde_win, add = TRUE, border = 1, lwd = 0.8)
  contour(res_sm, add = TRUE, drawlabels = FALSE)
  dev.off()
}

# ---------- función principal por año ----------
comparar_modelos_anio <- function(year, nx = 6, ny = 6,
                                  poly_deg = 2, r_strauss = 200,
                                  incluir_covar = FALSE,
                                  plot_all = FALSE) {
  stopifnot(year %in% unique(marks(X_all)$ANIO))
  Xy <- unmark(X_all[with(marks(X_all), ANIO == year)])
  
  message(sprintf("== Comparación de modelos | %d | n=%d ==", year, npoints(Xy)))
  
  modelos   <- list()
  etiquetas <- list()
  
  # (1) CSR
  modelos$csr <- ppm(Xy ~ 1); etiquetas$csr <- "CSR"
  
  # (2) Poisson no estacionario (tendencia)
  f_trend <- if (poly_deg <= 1) ~ x + y else ~ x + y + I(x^2) + I(y^2) + I(x*y)
  modelos$inhom <- ppm(Xy, trend = f_trend); etiquetas$inhom <- "Inhom(x,y)"
  
  # (3) Interacción (Strauss)
  modelos$strauss <- ppm(Xy, trend = ~ 1, interaction = Strauss(r = r_strauss))
  etiquetas$strauss <- sprintf("Strauss(r=%d)", r_strauss)
  
  # (4) Opcional: covariable (distancia a centroide)
  if (incluir_covar) {
    cxy <- as.numeric(st_coordinates(st_centroid(st_as_sf(borde_union))))
    dist_im <- as.im(function(x,y) sqrt((x - cxy[1])^2 + (y - cxy[2])^2), W = borde_win)
    modelos$covar <- ppm(Xy ~ dist_im + I(dist_im^2))
    etiquetas$covar <- "Covar(dist_centro)"
  }
  
  # --- Tabla AIC, pR2, BIC ---
  ll_csr <- as.numeric(logLik(modelos$csr))  # (pseudo)logLik para referencia
  aic_tab <- data.frame(
    Modelo   = names(modelos),
    Etiqueta = unlist(etiquetas),
    AIC      = sapply(modelos, AIC),
    logLik   = sapply(modelos, function(m) as.numeric(logLik(m)))
  )
  # pR2 de McFadden vs CSR (para Gibbs usa pseudo-LL, válido para comparación relativa)
  aic_tab$pR2 <- 1 - (aic_tab$logLik / ll_csr)
  aic_tab$pR2[aic_tab$Modelo == "csr"] <- 0
  
  # BIC = -2*logLik + k*log(n)  (k = #parámetros; n = #puntos)
  get_k <- function(m) {
    k <- attr(logLik(m), "df")
    if (is.null(k)) k <- length(coef(m))
    as.numeric(k)
  }
  n_pts <- npoints(Xy)
  aic_tab$BIC <- (-2 * aic_tab$logLik) + sapply(modelos, get_k) * log(n_pts)
  
  aic_tab <- aic_tab[order(aic_tab$AIC), ]
  print(aic_tab)
  
  # --- Métricas RECM/MAE/MAPE ---
  mets <- lapply(modelos, function(m) metrics_from_model(m, Xy, nx = nx, ny = ny))
  err_tab <- data.frame(
    Modelo   = names(modelos),
    Etiqueta = unlist(etiquetas),
    RECM     = sapply(mets, `[[`, "RMSE"),
    EAM      = sapply(mets, `[[`, "MAE"),
    MAPE     = sapply(mets, `[[`, "MAPE")
  )
  rownames(err_tab) <- NULL
  err_tab <- err_tab[match(aic_tab$Modelo, err_tab$Modelo), ]
  print(err_tab)
  
  # --- Elegir mejor por AIC y graficar 2 figuras grandes ---
  best_name  <- aic_tab$Modelo[1]
  best_fit   <- modelos[[best_name]]
  best_label <- aic_tab$Etiqueta[1]
  
  if (!plot_all) {
    # Solo mejor modelo
    plot_intensity_png(best_fit, Xy, year, best_label, sprintf("ppm_trend_%d.png", year))
    plot_residuals_png(best_fit, year, sprintf("ppm_residuos_%d.png", year))
  } else {
    # Guardar TODOS los modelos
    for (nm in names(modelos)) {
      fit_i   <- modelos[[nm]]
      label_i <- etiquetas[[nm]]
      suf     <- gsub("[^A-Za-z0-9]+", "_", label_i)
      plot_intensity_png(fit_i, Xy, year, label_i,
                         sprintf("ppm_trend_%s_%d.png", suf, year))
      plot_residuals_png(fit_i, year,
                         sprintf("ppm_residuos_%s_%d.png", suf, year))
    }
  }
}

# -----------------------
# Ejecutar la funcion de comparacion de modelos por año para 2009 y 2010 y guardar resultados 
# en la carpeta outputs_modelo
# -----------------------
res_2009 <- comparar_modelos_anio(2009, nx = 6, ny = 6, poly_deg = 2, r_strauss = 200, incluir_covar = TRUE)
res_2010 <- comparar_modelos_anio(2010, nx = 6, ny = 6, poly_deg = 2, r_strauss = 200, incluir_covar = TRUE)
