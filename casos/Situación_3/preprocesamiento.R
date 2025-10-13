#----------------------------------------------------
# LIMPIEZA DE DATOS PROBLEMÁTICA 3: CASOS DE DENGUE #
#----------------------------------------------------
install.packages("clipr")

library(sf)           
library(readxl)       
library(dplyr)        
library(stringi)      
library(stringr)      
library(tidyr)        
library(spdep) 
library(clipr)

# 1. Rutas y lectura de archivos ----------------------------------------------
setwd("~/project_2/datasets/datasets/Situación_3")

shp_candidates <- c("Cartografia/barrios_mo.shp", "./Cartografia/barrios_mo.shp", "/mnt/data/barrios_mo.shp", "barrios_mo.shp")
shp_path <- shp_candidates[file.exists(shp_candidates)][1]
if(is.na(shp_path) || is.null(shp_path)) stop("No se encontró el shapefile 'barrios_mo.shp' en las rutas esperadas. Ajusta shp_path.")

csv_path <- "DatosDengueCaliBarrios.csv"
if(!file.exists(csv_path)) stop("No se encontró DatosDengueCaliBarrios.csv en /mnt/data")

xls_path <- "InfomacionCaliCenso2005Barrio.xls"
if(!file.exists(xls_path)) stop("No se encontró InfomacionCaliCenso2005Barrio.xls en /mnt/data")

shp <- st_read(shp_path, quiet = TRUE)
plot(st_geometry(shp), axes = T)

dengue <- read.csv(csv_path, fileEncoding = "LATIN1")

# Leer Excel: detectamos la fila de encabezado y leemos desde ahí
censo_raw <- read_excel("InfomacionCaliCenso2005Barrio.xls", 
                        skip = 7)

names(censo_raw)

censo <- censo_raw %>%
  rename(Codigo = ...1,
         Barrio = ...2,
         Poblacion_Total = Total,
         Poblacion_H = H,
         Poblacion_M = M)

censo <- censo %>%
  filter(!is.na(Barrio)) %>%
  filter(!str_detect(Barrio, regex("Total|Comuna", ignore_case = TRUE)))

censo <- censo %>%
  mutate(across(starts_with("Poblacion"),
                ~ as.numeric(str_replace_all(., "[.]", ""))))


# 2. Estandarización de nombres y limpieza -----------------------------------
# Normalizamos nombres:

clean_name <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- stri_trans_general(x, "Latin-ASCII")
  x <- gsub("[^a-z0-9 ]+", " ", x)
  x <- gsub("\\s+", " ", trimws(x))
  x
}

# CSV dengue: primera columna con nombre extraño, por lo que se renombra a BARRIO
if(!"BARRIO" %in% names(dengue)){
  names(dengue)[1] <- "BARRIO"
}

# Asegurar columnas numéricas
dengue$CDENGUE <- as.numeric(dengue$CDENGUE)
dengue$NSUMIDERO <- as.numeric(dengue$NSUMIDERO)

dengue$BARRIO_clean <- clean_name(dengue$BARRIO)
censo$Barrio_clean <- clean_name(censo$Barrio)

# Shapefile: detectar columna que contiene nombre del barrio
possible_shp_names <- names(shp)[grepl("barrio|name|nom|label|NOMBRE|NAME", names(shp), ignore.case = TRUE)]
if(length(possible_shp_names) >= 1){
  shp_name_col <- possible_shp_names[1]
} else {
  # fallback: usar ID si existe
  shp_name_col <- names(shp)[1]
  message("No se detectó campo obvio de nombre en shapefile; se usará: ", shp_name_col)
}

shp$BARRIO_shape <- as.character(shp[[shp_name_col]])
shp$BARRIO_shape_clean <- clean_name(shp$BARRIO)

#Otras transformaciones para nombre de barrio por mapeo

diccionario <- data.frame(
  dengue = c(
    "ciudad de los alamos",
    "petecuy ii",
    "chiminangos ii",
    "olaya herrera",
    "chiminangos i",
    "petecuy iii",
    "petecuy i",
    "alfonso lopez i",
    "las delicas",
    "alfonso lopez ii",
    "alfonso lopez iii",
    "obrero",
    "villacolombia",
    "jose manuel marroquin ii",
    "alfonso barberena a",
    "los naranjos",
    "los comuneros ii",
    "eucaristico",
    "jose manuel marroquin i",
    "san carlos",
    "departamental",
    "los comuneros i",
    "u d a galindo plaza de toros",
    "camino real joaquin borrero sinisterra",
    "caldas",
    "cuarteles de napoles"
  ),
  censo = c(
    "ciudad los alamos",
    "petecuy segunda etapa",
    "chiminangos segunda etapa",
    "barrio olaya herrera",
    "chiminangos primera etapa",
    "petecuy tercera etapa",
    "petecuy primera etapa",
    "alfonso lopez 1a etapa",
    "las delicias",
    "alfonso lopez 2a etapa",
    "alfonso lopez 3a etapa",
    "barrio obrero",
    "villa colombia",
    "jose manuel marroquin ii etapa",
    "barrio alfonso barberena a",
    "los naranjos i",
    "los comuneros ii etapa",
    "barrio eucaristico",
    "jose manuel marroquin i etapa",
    "barrio san carlos",
    "barrio departamental",
    "los comuneros i etapa",
    "unidad deportiva a galindo plaza de toros",
    "camino real j borrero s",
    "barrio caldas",
    "cuarteles napoles"
  )
)

censo <- censo %>%
  left_join(diccionario, by = c("Barrio_clean" = "censo")) %>%
  mutate(
    Barrio_clean = ifelse(!is.na(dengue), dengue, Barrio_clean)
  )

# Si la columna "dengue" existe, la eliminamos después del mapeo
if("dengue" %in% names(censo)) {
  censo <- censo %>% dplyr::select(-dengue)
}

# 3. Unión de tablas ---------------------------------------------------------

#Se comprueba qué barrios no están en dengue que sí estén en censo:
barrios_faltantes <- setdiff(censo$Barrio_clean, dengue$BARRIO_clean)
barrios_faltantes
length(barrios_faltantes)

#Y caso contrario:
setdiff(dengue$BARRIO_clean, censo$Barrio_clean) #aparece el NA

# Unir dengue y censo por nombres limpios
dengue_censo <- dengue %>%
  dplyr::select(BARRIO_clean, CDENGUE, NSUMIDERO) %>%
  left_join(
    censo %>% dplyr::select(Barrio_clean, Poblacion_Total, Poblacion_H, Poblacion_M),
    by = c("BARRIO_clean" = "Barrio_clean")
  )

#Al unir se pierden 12 barrios que no existen en la tabla de casos de dengue y tampoco en el shapefile 

# Unir con shapefile por nombres limpios
shp2 <- shp %>%
  left_join(dengue_censo, 
            by = c("BARRIO_shape_clean" = "BARRIO_clean"))



# 4. ANÁLISIS EXPLORATORIO DE DATOS -----------------------
na_resumen <- shp2 %>%
  st_drop_geometry() %>%
  summarise(across(everything(), ~ sum(is.na(.)), .names = "na_{col}")) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "num_NA") %>%
  arrange(desc(num_NA))

na_resumen

na_barrios <- shp2 %>%
  st_drop_geometry() %>%
  dplyr::filter(is.na(CDENGUE) | is.na(NSUMIDERO) | is.na(Poblacion_Total)) %>%
  dplyr::select(BARRIO_shape_clean, CDENGUE, NSUMIDERO, Poblacion_Total)

na_barrios

#Confirmamos que hay una fila con valores nulos tanto en el nombre del barrio como con población. 
#Ahora pasamos a revisar si hay barrios con población total igual a 0

barrios_pob0 <- shp2 %>%
  st_drop_geometry() %>%
  dplyr::filter(Poblacion_Total == 0) %>%
  dplyr::select(BARRIO_shape_clean, CDENGUE, NSUMIDERO, Poblacion_Total)

barrios_pob0

#Aquí tenemos como resultado tres barrios con población igual a 0, que son industria de licores, Unidad Deportiva
#galindo plaza de toros y club campestre. No se encontró otro censo donde haya registro al menos de las poblaciones totales
#por lo que se opta por eliminarlas junto con el dato nulo. También, se opta por eliminarlos al no ser barrios como tal

shp2 <- shp2 %>%
  dplyr::filter(!is.na(Poblacion_Total) & Poblacion_Total > 0)

#Y confirmamos que haya funcionado con:
summary(shp2$Poblacion_Total)
sum(is.na(shp2$Poblacion_Total))
sum(shp2$Poblacion_Total == 0, na.rm = TRUE)
nrow(shp2)

#Ahora, revisaremos si hay barrios donde hayan más casos de dengue que población total:

barrios_inconsistentes <- shp2 %>%
  st_drop_geometry() %>%
  dplyr::filter(CDENGUE > Poblacion_Total) %>%
  dplyr::select(BARRIO_shape_clean, CDENGUE, Poblacion_Total, NSUMIDERO)

barrios_inconsistentes

#Se confirma que hay 15 barrios con esta condición. Los tratamos como datos inconsistentes,
# así que los volvemos NA y los imputamos usando el promedio de los barrios vecinos, usando 
#la vecindad tipo queen.

shp2 <- shp2 %>%
  mutate(CDENGUE = ifelse(CDENGUE > Poblacion_Total, NA, CDENGUE))

nb <- poly2nb(shp2)
listw <- nb2listw(nb, style="W", zero.policy=TRUE)

barrios_imputar <- which(is.na(shp2$CDENGUE))

#Aquí pasamos a aplicar el promedio de los barrios vecinos que tengan valores válidos, en caso de que no 
#Se usa el promedio global.

for (i in barrios_imputar) {
  vecinos <- nb[[i]]
  vecinos_validos <- vecinos[!is.na(shp2$CDENGUE[vecinos])]

  if (length(vecinos_validos) > 0) {
    shp2$CDENGUE[i] <- mean(shp2$CDENGUE[vecinos_validos], na.rm=TRUE)
  } else {
    # Promedio global
    shp2$CDENGUE[i] <- mean(shp2$CDENGUE, na.rm=TRUE)
  }
}

# Aquí nos damos cuenta que usar el promedio con vecinos cercanos no soluciona el problema.
# Usaremos la mediana de la tasa de dengue para poder imputar:

shp2 <- shp2 %>%
  mutate(tasa_dengue = CDENGUE / Poblacion_Total * 10000)

shp2$CDENGUE[shp2$CDENGUE > shp2$Poblacion_Total] <- NA

# Calcular la tasa de dengue solo para barrios válidos
shp2 <- shp2 %>%
  mutate(tasa_dengue_valid = ifelse(Poblacion_Total > 0 & !is.na(CDENGUE), CDENGUE / Poblacion_Total, NA))

# Se calcula la mediana de la tasa usando solo los casos válidos
mediana_tasa <- median(shp2$tasa_dengue_valid, na.rm = TRUE)

#Se imputa CDENGUE en los barrios inválidos con a mediana de la tasa

shp2 <- shp2 %>%
  mutate(
    CDENGUE = ifelse(
      is.na(CDENGUE),
      round(pmin(Poblacion_Total * mediana_tasa, Poblacion_Total)),  # Se supone que nunca habrá más casos que personas
      CDENGUE
    )
  )

shp2 <- shp2 %>% dplyr::select(-tasa_dengue_valid)
shp2 <- shp2 %>% dplyr::select(-tasa_dengue)

#Ahora sí calculamos la columna de tasa

shp2 <- shp2 %>%
  mutate(tasa_dengue = CDENGUE / Poblacion_Total * 10000)


st_write(shp2, "datos_dengue.shp", delete_dsn = TRUE)

#Finalmente el shp queda con 334 objetos. Ahora pasamos a crear la tasa de casos de dengue 
#por cada 10000 habitantes:

shp2 <- shp2 %>%
  mutate(
    tasa_10000 = (CDENGUE / Poblacion_Total) * 10000)

