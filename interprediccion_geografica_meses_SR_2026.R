# ===================================================================
# Geographic interprediction calculation between seasonal niches and migratory months niches: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
# ===================================================================

##This code works based on authors file names and inner file structure. To run it with other file names and structure, adjust code

library(raster)
library(sf)
library(dplyr)
library(purrr)
library(lubridate)

# -------------------------------
# Raster short names
# -------------------------------
nombres_cortos <- list(
  invierno = "invierno",
  reproduccion = "repro"
)

# Migratory season of each month
mes_a_temporada <- list(
  febrero = "migprim",
  marzo = "migprim",
  abril = "migprim",
  julio = "migfall",
  agosto = "migfall",
  septiembre = "migfall",
  octubre = "migfall"
)

# -------------------------------
# Main function
# -------------------------------

evaluar_interprediccion <- function(temporada_modelo, 
                                    mes_presencia, 
                                    dir_rasters, 
                                    dir_presencias, 
                                    dir_shapes = NULL, 
                                    usar_shape = FALSE) {
  
  modelo_corto <- nombres_cortos[[temporada_modelo]]
  raster_path <- file.path(dir_rasters, paste0(modelo_corto, "_to_", mes_presencia, "_bin.tif"))
  
  # Determine the source season of the month
  temporada_migratoria <- mes_a_temporada[[mes_presencia]]
  presencia_path <- file.path(dir_presencias, paste0("data_", temporada_migratoria, ".csv"))
  
  if (!file.exists(raster_path)) {
    warning(paste("Raster not found:", raster_path))
    return(NULL)
  }
  if (!file.exists(presencia_path)) {
    warning(paste("Records file not found:", presencia_path))
    return(NULL)
  }
  
  # Read and filter data per month 
  data <- read.csv(presencia_path)
  data$fecha <- as.Date(data$layer_dates)
  data_filtrada <- data[as.character(month(data$fecha, label = TRUE, abbr = FALSE)) == mes_presencia, ]
  if (nrow(data_filtrada) == 0) {
    warning(paste("No data for the month found:", mes_presencia))
    return(data.frame(
      modelo = temporada_modelo,
      mes = mes_presencia,
      n_total = 0,
      n_predicho = NA,
      porcentaje_predicho = NA
    ))
  }
  
  r <- raster(raster_path)
  pres_sf <- st_as_sf(data_filtrada, coords = c("lon.1", "lat.1"), crs = 4326)
  pres_sf <- st_transform(pres_sf, crs = crs(r))
  
  # Shapefile cutting (optional)
  if (usar_shape && !is.null(dir_shapes)) {
    shp_path <- file.path(dir_shapes, paste0("m_",temporada_migratoria, ".shp"))
    if (file.exists(shp_path)) {
      shp <- st_read(shp_path, quiet = TRUE) %>% st_transform(crs = crs(r))
      r <- raster::mask(raster::crop(r, shp), shp)
      pres_sf <- pres_sf[st_within(pres_sf, shp, sparse = FALSE), ]
    } else {
      warning(paste("Shapefile not found:", shp_path))
    }
  }
  
  if (nrow(pres_sf) == 0) {
    return(data.frame(
      modelo = temporada_modelo,
      mes = mes_presencia,
      n_total = 0,
      n_predicho = NA,
      porcentaje_predicho = NA
    ))
  }
  
  valores <- raster::extract(r, pres_sf)
  n_total <- length(valores)
  n_predicho <- sum(valores == 1, na.rm = TRUE)
  porcentaje <- round(n_predicho / n_total * 100, 2)
  
  return(data.frame(
    modelo = temporada_modelo,
    mes = mes_presencia,
    n_total = n_total,
    n_predicho = n_predicho,
    porcentaje_predicho = porcentaje
  ))
}

# -------------------------------
# Settings
# -------------------------------

modelos <- c("invierno", "reproduccion")
meses_migratorios <- c("febrero", "marzo", "abril", "julio", "agosto", "septiembre", "octubre")

dir_rasters <- "E:/maps/binarios_meses/"
dir_presencias <- "E:/data/"
dir_shapes <- "E:/m/"  

usar_shape <- TRUE # optional

# Run all comparisons
comparaciones <- expand.grid(modelo = modelos, mes = meses_migratorios)

resultados_meses <- purrr::pmap_dfr(
  list(comparaciones$modelo, comparaciones$mes),
  evaluar_interprediccion,
  dir_rasters = dir_rasters,
  dir_presencias = dir_presencias,
  dir_shapes = dir_shapes,
  usar_shape = usar_shape
)

print(resultados_meses)
