# ===================================================================
# Geographic interprediction calculation between seasonal niches: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
# ===================================================================

##This code works based on authors file names and inner file structure. To run it with other file names and structure, adjust code. 
library(raster)
library(sf)
library(dplyr)
library(purrr)
library(stringr)

# -------------------------------
# Equivalence dictionary
# -------------------------------
nombres_cortos <- list(
  invierno = "invi",
  reproduccion = "repro",
  migprim = "migprim",
  migfall = "migfall"
)

# -------------------------------
# Main function
# -------------------------------

evaluar_interprediccion <- function(temporada_modelo, temporada_presencia, 
                                    dir_rasters, dir_presencias, dir_shapes = NULL, 
                                    usar_shape = FALSE) {
  
  modelo_corto <- nombres_cortos[[temporada_modelo]]
  presencia_corto <- nombres_cortos[[temporada_presencia]]
  
  raster_path <- file.path(dir_rasters, paste0(modelo_corto, "_to_", presencia_corto, "_bin.tif"))
  presencia_path <- file.path(dir_presencias, paste0("data_", temporada_presencia, ".csv"))
  
  if (!file.exists(raster_path)) {
    warning(paste("Raster not found:", raster_path))
    return(NULL)
  }
  if (!file.exists(presencia_path)) {
    warning(paste("Records file not found:", presencia_path))
    return(NULL)
  }
  
  # Read raster
  r <- raster(raster_path)
  
  # Read records of presence and convert them into spatial points
  presencias <- read.csv(presencia_path)
  pres_sf <- st_as_sf(presencias, coords = c("lon.1", "lat.1"), crs = 4326)
  pres_sf <- st_transform(pres_sf, crs = crs(r))
  
  # Shapefile cut (opctional)
  if (usar_shape && !is.null(dir_shapes)) {
    shp_path <- file.path(dir_shapes, paste0("m_",temporada_presencia, ".shp"))
    
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
      presencia = temporada_presencia,
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
    presencia = temporada_presencia,
    n_total = n_total,
    n_predicho = n_predicho,
    porcentaje_predicho = porcentaje
  ))
}

# -------------------------------
# Settings
# -------------------------------

temporadas <- c("invierno", "reproduccion", "migprim", "migfall")

dir_rasters <- "E:/maps/binarios_recortados"     # Adjust 
dir_presencias <- "E:/data/"
dir_shapes <- "rutas/shapes"            

usar_shape <- FALSE # Optional

comparaciones <- expand.grid(modelo = temporadas, presencia = temporadas) %>%
  filter(modelo != presencia)

# -------------------------------
# Run all comparisons
# -------------------------------

resultados <- purrr::pmap_dfr(
  list(comparaciones$modelo, comparaciones$presencia),
  evaluar_interprediccion,
  dir_rasters = dir_rasters,
  dir_presencias = dir_presencias,
  dir_shapes = dir_shapes,
  usar_shape = usar_shape
)

print(resultados)
