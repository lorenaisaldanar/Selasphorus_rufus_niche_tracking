################################################################
#Data cleaning process: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
#This code is an adaptation of Prieto-Torres, 2024
################################################################

library(rgbif)
library(CoordinateCleaner)
library(dplyr)
library(terra)
library(sf)
library(sp)
library(stringr)
library(tenm)

#################
#1. DOWNLOAD DATA
#################

##Accoridng to Wallen (2021) tutorial: https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html 
taxonKey <- name_backbone("Selasphorus rufus")$usageKey
occ_download(pred("hasGeospatialIssue", FALSE),
             pred("hasCoordinate", TRUE),
             pred("occurrenceStatus","PRESENT"), 
             pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
             pred("taxonKey", 2476855), 
             user= "user", 
             pwd= "password", 
             email= "e-mail" , 
             format = "SIMPLE_CSV") # 5219534 is the taxonKey 

###Accoridng to: https://www.gbif.org/es/tool/81747/rgbif 
occ_download_wait('')
setwd("your_directory")
#After it finishes, use
d <- occ_download_get('occ_download_wait key') %>%
  occ_download_import()
my_download_metadata <- occ_download_meta("occ_download_wait_key")
gbif_citation(my_download_metadata)

################
#2.FIRST FILTER
################
d <- readr::read_csv("occ_download_wait_key.csv")

species1 <- d %>%  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000)# remove zoo and herbaria within 2km 

#Eliminate the data without information to lon, lay, year, month and locality
data1 <- subset(species1, !is.na(decimalLatitude) & !is.na(decimalLongitude))##eliminate the data without geographical coordinates Latitud and longitud
data2 <- subset(data1, !is.na(year))##eliminate the data without information to year of collected
data3 <- subset(data2, !is.na(stateProvince))##eliminate the data without information to state/province of collected
data4 <-subset(data3, !is.na(month))##eliminate the data without information to year of collected


gbifID<-data4$gbifID ###corresponding to "GBIF code" to each sample.
institutionCode<-data4$collectionCode###corresponding to "Institution" information to each sample.
catalogNumber<-data4$catalogNumber###corresponding to "catalogue number" to each sample..
kingdom<-data4$kingdom###corresponding to "taxonomic information" to species.
phylum<-data4$phylum###corresponding to "taxonomic information" to species.
class<-data4$class###corresponding to "taxonomic information" to species.
order<-data4$order###corresponding to "taxonomic information" to species.
family<-data4$family###corresponding to "taxonomic information" to species.
genus<-data4$genus###corresponding to "taxonomic information" to species.
species<-data4$species###corresponding to "taxonomic information" to species.
lat<-data4$decimalLatitude###corresponding to "geographical coordinates" to each sample.
lon<-data4$decimalLongitude###corresponding to "geographical coordinates" to each sample.
year<-data4$year###corresponding to "year" of collect for each sample.
month<-data4$month##corresponding to "month" of collect for each sample
country<-data4$countryCode###corresponding to "country" of colect to each sample.
region<-data4$stateProvince###corresponding to "province/state" of colect to each sample

##To combine all variable in a sample file
data_cleaned1 <-data.frame(gbifID,institutionCode,catalogNumber,kingdom,phylum,class,order,family,genus,species,lon,lat,year,month,country,region)

######################
##3.SPATIAL FILTERING
######################

data_cleaned2 <- data_cleaned1 %>%
  group_by(species, year, month) %>%
  group_modify(~ tenm::clean_dup(
    .x,
    longitude = "lon",
    latitude = "lat",
    threshold = 0.0083,
    by_mask = FALSE
  )) %>%
  ungroup()

###################
##4. ENVIRONMENTAL FILTERING PROCEDURE FOR 2001-2022 DATA
##################

data_1950_2000 <- data_cleaned2 %>%
  filter(year >= 1950 & year <= 2000)

data_2001_actual <- data_cleaned2 %>%
  filter(year >= 2001)

dir_layers <- "layers_directory"

##list all enivronmental layers, corresponding to each month

env_files <- list.files(
  dir_layers,
  pattern = "\\.tif$",
  full.names = TRUE
)

##Assuming a structure like: tmax_01.tif, tmin_01.tif, prec_01.tif...

layers_tmin <- env_files[str_detect(basename(env_files), "tmin")]
layers_tmax <- env_files[str_detect(basename(env_files), "tmax")]
layers_prec <- env_files[str_detect(basename(env_files), "prec")]

# Order by month
layers_tmin <- layers_tmin[order(layers_tmin)]
layers_tmax <- layers_tmax[order(layers_tmax)]
layers_prec <- layers_prec[order(layers_prec)]

# stack 
rast_tmin <- rast(layers_tmin)
rast_tmax <- rast(layers_tmax)
rast_prec <- rast(layers_prec)

c <- rast(capas_prec)

##extract enviornmental values according to the month of observation

extract_values <- function(df) {
  
  # to create empty columns
  df$tmin <- NA
  df$tmax <- NA
  df$prec <- NA
  
  # Extract by month
  for(m in 1:12) {
    
    # Filter month observations
    idx <- which(df$month == m)
    
    if(length(idx) > 0) {
      
      # Create spatial points
      spatial_points <- vect(
        df[idx, ],
        geom = c("lon", "lat"),
        crs = "EPSG:4326"
      )
      
      # Extract values
      df$tmin[idx] <- extract(rast_tmin[[m]], spatial_points)[,2]
      df$tmax[idx] <- extract(rast_tmax[[m]], spatial_points)[,2]
      df$prec[idx] <- extract(rast_prec[[m]], spatial_points)[,2]
    }
  }
  
  return(df)
}

data_1950_2000_env <- extract_values(data_1950_2000)

data_2001_actual_env <- extract_values(data_2001_actual)

##calculate envrionmental historic threshold
thresholds <- list(
  
  tmin_min = min(data_1950_2000_env$tmin, na.rm = TRUE),
  tmin_max = max(data_1950_2000_env$tmin, na.rm = TRUE),
  
  tmax_min = min(data_1950_2000_env$tmax, na.rm = TRUE),
  tmax_max = max(data_1950_2000_env$tmax, na.rm = TRUE),
  
  prec_min = min(data_1950_2000_env$prec, na.rm = TRUE),
  prec_max = max(data_1950_2000_env$prec, na.rm = TRUE)
)

print(thresholds)

##filter recent observations
data_2001_filtered <- data_2001_actual_env %>%
  
  filter(
    
    tmin >= thresholds$tmin_min &
      tmin <= thresholds$tmin_max,
    
    tmax >= thresholds$tmax_min &
      tmax <= thresholds$tmax_max,
    
    prec >= thresholds$prec_min &
      prec <= thresholds$prec_max
  )

data_cleaned3 <- bind_rows(
  data_1950_2000_env,
  data_2001_filtered
)

write.csv(data_cleaned3, file = "Selasphorus_rufus_data.csv")

##save occurences as points (.shp)
df <- data_cleaned3
crs <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

points_sp <- sp::SpatialPointsDataFrame(
  coords = df[, c("lon", "lat")],
  data = df,  
  proj4string = crs 
)

pts_sf <- sf::st_as_sf(points_sp)
sf::st_write(pts_sf, "Selasphorus_rufus_data.shp")
