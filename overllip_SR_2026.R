#################################################
#Ellipsoidal niche parameters, environmental and geographic projections, overlap, and unique areas per season: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
#################################################

library(overllip)
library(ntbox)
library(dplyr)
library(terra)
library(raster)

#read database with information of environmental values and training/testing identity of each observation (created with the process on "tenm_SR_2026.r)
temporal_df <- read.csv("temporal_df.csv")

train_df <- temporal_df |> dplyr::filter(trian_test == "Train") |> na.omit()
test_df<- temporal_df |> dplyr::filter(trian_test == "Test") |> na.omit()

#selection of records from each season
ids_winter <- grep("[0-9][0-9][0-9][0-9]-[1][1-2]-[0][1]|[0-9][0-9][0-9][0-9]-[0][1]-[0][1]",train_df$layer_dates,value = F)
winter <- train_df[ids_winter,]
d1 <- winter
write.csv(d1,file = "data_winter.csv")

ids_breeding <- grep("[0-9][0-9][0-9][0-9]-[0][5-6]-[0][1]",train_df$layer_dates, value = F)
breeding <- train_df[ids_breeding,]
d2 <- breeding
write.csv(d2,file = "data_breeding.csv")

ids_fall <- grep("[0-9][0-9][0-9][0-9]-[0][7-9]-[0][1]|[0-9][0-9][0-9][0-9]-[1][0]-[0][1]",train_df$layer_dates, value = F)
fall <- train_df[ids_fall,]
d3 <- fall
write.csv(d3,file = "data_migfall.csv")

ids_spring <- grep("[0-9][0-9][0-9][0-9]-[0][2-4]-[0][1]",train_df$layer_dates, value = F)
spring <- train_df[ids_spring,]
d4 <- spring
write.csv(d4,file = "data_migprim.csv")

#based on best model environmental variables, authors include only maximum and minimum temperature for each month

d1<- d1[complete.cases(d1[, c("wc2.1_30s_tmax", "wc2.1_30s_tmin")]) &
                 is.finite(d1$wc2.1_30s_tmax) &
                 is.finite(d1$wc2.1_30s_tmin), ]

d2<- d2[complete.cases(d2[, c("wc2.1_30s_tmax", "wc2.1_30s_tmin")]) &
          is.finite(d2$wc2.1_30s_tmax) &
          is.finite(d2$wc2.1_30s_tmin), ]

d3 <- d3[complete.cases(d3[, c("wc2.1_30s_tmax", "wc2.1_30s_tmin")]) &
          is.finite(d3$wc2.1_30s_tmax) &
          is.finite(d3$wc2.1_30s_tmin), ]

d4<- d4[complete.cases(d4[, c("wc2.1_30s_tmax", "wc2.1_30s_tmin")]) &
          is.finite(d4$wc2.1_30s_tmax) &
          is.finite(d4$wc2.1_30s_tmin), ]

#############################################
#Seasonal ellipsoids parameter estimation (centroid and covariance matrix)
#############################################

nicho1 <- ntbox::cov_center(d1,mve = T,level = 0.975,vars = c("wc2.1_30s_tmax","wc2.1_30s_tmin"))

nicho2 <- ntbox::cov_center(d2,mve = T,level = 0.975,vars = c("wc2.1_30s_tmax", "wc2.1_30s_tmin"))

nicho3 <- ntbox::cov_center(d3,mve = T,level = 0.975,vars = c("wc2.1_30s_tmax","wc2.1_30s_tmin"))

nicho4 <- ntbox::cov_center(d4,mve = T,level = 0.975,vars = c("wc2.1_30s_tmax", "wc2.1_30s_tmin"))

#########################################
#workflow to estimate niche overlap
#########################################

d1_edata <- overllip::ellipsoid_data(centroid = nicho1$centroid,
                                     sigma = nicho1$covariance,cf = 0.975)
d2_edata <- overllip::ellipsoid_data(centroid = nicho2$centroid,
                                     sigma = nicho2$covariance,cf = 0.975)
d3_edata <- overllip::ellipsoid_data(centroid = nicho3$centroid,
                                     sigma = nicho3$covariance,cf = 0.975)
d4_edata <- overllip::ellipsoid_data(centroid = nicho4$centroid,
                                     sigma = nicho4$covariance,cf = 0.975)

ellipsoid_stack <- overllip::stack(d1_edata,d2_edata,d3_edata,d4_edata,ellipsoid_names = c("Winter","Breeding","Fall migration","Spring migration"))

env_data <- overllip::hypercube(ellipsoids = ellipsoid_stack,
                                n = 10000)

rmat <- overllip::stack_overlap(ellipsoid_stack = ellipsoid_stack,
                                env_data = env_data,
                                parallel = F)

#explore "rmat" slots to obtain Jaccard index value, ellipsoid #volume and volume of intersection 
rmat

overllip::plot(rmat, xlab= "Tmax", ylab="Tmin")

###########################################
#projection of each ellipsoid (seasonal niche) to the geographic space of the season or a season of interest
############################################

env_rdata <- raster::stack(list.files("directory of the environmental layers of the season of interest",
                                      pattern = ".tif$",
                                      full.names = TRUE)[c(2,3)]) 

#use here the name of the object that corresponds to the ellipsoid of each season. E.g "nicho1" corresponds to "winter" in this code
mProj <- ntbox::ellipsoidfit(env_rdata,
                             centroid = nicho1$centroid,
                             covar = nicho1$covariance,
                             level = 0.99,size = 3)

raster::plot(mProj$suitRaster)
points(pg[,c("longitude","latitude")],pch=20,cex=0.5)
terra::writeRaster(mProj$suitRaster, filename="map_season1_to_season1gspace.tif", overwrite = T)
suit_season <- mProj$suitRaster

##binarize this map by following the next process:

#call species records for the season projected
data_spp <- read.csv("data_season1.csv")

data_spp$lat.1 <- as.numeric(as.character(data_spp$lat))
data_spp$lon.1 <- as.numeric(as.character(data_spp$lon))

#convert presence data to points of presence (SpatialPoints)
points_occ_fin <- sp::SpatialPointsDataFrame(data_spp[,3:4], data_spp)

#extract suitability value of each point 
suit_values <- na.omit(raster::extract(suit_season, points_occ_fin))

#calculate the 10% training presence threshold
threshold <- quantile(suit_values, prob = 0.09)

#transform suit_season continuous map to a binary map 
suit_season_bin <- suit_season >= threshold
terra::writeRaster(present_bin, filename="season1_bin")

x11()
plot(suit_season_bin)
points(data_spp[,c("lon.1","lat.1")],pch=20,cex=0.5)

###################################
#workflow to calculate unique areas suitable per season 
###################################

#used as example these lines are for the season of winter
env_rdataUA <- raster::stack(list.files("directory of the environmental layers of the season of interest (winter)",
                                      pattern = ".tif$",
                                      full.names = TRUE)[c(2,3)])

#it is important to check if variable names are correct, fix it if not, for example: 
names(env_rdataUA)
names(env_rdataUA)[names(env_rdataUA) == "mean.1"] <- "wc2.1_2.5m_tmax"
names(env_rdataUA)[names(env_rdataUA) == "mean.2"] <- "wc2.1_2.5m_tmin"

rras <- overllip::stack_overlap(ellipsoid_stack = ellipsoid_stack,
                                env_data = env_rdataUA,
                                suitability_differences =TRUE,
                                parallel = F)
x11()
raster::plot(rras@geographic_intersection)

# reclassify the map of "Global Intersection" to be a mask from where to calculate unique areas suitable per season 

map <- rras@geographic_intersection$Global_intersection
x11()
raster::plot(map)
map[map==1] <- 100
map[map==0] <- 1
map[map==100] <- 0

terra::writeRaster(map, filename="mask_winter.tif")

#repeat the previous process to create the mask for the missing seasons

#multiply this mask by the suitability binary map of the season of interest obtained previously (line 93)
map_result_1 <- map_mask* suit_season_bin

terra::writeRaster(map_result_1, filename="unique_areas_winter.tif")

