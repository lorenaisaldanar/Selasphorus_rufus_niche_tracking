######################################################################
#DBSCAN algorithm application: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
######################################################################

library(dbscan)
library(ggplot2)
library(dplyr)

#read observation data of each season
invierno <- read.csv("data_invierno_1.csv")
reproduccion <- read.csv("data_reproduccion_1.csv")
migprim <- read.csv("data_migprim.csv")
migfall <- read.csv("data_migfall.csv")

data <- rbind(invierno,migfall, reproduccion, invierno) 

#select enviornmental information of each record, contained on each variable column
data <- data |> 
  select(wc2.1_30s_tmax, wc2.1_30s_tmin)

#verify no NA
sum(is.na(data))
data <- na.omit(data)

##CALCULATE EPSILON FOR DBSCAN
# Define the number of neighbors (k) that correspond to MinPts - 1
k <- 4  # For example, if MinPts = 5, then k = 4

# Calculate the distances to the k-th nearest neighbor
distancias <- kNNdist(data, k = k)

# Order distances 
distancias_ordenadas <- sort(distancias)

# Create a data frma to help with visualization
df_distancias <- data.frame(
  índice = 1:length(distancias_ordenadas),
  distancia = distancias_ordenadas
)

# Plot distances
ggplot(df_distancias, aes(x = índice, y = distancia)) +
  geom_line() +
  labs(
    x = "Ordered Points Index",
    y = "Distances to the k-th nearest neighbor"
  ) +
  theme_minimal()

#DBSCAN
dbscan_res <- dbscan(data, eps = 0.3, minPts = 5) # epsilon based on the elbow point in the previous plot

plot(data, col=dbscan_res$cluster+1, main="DBSCAN") 
