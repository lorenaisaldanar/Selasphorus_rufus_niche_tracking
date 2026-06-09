###################################################################
##Calculate the distances to the niche centroid of the records for each season in the environmental space, and determine statistical significance of the difference observed between seasons: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026).
###################################################################

library(ggplot2)
#read databases of the records from each season, and filter the trainning data from the annual database

invierno <- read.csv("data_invierno.csv")
repro <- read.csv("data_reproduccion_1.csv")
migprim <- read.csv("data_migprim.csv")
migfall <- read.csv("data_migfall.csv")

invierno_train <- dplyr::filter(invierno, trian_test == "Train") |> na.omit()

repro_train <- dplyr::filter(repro, trian_test == "Train") |> na.omit()

migprim_train <- dplyr::filter(migprim, trian_test == "Train") |> na.omit()

migfall_train <- dplyr::filter(migfall, trian_test == "Train") |> na.omit()

train_df <- read.csv("temporal_df.csv") |> dplyr::filter(trian_test == "Train") |> na.omit()

#build the annual ellipsoid
nicho <- ntbox::cov_center(train_df,mve = T,level = 0.975,vars = c("wc2.1_30s_tmax",
                                                                   "wc2.1_30s_tmin"))

#calculate Mahalanobis distances of each record to the niche centroid
#wintering
dcentroid_invi <- mahalanobis(x = invierno[,names(nicho$centroid)],
                           center = nicho$centroid,cov = nicho$covariance)

## Data were randomly sampled within the ranges for each season so that all seasons would have the same amount of data for statistical analysis, based on the season with the fewest records, which is wintering (n=4165)
muestra_invierno <- sample(dcentroid_invi, size = 4165, replace = FALSE)
#breeding
dcentroid_repro <- mahalanobis(x = repro[,names(nicho$centroid)],
                           center = nicho$centroid,cov = nicho$covariance)

muestra_repro <- sample(dcentroid_repro, size = 4165, replace = FALSE)


#spring migration
dcentroid_migprim <- mahalanobis(x = migprim[,names(nicho$centroid)],
                               center = nicho$centroid,cov = nicho$covariance)

muestra_migprim <- sample(dcentroid_migprim, size = 4165, replace = FALSE)


#fall migration
dcentroid_migfall <- mahalanobis(x = migfall[,names(nicho$centroid)],
                                 center = nicho$centroid,cov = nicho$covariance)

muestra_migfall <- sample(dcentroid_migfall, size = 4165, replace = FALSE)



#Join all distances in one data frame
df_distnacias <- data.frame(distancia = c(muestra_invierno, muestra_repro, muestra_migprim, muestra_migfall),
                            temporada = c(rep("winter",length(muestra_invierno)),
                                          rep("breeding",length(muestra_repro)),
                                          rep("spring migration",length(muestra_migprim)),
                                          rep("fall migration",length(muestra_migfall))))

#########################
##Testing data normality 

# Iterate through each column of the dataframe
df_distnacias_ancho <- data.frame(
  winter = muestra_invierno,
  breeding = muestra_repro,
  spring_migration = muestra_migprim,
  fall_migration = muestra_migfall
)

#change the name each time, depending on the season to evaluate. For example, here fall migration
nombre_columna <- "fall_migration"

# Delete NA
df_sin_na <- df_distnacias_ancho[!is.na(df_distnacias_ancho[[nombre_columna]]), ]

# Kolmogorov-Smirnov test 
resultado_ks <- ks.test(df_sin_na[[nombre_columna]], "pnorm",
                        mean = mean(df_sin_na[[nombre_columna]]),
                        sd = sd(df_sin_na[[nombre_columna]]))

cat(paste("Kolmogorov-Smirnov test for:", nombre_columna, "/n"))
print(resultado_ks)

###########################
# Kruskall-Walis test
###########################
library(ggplot2)

head(df_distnacias)

# Verify that column "temporada" is a factor
df_distnacias$temporada <- as.factor(df_distnacias$temporada)

# Kruskal-Wallis test
kruskal_test <- kruskal.test(distancia ~ temporada, data = df_distnacias)

print(kruskal_test)

# If the test results significant, perform a post-hoc analysis 
if (kruskal_test$p.value < 0.05) {
  # Use Dunn test for multiple comparisons (with Bonferroni correction)
  library(FSA)
  dunn_result <- dunnTest(distancia ~ temporada, data = df_distnacias, method = "bonferroni")
  print(dunn_result)
}

# Visualice distributions per season
orden_temporadas <- c("breeding", "fall migration", "winter", "spring migration")
df_distnacias$temporada <- factor(df_distnacias$temporada,
                                  levels = orden_temporadas)

dist_niche_centroid <- ggplot(df_distnacias, aes(x = temporada, y = distancia, fill = temporada)) +
  geom_boxplot() +
  scale_fill_manual(values = c("winter" = "#79b7f2", 
                               "breeding" = "#d69dbd",
                               "spring migration" = "#f8d800", 
                               "fall migration" = "#f6a721")) +  
  theme_minimal() +
  labs(x = "Season",
       y = "Mahalanobis distance") +
  guides(fill = "none")

plot(dist_niche_centroid)

#### to visualize them in the environmental space 

plot(x = train_df$wc2.1_30s_tmax, y = train_df$wc2.1_30s_tmin, type = "n",
     xlab = "T max", ylab = "T min")

tenm::plot_ellipsoid(x = train_df$wc2.1_30s_tmax, y = train_df$wc2.1_30s_tmin, semiaxes = FALSE, col="#000")

points(invierno_train$wc2.1_30s_tmax, invierno_train$wc2.1_30s_tmin, pch=19, col= "#79b7f2" )
points(migfall_train$wc2.1_30s_tmax, migfall_train$wc2.1_30s_tmin,pch=19, col= "#f6a721" )
points(migprim_train$wc2.1_30s_tmax, migprim_train$wc2.1_30s_tmin,pch=19, col= "#f8d800" )
points(repro_train$wc2.1_30s_tmax, repro_train$wc2.1_30s_tmin, pch=19, col= "#d69dbd" )
points(nicho[[1]][1],nicho[[1]][2],pch=19, col="#000")#plot the niche centroid

#plot the densities
df_distancias_plot <- ggplot(df_distnacias, aes(x = distancia, fill = temporada, color = temporada)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(x = "Distancia al centroide",
       y = "Densidad")

plot(df_distancias_plot)

