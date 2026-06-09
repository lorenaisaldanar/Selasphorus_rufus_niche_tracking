######################################################################
#GMM algorithm application: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
######################################################################

library(mclust)
library(dplyr)

#read observation data of each season
invierno <- read.csv("data_invierno.csv")
reproduccion <- read.csv("data_reproduccion.csv")
migprim <- read.csv("data_migprim.csv")
migfall <- read.csv("data_migfall.csv")

data <- rbind(invierno,migfall, reproduccion, migprim)|> na.omit()

data_selected <- data |> 
  select(wc2.1_30s_tmax, wc2.1_30s_tmin)

#apply gmm algorithm
modelo_gmm <- Mclust(data_selected)

#see results
summary<- capture.output(summary(modelo_gmm, parameters = TRUE))

x11()
plot(modelo_gmm, what = "classification")

#to calculate proportion of each season on each cluster detected
data$layer_dates <- as.Date(data$layer_dates)
data$mes <- format(data$layer_dates, "%m")

# Asign season according to the month
data$temporada <- with(data, ifelse(mes %in% c("11", "12", "01"), "invierno",
                                            ifelse(mes %in% c("05", "06"), "reproduccion",
                                            ifelse(mes %in% c("02", "03", "04"), "migracion_primavera",
                                            ifelse(mes %in% c("07", "08", "09", "10"), "migracion_otono", NA)))))


data$cluster_gmm <- modelo_gmm$classification

# Cross table and proportion
tabla_cruce <- table(data$temporada, data$cluster_gmm)
prop_cluster <- prop.table(tabla_cruce, margin = 2)

# to visualize
library(ggplot2)
library(reshape2)

df_barras <- reshape2::melt(prop_cluster)
colnames(df_barras) <- c("Temporada", "Cluster", "Proporcion")

#asign colors to each season
colores_temporadas <- c(
  "invierno" = "#79b7f2",           # Azul (winter)
  "reproduccion" = "#d69dbd",       # Rosa (breeding)
  "migracion_primavera" = "#f8d800",# Verde (spring migration)
  "migracion_otono" = "#f6a721"     # Naranja (fall migration)
)

ggplot(df_barras, aes(x = factor(Cluster), y = Proporcion, fill = Temporada)) +
  geom_bar(stat = "identity") +
  labs(title = "Proporción de temporadas por cluster GMM",
       x = "Cluster GMM", y = "Proporción") +
  scale_fill_manual(values = colores_temporadas) +
  theme_minimal()


# Create contingency table
tabla_cruce <- table(data$temporada,data$cluster_gmm)

# Independence test
test_chi <- chisq.test(tabla_cruce)
print(test_chi)

residuos <- chisq.test(tabla_cruce)$stdres
print(residuos)

#visualize person standard residuals
library(vcd)
dir_out <- ""
file_out <- file.path(dir_out, "residuos_estandarizados_proporcion_temporadas.png")
jpeg(filename = file_out, width = 6, height = 6, units = "in", res = 300)  
assoc(tabla_cruce, shade = TRUE)
dev.off()
