#######################################################################Time-specific niche modelling code: Rufus hummingbird niche tracking (Saldaña-Reyes et al., 2026)
#This code is based on Osorio-Olvera's tenm package tutorial (available on:https://github.com/luismurao/tenm)
######################################################################

library(tenm)
librayr(terra)

datos<- read.csv("Selasphorus_rufus_data.csv")

#extraction of environmental values for each observation
abt <- tenm::sp_temporal_data(occs = datos,
                              longitude = "lon.1",
                              latitude = "lat.1",
                              sp_date_var = "date_obs",
                              occ_date_format="ymd",
                              layers_date_format= "ym",
                              layers_by_date_dir = "enviornmental layers directory",
                              layers_ext="*.tif$")

#random selection of observations used to train and test the models
abex <- tenm::ex_by_date(this_species = abt,
                         train_prop = 0.7)

#save database with environmental values extracted from each register and training or testing identity of each observation
write.csv(abex$temporal_df, file = "temporal_df.csv")

#background generation
abbg <- tenm::bg_by_date(this_species = abex,
                         buffer_ngbs = 10,n_bg = 10000)

#### variables selection 
varcorrs <- tenm::correlation_finder(environmental_data =
                                       abex$env_data[,-ncol(abex$env_data)],
                                     method = "spearman",
                                     threshold = 0.8,
                                     verbose = FALSE)

#(Saldaña-Reyes et al., 2026 used monthly tmin, tmax and prec, so naturally there is variable correlation between tmax and tmin but authors wanted to keep these three variables following Peña-Peniche et al.,2018)

nombres_list_cor <- names(varcorrs$list_cor)
varcorrs$descriptors <- nombres_list_cor
print(varcorrs)
vars2fit <- varcorrs$descriptors
print(vars2fit)

##model calibration and evaluation
mod_sel <- tenm::tenm_selection(this_species = abbg,
                                omr_criteria = 0.1,
                                ellipsoid_level=0.975,
                                vars2fit = vars2fit,
                                nvars_to_fit=c(2,3,4),
                                proc = T,
                                RandomPercent = 50,
                                NoOfIteration= 1000,
                                parallel= T,
                                n_cores= n)
##see and save selected models
head(mod_sel$mods_table, 11)
setwd("directory")
write.csv(mod_sel$mods_table, "selected_models.csv", row.names=F)

##project selected model
temporal_layers_dir <-"directory of enviornmental layers"

proj_env_layers <- list.dirs(temporal_layers_dir,
                             recursive = F)[9] ## number in square brackets corresponds to the month-year folder used to project the model to e-space, in this example: September. 
suit <- predict(mod_sel,
                model_variables = c("wc2.1_2.5m_tmax",
                                    "wc2.1_2.5m_tmin", 
                                    "wc2.1_2.5m_prec"),
                layers_path =proj_env_layers,
                layers_ext = ".tif$")

##projection of the selected model to g-space
par(mfrow=c(1,2), mar=c(4,4,2,2))
terra::plot(suit, main = "Projection for XX month")
