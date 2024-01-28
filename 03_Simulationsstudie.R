#### Simulationsstudie
rm(list = ls())

#Packages
if (!require('LNIRT')) install.packages('LNIRT')
if (!require('xlsx')) install.packages('xlsx')

library(xlsx)
library(LNIRT)

#Funktionen
load("functions/functions.RData")
load("functions/mml_hm.RData")

#Set Seed
set.seed(1234)

#### Vorbereitung ####
#Df mit Informationen für die Erstellung der Simulationsdaten
df_info_sim <- as.data.frame(matrix(NA, 4, 3))
names(df_info_sim) <- c("P", "N", "rho")
df_info_sim[1] <- c(rep(500, 2), rep(1000, 2))
df_info_sim[2] <- rep(20, 4)
df_info_sim[3] <- rep(c(0.00, 0.50), 2)

#Df für Parameter-Recovery
df_mean_para <- as.data.frame(matrix(NA, 4, 5 * 3))
names(df_mean_para) <-
  paste(rep(c("alpha", "beta", "omega", "ny", "rho"), each = 3),
        rep(c("True", "MML", "LNIRT"), 5))

#Df für Zeitmessung der Funktionen
df_time <- as.data.frame(matrix(NA, 4, 2))
names(df_time) <- c("MML", "LNIRT")

#Liste um Daten und Modelle zu speichern
list_mod <- list()


#### Ergebnisse für Parameter-Recovery ####
for (i in 1:nrow(df_info_sim)) {
  # Simulationsdaten erstellen
  data_sim <-
    simLNIRT(
      N = df_info_sim[i, 1],
      K = df_info_sim[i, 2],
      rho = df_info_sim[i, 3],
      WL = T
    )
  list_mod[[paste("Simulation", i, ":", "Data")]] <- data_sim
  
  # MML-Funktion anwenden
  start_mml <- Sys.time()
  # die Konstante D wird auf den Wert
  #1,702 gesetzt, um das zweiparametrische logistische Modell an das zweiparametrische
  #normal-ogive Modell anzunähern
  res_mml <- mml.hm(data_sim$Y, exp(data_sim$RT), D = 1.702)
  end_mml <- Sys.time()
  df_time[i, 1] <- difftime(end_mml, start_mml, units = "mins")
  list_mod[[paste("Simulation", i, ":", "MML")]] <- res_mml
  
  # LNIRT-Funktion anwenden
  start_lnirt <- Sys.time()
  res_lnirt <- LNIRT(data_sim$RT, data_sim$Y, WL = T, XG=5000)
  end_lnirt <- Sys.time()
  df_time[i, 2] <- difftime(end_lnirt, start_lnirt, units = "mins")
  list_mod[[paste("Simulation", i, ":", "LNIRT")]] <- res_lnirt
  
  ### True befüllen
  df_mean_para[i, grep("True", names(df_mean_para))] <-
    c(
      mean(data_sim$ab[, 1]),
      mean(data_sim$ab[, 2]),
      mean(data_sim$ab[, 3]),
      mean(data_sim$ab[, 4]),
      df_info_sim[i, 3]
    )
  ### MML befüllen
  df_mean_para[i, grep("MML", names(df_mean_para))] <-
    c(
      mean(res_mml$Item_Parameter[, 1]),
      # Die Beta-Werte werden nachträglich
      # mit -1 multipliziert, um die gleiche Parametrisierung zu erhalten, die auch im LNIRT-Package
      # verwendet wird
      -mean(res_mml$Item_Parameter[, 2]),
      mean(res_mml$Item_Parameter[, 3]),
      mean(res_mml$Item_Parameter[, 4]),
      res_mml$Kovarianzmatrix[1, 2] / (sqrt(res_mml$Kovarianzmatrix[2, 2]))
    )
  ### LNIRT befüllen
  df_mean_para[i, grep("LNIRT", names(df_mean_para))] <-
    c(
      mean(res_lnirt$Post.Means$Item.Discrimination),
      mean(res_lnirt$Post.Means$Item.Difficulty),
      mean(res_lnirt$Post.Means$Sigma2),
      mean(res_lnirt$Post.Means$Time.Intensity),
      res_lnirt$Post.Means$Cov.Person.Ability.Speed / (
        sqrt(res_lnirt$Post.Means$Var.Person.Speed) * sqrt(res_lnirt$Post.Means$Var.Person.Ability)
      )
    )
}

df_mean_para_rd <- round(df_mean_para, 2)
df_time


### Ergebnisse speichern
save(df_mean_para_rd, df_mean_para, df_time, list_mod, file="output/Simulation_Parameter_Recovery.RData")

write.xlsx(df_mean_para_rd, file="output/Simulation_Parameter_Recovery.xlsx")

######### End of Skript ###