#########################
# Importation des packages
#########################
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)


#########################
# Importation des données
#########################

culture1 <- list(RUE = 1.5,
                 TEc = 9,
                 VPR = 20, # mm/jour
                 CroissancePotLAI = 0.1,
                 k = 0.45,
                 LAI_initial = 1.5,
                 Biomasse_initial = 45) # g/m²

#optionnel
culture2 <- list(RUE = 1.5,  
                 TEc = 9,
                 VPR = 20,
                 CroissancePotLAI = 0.1,
                 k = 0.45,
                 LAI_initial = 1.5,
                 Biomasse_initial = 45)

sol_param <- data_frame(
  horizon = c(1, 2, 3),
  épaisseur = c(300, 300, 300),
  li = c(50, 50, 50),
  ls = c(100, 100, 100),
  es = c(100, 100, 100),
  kl = c(0.07, 0.07, 0.07)
)
