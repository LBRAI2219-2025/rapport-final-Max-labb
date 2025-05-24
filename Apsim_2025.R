# LBRAI2219 - Modélisation des systèmes biologiques (2025)
#Code inspiré par le modèle APSIM, fait par Alice Falzon et Maxime Cornez 
###############################################################################
# Chargement des packages
###############################################################################
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

###############################################################################
# Chargement des paramètres
###############################################################################

# Paramètres de la culture ----------------------------------------------------

culture_params <- list(
  RUE               = 1.25,  # Radiation Use Efficiency [g/MJ]
  TEc               = 9,     # Coefficient d'efficience de la transpiration
  VPR               = 20,    # Vitesse production de racines [mm/jour]
  CroissPotLAI      = 0.1,   # Croissance potentielle du LAI 
  k                 = 0.45,  # Coefficient d'extinction de la lumière
  LAI_initial       = 1.5,   
  Biomasse_initiale = 45     
)

# Paramètres du sol (3 horizons, épaisseur, etc.) ------------------------------

soil_params <- data.frame(
  Horizon     = c(1, 2, 3),
  Epaisseur   = c(300, 300, 300), # mm
  li          = c(50,  50,  50),  # Limite inférieure d'eau
  ls          = c(100, 100, 100), # Limite supérieure d'eau
  #es         = c(100, 100, 100), # Eau disponible dans sol (= sw)
  es_h        = c(40, 30, 30), # Eau disponible dans sol par horizon ; es_h1+es_h2+es_h3 = es
  kl          = c(0.07, 0.07, 0.07) # Taux d'extraction [mm/jour]
)



# Expansion foliaire (table O/D - CEF)------------------------------------------

expansion_foliaire <- data.frame(
  OD  = c(0.5, 1.5, 4),
  CEF = c(0,   1,   1)
)


# Paramètres météo (ex. VPD Frac)-----------------------------------------------

meteo_params <- list(
  VPD_Frac = 1.0  
)


# Vérification------------------------------------------------------------------

# Affiche la liste des paramètres de la culture
culture_params

# Affiche le data frame du sol
soil_params

# Affiche la table d’expansion foliaire
expansion_foliaire

###############################################################################
# Facteurs externes 
###############################################################################

facteur_externe <- data.frame(
  # Jours (30 à 60)
  Jours = c(
    30, 31, 32, 33, 34, 35, 
    36, 37, 38, 39, 40, 41, 
    42, 43, 44, 45, 46, 47, 
    48, 49, 50, 51, 52, 53, 
    54, 55, 56, 57, 58, 59, 60
  ),
  
  # Radiation (MJ/m2)
  Radiation = c(
    27, 27, 14, 24, 23, 21, 
    23, 25, 17, 14, 26, 26, 
    10, 26, 30, 27, 27, 29, 
    27, 26, 25, 23, 14, 25, 
    22, 20, 22, 24, 28, 25, 25
  ),
  
  # Tmax (°C)
  Tmax = c(
    32.3, 31.0, 26.6, 26.0, 26.6, 29.5, 
    30.8, 32.5, 32.3, 25.2, 27.8, 28.2, 
    27.3, 28.6, 28.6, 28.3, 27.6, 31.0, 
    35.0, 34.3, 31.2, 32.7, 29.9, 30.8, 
    31.2, 28.4, 27.7, 31.4, 33.0, 33.2, 32.7
  ),
  
  # Tmin (°C)
  Tmin = c(
    16.4, 15.8, 15.6, 10.0, 11.7, 13.0, 
    16.5, 13.8, 16.7, 16.7, 15.8, 12.8, 
    17.3, 11.3, 13.7, 13.4, 13.7, 12.2, 
    14.7, 20.4, 15.5, 17.9, 18.4, 16.0, 
    16.1, 18.0, 14.9, 16.0, 16.2, 17.2, 18.3
  ),
  
  # VPDobs (hPa)
  VPDobs = c(
    19, 17, 16, 14, 12, 15, 
    18, 19, 23, 21, 15, 14, 
    20, 16, 17, 17, 15, 15, 
    16, 17, 18, 21, 22, 17, 
    18, 20, 18, 19, 19, 20, 21
  )
)


# Affiche la liste pour vérifier
print(facteur_externe)

# Fonction pour calculer saturated vapour pressure
svp <- function(T) { # satured vapour pressure [kPa]
  6.1078 * exp(17.269 * T / (237.3 + T)) * 0.10
}

# Application à Tmax
facteur_externe$svpTmax <- svp(facteur_externe$Tmax) # pression saturante à Tmax

# Application à Tmin
facteur_externe$svpTmin <- svp(facteur_externe$Tmin) # pression saturante à Tmin

# Choisir la valeur de VPDfrac (par défaut on prend 0.75)
VPDfrac <- 0.75

# Fonction pour calculer VPDcalc
calc_VPDcalc <- function(svpTmax, svpTmin, VPDfrac) {
  VPDfrac * (svpTmax - svpTmin)*10
}

# Ajout de la colonne VPDcalc
facteur_externe$VPDcalc <- calc_VPDcalc(
  svpTmax  = facteur_externe$svpTmax,
  svpTmin  = facteur_externe$svpTmin,
  VPDfrac  = VPDfrac
)


###############################################################################
# Modèle
###############################################################################

simulate <- function(facteur_externe, soil_params, culture_params, expansion_foliaire) {
  n_days <- nrow(facteur_externe)
  
  # Vecteurs pour l'eau disponible par horizon
  ES1 <- numeric(n_days)
  ES2 <- numeric(n_days)
  ES3 <- numeric(n_days)
  
  # Initialisation des réserves d'eau pour chaque horizon
  ES1[1] <- soil_params$es_h[1]
  ES2[1] <- soil_params$es_h[2]
  ES3[1] <- soil_params$es_h[3]
  
  # Initialisation de la LAI dynamique
  LAI_dyn <- numeric(n_days)
  LAI_dyn[1] <- culture_params$LAI_initial
  
  # Préparation d'un data frame pour stocker les résultats totaux
  results <- data.frame(
    Jour = facteur_externe$Jours,
    Tot_ES = numeric(n_days),
    Pot_Supply = numeric(n_days),
    Pot_Demand = numeric(n_days),
    Transpiration = numeric(n_days),
    LAI = numeric(n_days)
  )
  results$LAI[1] <- LAI_dyn[1]
  
  # Fonction pour calculer l'effet d'expansion foliaire
  leaf_exp_effect <- function(sdratio, expansion_foliaire) {
    # Si le ratio est en-deçà du premier seuil, on renvoie la première valeur,
    # si au-dessus du dernier, on renvoie la dernière,
    # sinon on effectue une interpolation linéaire.
    if(sdratio <= expansion_foliaire$OD[1]) {
      return(expansion_foliaire$CEF[1])
    } else if (sdratio >= tail(expansion_foliaire$OD,1)) {
      return(tail(expansion_foliaire$CEF,1))
    } else {
      for(j in 1:(nrow(expansion_foliaire)-1)){
        if(sdratio >= expansion_foliaire$OD[j] && sdratio < expansion_foliaire$OD[j+1]){
          return(expansion_foliaire$CEF[j] + 
                   (sdratio - expansion_foliaire$OD[j]) * 
                   (expansion_foliaire$CEF[j+1] - expansion_foliaire$CEF[j]) / 
                   (expansion_foliaire$OD[j+1] - expansion_foliaire$OD[j]))
        }
      }
    }
  }
  
  # Boucle de simulation pour chaque jour
  for (i in 1:n_days) {
    DAS <- facteur_externe$Jours[i]
    # La profondeur racinaire effective est limitée par la somme des épaisseurs
    profondeur_totale <- sum(soil_params$Epaisseur)
    rdepth <- min(DAS * culture_params$VPR, profondeur_totale)
    
    # Calcul des offres potentielles pour chaque horizon
    of1 <- ifelse(rdepth >= soil_params$Epaisseur[1], 1, rdepth / soil_params$Epaisseur[1]) * ES1[i] * soil_params$kl[1]
    of2 <- ifelse(rdepth <= soil_params$Epaisseur[1], 0,
                  ifelse(rdepth > soil_params$Epaisseur[1] + soil_params$Epaisseur[2],
                         1, (rdepth - soil_params$Epaisseur[1]) / soil_params$Epaisseur[2])) * ES2[i] * soil_params$kl[2]
    of3 <- ifelse(rdepth <= (soil_params$Epaisseur[1] + soil_params$Epaisseur[2]), 0, 
                  (rdepth - soil_params$Epaisseur[1] - soil_params$Epaisseur[2]) / soil_params$Epaisseur[3]) * ES3[i] * soil_params$kl[3]
    
    Pot_Supply <- of1 + of2 + of3
    
    # Calcul de l'effet lumineux : LI = 1 - exp(-k * LAI_dyn)
    li <- 1 - exp(-culture_params$k * LAI_dyn[i])
    
    # Calcul de la demande potentielle
    rad <- facteur_externe$Radiation[i]
    VPDcalc <- facteur_externe$VPDcalc[i]
    Pot_Demand <- rad * li * culture_params$RUE * (VPDcalc / 10) / culture_params$TEc
    
    # Transpiration journalière
    transpiration <- min(Pot_Supply, Pot_Demand)
    
    # Stockage des résultats du jour
    results$Tot_ES[i] <- ES1[i] + ES2[i] + ES3[i]
    results$Pot_Supply[i] <- Pot_Supply
    results$Pot_Demand[i] <- Pot_Demand
    results$Transpiration[i] <- transpiration
    results$LAI[i] <- LAI_dyn[i]
    
    # Calcul du ratio offre/demande (si demande > 0)
    sdratio <- if (Pot_Demand > 0) Pot_Supply / Pot_Demand else 0
    # Calcul de l'effet d'expansion foliaire
    leaf_effect <- leaf_exp_effect(sdratio, expansion_foliaire)
    # delta LAI = leaf_effect * CroissPotLAI
    delta_LAI <- leaf_effect * culture_params$CroissPotLAI
    
    # Mise à jour de LAI pour le jour suivant
    if(i < n_days) {
      LAI_dyn[i + 1] <- LAI_dyn[i] + delta_LAI
    }
    
    # Mise à jour des réserves d'eau pour le jour suivant
    if (i < n_days) {
      if(Pot_Supply > 0) {
        ES1[i + 1] <- ES1[i] - (of1 / Pot_Supply) * transpiration
        ES2[i + 1] <- ES2[i] - (of2 / Pot_Supply) * transpiration
        ES3[i + 1] <- ES3[i] - (of3 / Pot_Supply) * transpiration
      } else {
        ES1[i + 1] <- ES1[i]
        ES2[i + 1] <- ES2[i]
        ES3[i + 1] <- ES3[i]
      }
    }
  }
  
  return(results)
}

# Appel de la fonction simulate
resultats <- simulate(facteur_externe, soil_params, culture_params, expansion_foliaire)
print(resultats)

###############################################################################
# Calcul de la biomasse
###############################################################################

resultats$I <- 1-exp(-culture_params$k*resultats$LAI)

# Calcul de la biomasse 
resultats$Biomasse_1 <- (resultats$Pot_Supply * culture_params$TEc)/ (facteur_externe$VPDcalc/10) # eau limitante 

resultats$Biomasse_2 <- facteur_externe$Radiation * culture_params$RUE * resultats$I # lumière limitante

resultats$O_D <- ifelse(resultats$Pot_Demand == 0, 0,
                        resultats$Pot_Supply / resultats$Pot_Demand)

resultats$Biomasse_reelle <- ifelse(
  test = resultats$O_D > 1,
  yes  = resultats$Biomasse_2,  # BioRadn
  no   = resultats$Biomasse_1   # BioEau
)

resultats$Biomasse_cumulee <- culture_params$Biomasse_initiale + cumsum(resultats$Biomasse_reelle)

###############################################################################
# Graphiques 
###############################################################################

plot(resultats$Jour, resultats$O_D, 
     type = "o", pch = 15, col = "magenta",
     xlab = "Jours", ylab = "O/D",
     main = "O/D")
abline(h = 1, lty = 2) 

plot(resultats$Jour, resultats$Tot_ES,
     type = "o", pch = 17, col = "cyan3",
     xlab = "Jours", ylab = "Eau disponible (mm)",
     main = "Eau disponible")

transpi_cum <- cumsum(resultats$Transpiration)

plot(resultats$Jour, transpi_cum,
     type = "o", pch = 18, col = "blue",
     xlab = "Jours", ylab = "Transpiration cumulée (mm)",
     main = "Transpiration cumulée")

plot(resultats$Jour, resultats$LAI,
     type = "o", pch = 18, col = "blue4",
     xlab = "Jours", ylab = "LAI",
     main = "LAI")

plot(resultats$Jour, resultats$Biomasse_cumulee,
     type = "o", pch = 15, col = "magenta",
     xlab = "Jours", ylab = "Biomasse cumulée (g/m2)",
     main = "Biomasse cumulée")


ggplot() +
  # Radiation
  geom_line(data = data.frame(Jour = resultats$Jour, Radiation = facteur_externe$Radiation),
            aes(x = Jour, y = Radiation, color = "Radiation")) +
  geom_point(data = data.frame(Jour = resultats$Jour, Radiation = facteur_externe$Radiation),
             aes(x = Jour, y = Radiation, color = "Radiation")) +
  # VPDcalc
  geom_line(data = data.frame(Jour = resultats$Jour, VPDcalc = facteur_externe$VPDcalc),
            aes(x = Jour, y = VPDcalc, color = "VPDcalc")) +
  geom_point(data = data.frame(Jour = resultats$Jour, VPDcalc = facteur_externe$VPDcalc),
             aes(x = Jour, y = VPDcalc, color = "VPDcalc")) +
  # Légende 
  scale_color_manual(name = "Légende",
                     values = c("Radiation" = "blue", "VPDcalc" = "red")) +
  labs(x = "Jours", y = "Valeurs") +
  theme_minimal() +
  theme(legend.position = "bottom")
