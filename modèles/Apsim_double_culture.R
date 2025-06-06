# LBRAI2219 - Modélisation des systèmes biologiques (2025)
#Code inspiré par le modèle APSIM fait par Alice Falzon et Maxime Cornez 
###############################################################################
# Chargement des packages
###############################################################################
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jsonlite)
library(httr)

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
# Passer kl dans la plante car il varie entre les espèces de plantes 

culture2_params <- list(RUE               = 1.5,    # Radiation Use Efficiency [g/MJ]
                        TEc               = 8,      # Coefficient d'efficience de la transpiration
                        VPR               = 10,     # Vitesse production de racines [mm/jour]
                        CroissPotLAI      = 0.2,    # Croissance potentielle du LAI 
                        k                 = 0.45,   # Coefficient d'extinction de la lumière
                        LAI_initial       = 1.8,   
                        Biomasse_initiale = 30
)
# Paramètres du sol (3 horizons, épaisseur, etc.) ------------------------------

soil_params <- data.frame(
  Horizon     = c(1, 2, 3),
  Epaisseur   = c(300, 300, 300),    # mm
  li          = c(50,  50,  50),     # Limite inférieure d'eau
  ls          = c(100, 100, 100),    # Limite supérieure d'eau
  es          = c(100, 100, 100),    # Eau disponible dans sol (= sw)
  es_h        = c(40, 30, 30),       # Eau disponible dans sol par horizon ; es_h1 + es_h2 + es_h3 = es
  kl          = c(0.07, 0.07, 0.07)  # Taux d'extraction [mm/jour]
)



# Expansion foliaire (table O/D - CEF)------------------------------------------
# Culture 1
expansion_foliaire <- data.frame(
  OD  = c(0.5, 1.5, 4), # Offre sur demande
  CEF = c(0,   1,   1)  # Coefficient d'expansion foliaire
)

# Culture 2
expansion_foliaire2 <- data.frame(
  OD  = c(0.4, 1.0, 3.5),
  CEF = c(0,   0.8, 1)
)

# Paramètres météo (ex. VPD Frac)-----------------------------------------------
# meteo_params <- list(
#   VPD_Frac = 1.0  
# )

# Vérification------------------------------------------------------------------

# Affiche les listes des paramètres de la culture
culture_params
culture2_params

# Affiche le data frame du sol
soil_params

# Affiche la table d’expansion foliaire
expansion_foliaire
expansion_foliaire2

###############################################################################
# Facteurs externes 
###############################################################################
# Données météo réelles 
lat <- 50.666265; lon <- 4.622322                              # Localisation de la zone d'étude : 50.666265 ; 4.622322 = batiment De serres
start_date <- "2024-06-01"; end_date <- "2024-06-30"           # Période de simulation
res <- GET("https://archive-api.open-meteo.com/v1/archive",    # Obtention des données de la station météo belge
           query = list(latitude = lat, longitude = lon,
                        start_date = start_date, end_date = end_date,
                        daily = "temperature_2m_max,temperature_2m_min,shortwave_radiation_sum",
                        timezone = "Europe/Brussels"))
meteo_data <- fromJSON(content(res, "text"))$daily             # Importation du fichier de données JSON

meteo_df <- as.data.frame(meteo_data) %>%
  rename(Date = time, Tmax = temperature_2m_max,
         Tmin = temperature_2m_min, Radiation = shortwave_radiation_sum) %>%
  mutate(Jours = 1:n(),
         svpTmax = 6.1078 * exp(17.269 * Tmax / (237.3 + Tmax)) * 0.10,
         svpTmin = 6.1078 * exp(17.269 * Tmin / (237.3 + Tmin)) * 0.10,
         VPDcalc = 0.75 * (svpTmax - svpTmin) * 10)

# Génération de données artificielles

facteur_externe <- data.frame(
  # Durée de la simulation (jours)
  Jours = 30:60, 
  
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


# Affiche la liste pour vérification
print(facteur_externe)

# Fonction pour calculer saturated vapour pressure (SVP)
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
# Choix de la source de données : "reel" ou "artificiel" (switch)
###############################################################################
data_source <- "artificiel"    # mettre "artificiel" si on veux le jeu généré

if (data_source == "reel") {
  facteur_externe <- meteo_df %>%
    select(Jours, Radiation, Tmax, Tmin, VPDcalc)
} else if (data_source == "artificiel") {
  facteur_externe <- facteur_externe %>%
    select(Jours, Radiation, Tmax, Tmin, VPDcalc)
} else {
  stop("data_source doit être 'reel' ou 'artificiel'")
}

###############################################################################
# Modèle
###############################################################################
# Densite des deux cultures
Densite1 <- 0.5 # Densité de la culture 1 
Densite2 <- 1-Densite1 # Densité de la culture 2 

simulate_two_cultures <- function(facteur_externe, soil_params, culture_params, culture2_params, expansion_foliaire, expansion_foliaire2, Densite1, Densite2) {
  n_days <- nrow(facteur_externe)
  
  # Initialisation du sol (réserves d'eau pour chaque horizon)
  ES1 <- numeric(n_days)
  ES2 <- numeric(n_days)
  ES3 <- numeric(n_days)
  ES1[1] <- soil_params$es_h[1]
  ES2[1] <- soil_params$es_h[2]
  ES3[1] <- soil_params$es_h[3]
  
  # Initialisation des profondeurs racinaires
  rdepth <- numeric(n_days)
  rdepth2 <- numeric(n_days)
  
  # Initialisation des LAI dynamiques pour chaque culture
  LAI1 <- numeric(n_days)
  LAI2 <- numeric(n_days)
  LAI1[1] <- culture_params$LAI_initial
  LAI2[1] <- culture2_params$LAI_initial
  
  # Initialisation de la biomasse pour chaque culture
  Biomasse1 <- numeric(n_days)
  Biomasse2 <- numeric(n_days)
  Biomasse1[1] <- culture_params$Biomasse_initiale
  Biomasse2[1] <- culture2_params$Biomasse_initiale
  
  # Initialisation des vecteurs de biomasse cumulée
  Biomasse_cum1 <- numeric(n_days)
  Biomasse_cum2 <- numeric(n_days)
  Biomasse_cum1[1] <- culture_params$Biomasse_initiale
  Biomasse_cum2[1] <- culture2_params$Biomasse_initiale
  
  # Préparation d'un data frame pour stocker les résultats quotidiens
  results <- data.frame(
    Jour = facteur_externe$Jours,
    Tot_ES = numeric(n_days),         # Eau totale dans le sol
    Pot_Supply = numeric(n_days),     # Offre potentielle globale du sol
    Pot_Demand1 = numeric(n_days),    # Demande potentielle culture 1
    Pot_Demand2 = numeric(n_days),    # Demande potentielle culture 2
    sdratio1 = numeric(n_days),       # Offre/Demande culture 1
    sdratio2 = numeric(n_days),       # Offre/Demande culture 2
    Transpiration1 = numeric(n_days), # Transpiration culture 1
    Transpiration2 = numeric(n_days), # Transpiration culture 2
    LAI1 = numeric(n_days),           # LAI culture 1
    LAI2 = numeric(n_days),           # LAI culture 2
    rdepth = numeric(n_days),         # Profondeur racinaire culture 1
    rdepth2 = numeric(n_days),        # Profondeur racinaire culture 2
    Biomasse1 = numeric(n_days),      # Biomasse journalière culture 1
    Biomasse2 = numeric(n_days)       # Biomasse journalière culture 2
  )
  
  results$LAI1[1] <- LAI1[1]            # Initialisation du LAI1
  results$LAI2[1] <- LAI2[1]            # Initialisation du LAI2
  results$Biomasse1[1] <- Biomasse1[1]  # Initialisation de la biomasse journalière culture 1
  results$Biomasse2[1] <- Biomasse2[1]  # Initialisation de la biomasse journalière culture 2
  
  # Fonction pour calculer l'effet d'expansion foliaire (interpolation sur la table OD/CEF)
  leaf_exp_effect <- function(sdratio, expansion_foliaire) {
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
  
  # Boucle de simulation quotidienne
  for (i in 1:n_days) {
    ## 1) Profondeurs racinaires
    DAS      <- facteur_externe$Jours[i]
    max_root <- sum(soil_params$Epaisseur)
    rdepth   <- min(DAS * culture_params$VPR,  max_root)
    rdepth2  <- min(DAS * culture2_params$VPR, max_root)
    results$rdepth[i]  <- rdepth
    results$rdepth2[i] <- rdepth2
    
    ## 2) Offres potentielles par horizon [1 à 3], par culture [c1 ; c2]
    of1_c1 <- ifelse(rdepth  >= soil_params$Epaisseur[1], 1, rdepth  / soil_params$Epaisseur[1]) * ES1[i] * soil_params$kl[1]
    of2_c1 <- ifelse(rdepth  <= soil_params$Epaisseur[1], 0,
                     ifelse(rdepth > soil_params$Epaisseur[1]+soil_params$Epaisseur[2],
                            1, (rdepth - soil_params$Epaisseur[1]) / soil_params$Epaisseur[2])) * ES2[i] * soil_params$kl[2]
    of3_c1 <- ifelse(rdepth  <= sum(soil_params$Epaisseur[1:2]), 0,
                     (rdepth - sum(soil_params$Epaisseur[1:2])) / soil_params$Epaisseur[3]) * ES3[i] * soil_params$kl[3]
    
    of1_c2 <- ifelse(rdepth2 >= soil_params$Epaisseur[1], 1, rdepth2 / soil_params$Epaisseur[1]) * ES1[i] * soil_params$kl[1]
    of2_c2 <- ifelse(rdepth2 <= soil_params$Epaisseur[1], 0,
                     ifelse(rdepth2 > soil_params$Epaisseur[1]+soil_params$Epaisseur[2],
                            1, (rdepth2 - soil_params$Epaisseur[1]) / soil_params$Epaisseur[2])) * ES2[i] * soil_params$kl[2]
    of3_c2 <- ifelse(rdepth2 <= sum(soil_params$Epaisseur[1:2]), 0,
                     (rdepth2 - sum(soil_params$Epaisseur[1:2])) / soil_params$Epaisseur[3]) * ES3[i] * soil_params$kl[3]
    
    Pot_Supply1 <- Densite1 * (of1_c1 + of2_c1 + of3_c1)  # offre potentielle culture 1
    Pot_Supply2 <- Densite2 * (of1_c2 + of2_c2 + of3_c2)  # offre potentielle culture 2
    Pot_Supply  <- Pot_Supply1 + Pot_Supply2  # totale
    
    ## 3) Effet lumineux
    li1 <- 1 - exp(-culture_params$k  * LAI1[i])
    li2 <- 1 - exp(-culture2_params$k * LAI2[i])
    
    results$li1[i] <- li1;  results$li2[i] <- li2
    
    ## 4) Demandes potentielles
    rad    <- facteur_externe$Radiation[i]
    VPDcalc <- facteur_externe$VPDcalc[i]
    Pot_Demand1 <- rad * li1 * culture_params$RUE  * (VPDcalc/10) / culture_params$TEc
    Pot_Demand2 <- rad * li2 * culture2_params$RUE * (VPDcalc/10) / culture2_params$TEc
    
    results$Pot_Demand1[i] <- Pot_Demand1
    results$Pot_Demand2[i] <- Pot_Demand2
    
    ## 5) Transpiration journalière par culture
    TT1     <- min(Pot_Supply1, Pot_Demand1)
    TT2     <- min(Pot_Supply2, Pot_Demand2)
    transp1 <- if (Pot_Demand1 > 0) TT1 * (Pot_Demand1 * Densite1) / (Pot_Demand1 * Densite1 + Pot_Demand2 * Densite2) else 0
    transp2 <- if (Pot_Demand2 > 0) TT2 * (Pot_Demand2 * Densite2) / (Pot_Demand1 * Densite1 + Pot_Demand2 * Densite2) else 0
    
    results$Transpiration1[i] <- TT1
    results$Transpiration2[i] <- TT2
    
    
    # Stockage de l'eau totale du sol
    results$Tot_ES[i] <- ES1[i] + ES2[i] + ES3[i]
    results$Pot_Supply[i] <- Pot_Supply
    
    # Mise à jour des LAI pour chaque culture
    # On calcule le ratio offre/demande pour chaque culture
    sdratio1 <- if (Pot_Demand1 > 0) Pot_Supply1 / Pot_Demand1 else 0 # =O/D culture 1
    sdratio2 <- if (Pot_Demand2 > 0) Pot_Supply2 / Pot_Demand2 else 0 # =O/D culture 2
    leaf_effect1 <- leaf_exp_effect(sdratio1, expansion_foliaire)
    leaf_effect2 <- leaf_exp_effect(sdratio2, expansion_foliaire2)
    delta_LAI1 <- leaf_effect1 * culture_params$CroissPotLAI
    delta_LAI2 <- leaf_effect2 * culture2_params$CroissPotLAI
    
    # Mise à jour du LAI pour chaque culture
    if(i < n_days) {
      LAI1[i + 1] <- LAI1[i] + delta_LAI1
      LAI2[i + 1] <- LAI2[i] + delta_LAI2
    }
    
    # Mise à jour des résultats
    results$sdratio1[i] <- sdratio1
    results$sdratio2[i] <- sdratio2
    results$LAI1[i] <- LAI1[i]
    results$LAI2[i] <- LAI2[i]
    
    # Mise à jour de la biomasse pour chaque culture
    # Calcul de la biomasse journalière pour chaque culture
    Biomasse_calc1 <- if (sdratio1 > 1) rad * li1 * culture_params$RUE else transp1 * (culture_params$TEc / (VPDcalc/10))
    Biomasse_calc2 <- if (sdratio2 > 1) rad * li2 * culture2_params$RUE else transp2 * (culture2_params$TEc / (VPDcalc/10))
    
    # Mise à jour de la biomasse cumulée pour chaque culture
    if(i < n_days) {
      Biomasse_cum1[i + 1] <- Biomasse_cum1[i] + Biomasse_calc1
      Biomasse_cum2[i + 1] <- Biomasse_cum2[i] + Biomasse_calc2
    }
    results$Biomasse1[i] <- Biomasse_calc1
    results$Biomasse2[i] <- Biomasse_calc2
    results$Biomasse_cum1[i] <- Biomasse_cum1[i]
    results$Biomasse_cum2[i] <- Biomasse_cum2[i]
    
    ## 9) Mise à jour des ES horizon par horizon
    if (i < n_days) {
      if (Pot_Supply1 + Pot_Supply2 > 0) {
        ES1[i+1] <- ES1[i] - ((of1_c1 / Pot_Supply1) * transp1 + (of1_c2 / Pot_Supply2) * transp2)
        ES2[i+1] <- ES2[i] - ((of2_c1 / Pot_Supply1) * transp1 + (of2_c2 / Pot_Supply2) * transp2)
        ES3[i+1] <- ES3[i] - ((of3_c1 / Pot_Supply1) * transp1 + (of3_c2 / Pot_Supply2) * transp2)
        } 
      else {
        ES1[i+1] <- ES1[i]
        ES2[i+1] <- ES2[i]
        ES3[i+1] <- ES3[i]
      }
    }
  }  
  
  return(results)
}

resultats_two <- simulate_two_cultures(facteur_externe, soil_params, culture_params, culture2_params, expansion_foliaire, expansion_foliaire2, Densite1, Densite2)
print(resultats_two)

###############################################################################
# Graphiques 
###############################################################################

# Graphe O/D
# On transforme le data frame en format long pour faciliter la visualisation avec ggplot2
data_long <- resultats_two %>% 
  select(Jour, sdratio1, sdratio2) %>% 
  pivot_longer(cols = c(sdratio1, sdratio2), 
               names_to = "variable", 
               values_to = "value")

ggplot(data_long, aes(x = Jour, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Jours", y = "O/D", title = "Courbes de sdratio1 et sdratio2") +
  theme_minimal()


# Graphe de l'eau totale dans le sol (ES)
plot(resultats_two$Jour, resultats_two$Tot_ES,
     type = "o", pch = 17, col = "cyan3",
     xlab = "Jours", ylab = "Eau disponible (mm)",
     main = "Eau disponible")


# Graphe des transpirations cumulées
# Calcul des cumuls pour chaque culture
resultats_two <- resultats_two %>%
  mutate(Transpiration1_cum = cumsum(Transpiration1),
         Transpiration2_cum = cumsum(Transpiration2))

# Conversion au format long
data_long_transpi <- resultats_two %>%
  select(Jour, Transpiration1_cum, Transpiration2_cum) %>%
  pivot_longer(cols = c(Transpiration1_cum, Transpiration2_cum),
               names_to = "Culture",
               values_to = "Transpiration_cum")

# Tracé du graphique avec ggplot2
ggplot(data_long_transpi, aes(x = Jour, y = Transpiration_cum, color = Culture)) +
  geom_line() +
  geom_point() +
  labs(x = "Jours",
       y = "Transpiration cumulée (mm)",
       title = "Transpiration cumulée") +
  theme_minimal()


# Graphe du LAI
data_long_LAI <- resultats_two %>% 
  select(Jour, LAI1, LAI2) %>% 
  pivot_longer(cols = c(LAI1, LAI2), names_to = "Culture", values_to = "LAI")

ggplot(data_long_LAI, aes(x = Jour, y = LAI, color = Culture)) +
  geom_line() +
  geom_point() +
  labs(x = "Jours", y = "LAI", title = "Evolution du LAI") +
  theme_minimal()


# Graphe des biomasses cumulées
data_long_biom <- resultats_two %>% 
  select(Jour, Biomasse_cum1, Biomasse_cum2) %>% 
  pivot_longer(cols = c(Biomasse_cum1, Biomasse_cum2), 
               names_to = "Culture", 
               values_to = "Biomasse_cumulee")

ggplot(data_long_biom, aes(x = Jour, y = Biomasse_cumulee, color = Culture)) +
  geom_line() +
  geom_point() +
  labs(x = "Jours", y = "Biomasse cumulée (g/m2)", title = "Biomasse cumulée") +
  theme_minimal()


# Graphe d'évolution la radiation et du VPDcalc
ggplot() +
  # Radiation
  geom_line(data = data.frame(Jour = resultats_two$Jour, Radiation = facteur_externe$Radiation),
            aes(x = Jour, y = Radiation, color = "Radiation")) +
  geom_point(data = data.frame(Jour = resultats_two$Jour, Radiation = facteur_externe$Radiation),
             aes(x = Jour, y = Radiation, color = "Radiation")) +
  # VPDcalc
  geom_line(data = data.frame(Jour = resultats_two$Jour, VPDcalc = facteur_externe$VPDcalc),
            aes(x = Jour, y = VPDcalc, color = "VPDcalc")) +
  geom_point(data = data.frame(Jour = resultats_two$Jour, VPDcalc = facteur_externe$VPDcalc),
             aes(x = Jour, y = VPDcalc, color = "VPDcalc")) +
  # Légende 
  scale_color_manual(name = "Légende",
                     values = c("Radiation" = "blue", "VPDcalc" = "red")) +
  labs(x = "Jours", y = "Valeurs") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Pot_supply == Transpiration1 + Transpiration2 ?
resultats_two$Pot_Supply >= resultats_two$Transpiration1 + resultats_two$Transpiration2
#différence ?
resultats_two$Pot_Supply - (resultats_two$Transpiration1 + resultats_two$Transpiration2)
