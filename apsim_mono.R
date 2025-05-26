# LBRAI2219 - Modélisation des systèmes biologiques (2025)
#Code inspiré par le modèle APSIM fait par Alice Falzon et Maxime Cornez
###############################################################################
# Modèle
###############################################################################

# Ce script simule la croissance d'une culture en fonction de l'eau disponible et de la radiation solaire pour une seule culture

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
    LAI = numeric(n_days),
    I = numeric(n_days),
    Biomasse_1 = numeric(n_days),
    Biomasse_2 = numeric(n_days),
    O_D = numeric(n_days),
    Biomasse_reelle = numeric(n_days),
    Biomasse_cum = numeric(n_days),
    ES1 = numeric(n_days),
    ES2 = numeric(n_days),
    ES3 = numeric(n_days),
    rdepth = numeric(n_days)
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
    results$rdepth[i] <- rdepth
    
    
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
  # Stockage des réserves d'eau dans le data frame des résultats
  results$ES1 <- ES1
  results$ES2 <- ES2
  results$ES3 <- ES3
  
  # Calcul des indices et biomasses après simulation
  results$I <- 1 - exp(-culture_params$k * results$LAI)
  results$Biomasse_1 <- (results$Pot_Supply * culture_params$TEc) / (facteur_externe$VPDcalc / 10)
  results$Biomasse_2 <- facteur_externe$Radiation * culture_params$RUE * results$I
  results$O_D <- ifelse(results$Pot_Demand == 0, 0, results$Pot_Supply / results$Pot_Demand)
  results$Biomasse_reelle <- ifelse(results$O_D > 1, results$Biomasse_2, results$Biomasse_1)
  results$Biomasse_cum <- culture_params$Biomasse_initiale + cumsum(results$Biomasse_reelle)
  
  return(results)
}

# Appel de la fonction simulate
# resultats <- simulate(facteur_externe, soil_params, culture_params, expansion_foliaire)
# print(resultats)