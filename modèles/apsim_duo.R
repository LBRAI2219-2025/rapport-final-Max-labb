# LBRAI2219 - Modélisation des systèmes biologiques (2025)
#Code inspiré par le modèle APSIM fait par Alice Falzon et Maxime Cornez 
# Ce script simule la croissance de deux cultures en interculture en tenant compte de l'eau disponible et de la radiation solaire

###############################################################################
# Modèle
###############################################################################

# Densite des deux cultures
Densite1 <- 0.5 # Densité de la culture 1 
Densite2 <- 1-Densite1 # Densité de la culture 2 

simulate_two <- function(facteur_externe, soil_params, culture_params, culture2_params, expansion_foliaire, expansion_foliaire2, Densite1, Densite2) {
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
    Tot_ES = numeric(n_days),              # Eau totale dans le sol
    Pot_Supply = numeric(n_days),          # Offre potentielle globale du sol
    Pot_Demand1 = numeric(n_days),         # Demande potentielle culture 1
    Pot_Demand2 = numeric(n_days),         # Demande potentielle culture 2
    sdratio1 = numeric(n_days),            # Offre/Demande culture 1
    sdratio2 = numeric(n_days),            # Offre/Demande culture 2
    Transpiration1_loc = numeric(n_days),  # Transpiration culture 1 si elle était seule
    Transpiration2_loc = numeric(n_days),  # Transpiration culture 2 si elle était seule
    Transpiration1 = numeric(n_days),      # Transpiration journalière culture 1 avec arbitrage
    Transpiration2 = numeric(n_days),      # Transpiration journalière culture 2 avec arbitrage
    LAI1 = numeric(n_days),                # LAI culture 1
    LAI2 = numeric(n_days),                # LAI culture 2
    rdepth = numeric(n_days),              # Profondeur racinaire culture 1
    rdepth2 = numeric(n_days),             # Profondeur racinaire culture 2
    Biomasse1 = numeric(n_days),           # Biomasse journalière culture 1
    Biomasse2 = numeric(n_days),           # Biomasse journalière culture 2
    ES1 = numeric(n_days),                 # Eau disponible horizon 1
    ES2 = numeric(n_days),                 # Eau disponible horizon 2
    ES3 = numeric(n_days)                  # Eau disponible horizon 3
  )
  
  results$LAI1[1] <- LAI1[1]               # Initialisation du LAI1
  results$LAI2[1] <- LAI2[1]               # Initialisation du LAI2
  results$Biomasse1[1] <- Biomasse1[1]     # Initialisation de la biomasse journalière culture 1
  results$Biomasse2[1] <- Biomasse2[1]     # Initialisation de la biomasse journalière culture 2
  
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
    
    Pot_Supply1 <-  (of1_c1 + of2_c1 + of3_c1)  # offre potentielle culture 1
    Pot_Supply2 <-  (of1_c2 + of2_c2 + of3_c2)  # offre potentielle culture 2
    Pot_Supply  <- Pot_Supply1 + Pot_Supply2    # totale
    
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
    TT1     <- min(Pot_Supply1, Pot_Demand1) # Transpiration locale de la culture 1
    TT2     <- min(Pot_Supply2, Pot_Demand2) # Transpiration locale de la culture 2
    transp1 <- if (Pot_Demand1 > 0) TT1 * (Pot_Demand1 * Densite1) / (Pot_Demand1 * Densite1 + Pot_Demand2 * Densite2) else 0    # Transpiration culture 1 avec arbitrage de l'eau
    transp2 <- if (Pot_Demand2 > 0) TT2 * (Pot_Demand2 * Densite2) / (Pot_Demand1 * Densite1 + Pot_Demand2 * Densite2) else 0    # Transpiration culture 2 avec arbitrage de l'eau
    
    results$Transpiration1_loc[i] <- TT1
    results$Transpiration2_loc[i] <- TT2
    results$Transpiration1[i] <- transp1
    results$Transpiration2[i] <- transp2
    
    
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
    results$Biomasse1[i] <- Biomasse_calc1        # biomasse produite au jour i par la culture 1
    results$Biomasse2[i] <- Biomasse_calc2        # biomasse produite au jour i par la culture 2
    results$Biomasse_cum1[i] <- Biomasse_cum1[i]  # biomasse cumulée au jour i par la culture 1
    results$Biomasse_cum2[i] <- Biomasse_cum2[i]  # biomasse cumulée au jour i par la culture 2
    
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
  
  # Stockage des réserves d'eau dans le data frame des résultats
  results$ES1 <- ES1
  results$ES2 <- ES2
  results$ES3 <- ES3
  
  return(results)
}

