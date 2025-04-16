# Calcul des offres potentielles pour chaque horizon, pour chaque culture
of1_c1 <- ifelse(rdepth >= soil_params$Epaisseur[1], 1, rdepth / soil_params$Epaisseur[1]) * ES1[i] * soil_params$kl[1]
of2_c1 <- ifelse(rdepth <= soil_params$Epaisseur[1], 0,
                 ifelse(rdepth > soil_params$Epaisseur[1] + soil_params$Epaisseur[2],
                        1, (rdepth - soil_params$Epaisseur[1]) / soil_params$Epaisseur[2])) * ES2[i] * soil_params$kl[2]
of3_c1 <- ifelse(rdepth <= (soil_params$Epaisseur[1] + soil_params$Epaisseur[2]), 0, 
                 (rdepth - soil_params$Epaisseur[1] - soil_params$Epaisseur[2]) / soil_params$Epaisseur[3]) * ES3[i] * soil_params$kl[3]

of1_c2 <- ifelse(rdepth2 >= soil_params$Epaisseur[1], 1, rdepth2 / soil_params$Epaisseur[1]) * ES1[i] * soil_params$kl[1]
of2_c2 <- ifelse(rdepth2 <= soil_params$Epaisseur[1], 0,
                 ifelse(rdepth2 > soil_params$Epaisseur[1] + soil_params$Epaisseur[2],
                        1, (rdepth2 - soil_params$Epaisseur[1]) / soil_params$Epaisseur[2])) * ES2[i] * soil_params$kl[2]
of3_c2 <- ifelse(rdepth2 <= (soil_params$Epaisseur[1] + soil_params$Epaisseur[2]), 0, 
                 (rdepth2 - soil_params$Epaisseur[1] - soil_params$Epaisseur[2]) / soil_params$Epaisseur[3]) * ES3[i] * soil_params$kl[3]

Pot_Supply1 <- of1_c1 + of2_c1 + of3_c1
Pot_Supply2 <- of1_c2 + of2_c2 + of3_c2


sdratio1 <- if (Pot_Demand1 > 0) Pot_Supply1 / Pot_Demand1 else 0
sdratio2 <- if (Pot_Demand2 > 0) Pot_Supply2 / Pot_Demand2 else 0

Total_Supply <- Densite1 * Pot_Supply1 + Densite2 * Pot_Supply2
Total_Demand <- Densite1 * Pot_Demand1  + Densite2 * Pot_Demand2
Total_Transpiration <- min(Total_Supply, Total_Demand)
transp1 <- if(Total_Demand > 0) Total_Transpiration * (Densite1 * Pot_Demand1 / Total_Demand) else 0
transp2 <- if(Total_Demand > 0) Total_Transpiration * (Densite2 * Pot_Demand2 / Total_Demand) else 0

Transpi1 <- min(Pot_Supply1, Total_Demand)
Transpi2 <- min(Total_Supply2, Total_Demand)


ratio1 <- if (Pot_Supply1 > 0) of1_c1 / Pot_Supply1 else 0
ratio2 <- if (Pot_Supply2 > 0) of1_c2 / Pot_Supply2 else 0
ES1[i + 1] <- ES1[i] - (ratio1 * Transpi1 + ratio2 * Transpi2)
if(i < n_days) {
  if(Total_Supply > 0) {
    ES1[i + 1] <- ES1[i] - ((((of1_c1) / Pot_Supply1) * Transpi1)+(((of1_c2) / Pot_Supply2) * Transpi2)) # transp ou Transpi ?
    ES2[i + 1] <- ES2[i] - ((of2_c1 * Densite1 + of2_c2 * Densite2) / Total_Supply) * Total_Transpiration
    ES3[i + 1] <- ES3[i] - ((of3_c1 * Densite1 + of3_c2 * Densite2) / Total_Supply) * Total_Transpiration
  } else {
    ES1[i + 1] <- ES1[i]
    ES2[i + 1] <- ES2[i]
    ES3[i + 1] <- ES3[i]
  }
}

#Idée de Draye : calculer une demande, voir si elle fit avec le supply : si oui la racine prend tout et sinon elle est limitée par le supply ça revient à calculer une transpiration et on a plus besoin des offres 