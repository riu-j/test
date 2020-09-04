library("stargazer")

data = read.csv("../data/data_trim.csv", na.strings=".")
cols = c("Sample size", "Age", "female", "rs76904798", "G2019S", "rs356181", "rs3910105", "DaT scan", "MDS-UPDRS Part1", 
        "MDS-UPDRS Part2", "MDS-UPDRS Part3", "SCOPA-AUT", "MoCA", "ESS", "α-synuclein", "Aβ1-42", "p-tau", "T-tau")
PD = data[data$Diagnosis==1,]
GenPD = data[data$Diagnosis==5,]
RegPD = data[data$Diagnosis==7,]

result = data.frame(value=cols, PD = rep(0, length(cols)), PD2 = rep(0, length(cols)), 
                    GenPD = rep(0, length(cols)), GenPD2 = rep(0, length(cols)),
                    RegPD = rep(0, length(cols)), RegPD2 = rep(0, length(cols)), 
                    total = rep(0, length(cols)), total2 = rep(0, length(cols)))

result[1, 2] = length(unique(PD$PATNO))
result[1, 4] = length(unique(GenPD$PATNO))
result[1, 6] = length(unique(RegPD$PATNO))
result[1, 8] = length(unique(data$PATNO))

data2 = data[!duplicated(data$PATNO),]
PD2 = PD[!duplicated(PD$PATNO),]
GenPD2 = GenPD[!duplicated(GenPD$PATNO),]
RegPD2 = RegPD[!duplicated(RegPD$PATNO),]
result[2, 2] = mean(PD2$age)
result[2, 3] = sd(PD2$age)
result[2, 4] = mean(GenPD2$age)
result[2, 5] = sd(GebPD2$age)
result[2, 6] = mean(RegPD2$age)
result[2, 7] = sd(RegPD2$age)
result[2, 8] = mean(data2$age)
result[2, 9] = sd(data2$age)

result[3, 2] = nrow(PD2[PD2$Sex==1,])
result[3, 3] = result[3, 2] / result[1, 2]
result[3, 4] = nrow(GenPD2[GenPD2$Sex==1,])
result[3, 5] = result[3, 4] / result[1, 4]
result[3, 6] = nrow(RegPD2[RegPD2$Sex==1,])
result[3, 7] = result[3, 6] / result[1, 6]
result[3, 8] = nrow(data2[data2$Sex==1,])
result[3, 9] = result[3, 8] / result[1, 8]

result[4, 2] = nrow(PD2[PD2$rs76904798==1,])
result[4, 3] = result[4, 2] / result[1, 2]
result[4, 4] = nrow(GenPD2[GenPD2$rs76904798==1,])
result[4, 5] = result[4, 4] / result[1, 4]
result[4, 6] = nrow(RegPD2[RegPD2$rs76904798==1,])
result[4, 7] = result[4, 6] / result[1, 6]
result[4, 8] = nrow(data2[data2$rs76904798==1,])
result[4, 9] = result[4, 8] / result[1, 8]

result[5, 2] = nrow(PD2[PD2$G2019S==1,])
result[5, 3] = result[5, 2] / result[1, 2]
result[5, 4] = nrow(GenPD2[GenPD2$G2019S==1,])
result[5, 5] = result[5, 4] / result[1, 4]
result[5, 6] = nrow(RegPD2[RegPD2$G2019S==1,])
result[5, 7] = result[5, 6] / result[1, 6]
result[5, 8] = nrow(data2[data2$G2019S==1,])
result[5, 9] = result[5, 8] / result[1, 8]

result[6, 2] = nrow(PD2[PD2$rs356181==1,])
result[6, 3] = result[6, 2] / result[1, 2]
result[6, 4] = nrow(GenPD2[GenPD2$rs356181==1,])
result[6, 5] = result[6, 4] / result[1, 4]
result[6, 6] = nrow(RegPD2[RegPD2$rs356181==1,])
result[6, 7] = result[6, 6] / result[1, 6]
result[6, 8] = nrow(data2[data2$rs356181==1,])
result[6, 9] = result[6, 8] / result[1, 8]

result[7, 2] = nrow(PD2[PD2$rs3910105==1,])
result[7, 3] = result[7, 2] / result[1, 2]
result[7, 4] = nrow(GenPD2[GenPD2$rs3910105==1,])
result[7, 5] = result[7, 4] / result[1, 4]
result[7, 6] = nrow(RegPD2[RegPD2$rs3910105==1,])
result[7, 7] = result[7, 6] / result[1, 6]
result[7, 8] = nrow(data2[data2$rs3910105==1,])
result[7, 9] = result[7, 8] / result[1, 8]

tmpdata = data[!is.na(data$mean_striatum),]
tmpPD = PD[!is.na(PD$mean_striatum),]
tmpGenPD = GenPD[!is.na(GenPD$mean_striatum),]
tmpRegPD = RegPD[!is.na(RegPD$mean_striatum),]
result[8, 2] = length(unique(tmpPD$PATNO))
result[8, 3] = nrow(tmpPD)
result[8, 4] = length(unique(tmpGenPD$PATNO))
result[8, 5] = nrow(tmpGenPD)
result[8, 6] = length(unique(tmpRegPD$PATNO))
result[8, 7] = nrow(tmpRegPD)
result[8, 8] = length(unique(data$PATNO))
result[8, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$mds.updrs1_odds),]
tmpPD = PD[!is.na(PD$mds.updrs1_odds),]
tmpGenPD = GenPD[!is.na(GenPD$mds.updrs1_odds),]
tmpRegPD = RegPD[!is.na(RegPD$mds.updrs1_odds),]
result[9, 2] = length(unique(tmpPD$PATNO))
result[9, 3] = nrow(tmpPD)
result[9, 4] = length(unique(tmpGenPD$PATNO))
result[9, 5] = nrow(tmpGenPD)
result[9, 6] = length(unique(tmpRegPD$PATNO))
result[9, 7] = nrow(tmpRegPD)
result[9, 8] = length(unique(data$PATNO))
result[9, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$mds.updrs2_odds),]
tmpPD = PD[!is.na(PD$mds.updrs2_odds),]
tmpGenPD = GenPD[!is.na(GenPD$mds.updrs2_odds),]
tmpRegPD = RegPD[!is.na(RegPD$mds.updrs2_odds),]
result[10, 2] = length(unique(tmpPD$PATNO))
result[10, 3] = nrow(tmpPD)
result[10, 4] = length(unique(tmpGenPD$PATNO))
result[10, 5] = nrow(tmpGenPD)
result[10, 6] = length(unique(tmpRegPD$PATNO))
result[10, 7] = nrow(tmpRegPD)
result[10, 8] = length(unique(data$PATNO))
result[10, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$mds.updrs3_odds),]
tmpPD = PD[!is.na(PD$mds.updrs3_odds),]
tmpGenPD = GenPD[!is.na(GenPD$mds.updrs3_odds),]
tmpRegPD = RegPD[!is.na(RegPD$mds.updrs3_odds),]
result[11, 2] = length(unique(tmpPD$PATNO))
result[11, 3] = nrow(tmpPD)
result[11, 4] = length(unique(tmpGenPD$PATNO))
result[11, 5] = nrow(tmpGenPD)
result[11, 6] = length(unique(tmpRegPD$PATNO))
result[11, 7] = nrow(tmpRegPD)
result[11, 8] = length(unique(data$PATNO))
result[11, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$scopa.aut_odds),]
tmpPD = PD[!is.na(PD$scopa.aut_odds),]
tmpGenPD = GenPD[!is.na(GenPD$scopa.aut_odds),]
tmpRegPD = RegPD[!is.na(RegPD$scopa.aut_odds),]
result[12, 2] = length(unique(tmpPD$PATNO))
result[12, 3] = nrow(tmpPD)
result[12, 4] = length(unique(tmpGenPD$PATNO))
result[12, 5] = nrow(tmpGenPD)
result[12, 6] = length(unique(tmpRegPD$PATNO))
result[12, 7] = nrow(tmpRegPD)
result[12, 8] = length(unique(data$PATNO))
result[12, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$moca_odds),]
tmpPD = PD[!is.na(PD$moca_odds),]
tmpGenPD = GenPD[!is.na(GenPD$moca_odds),]
tmpRegPD = RegPD[!is.na(RegPD$moca_odds),]
result[13, 2] = length(unique(tmpPD$PATNO))
result[13, 3] = nrow(tmpPD)
result[13, 4] = length(unique(tmpGenPD$PATNO))
result[13, 5] = nrow(tmpGenPD)
result[13, 6] = length(unique(tmpRegPD$PATNO))
result[13, 7] = nrow(tmpRegPD)
result[13, 8] = length(unique(data$PATNO))
result[13, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$ess_odds),]
tmpPD = PD[!is.na(PD$ess_odds),]
tmpGenPD = GenPD[!is.na(GenPD$ess_odds),]
tmpRegPD = RegPD[!is.na(RegPD$ess_odds),]
result[14, 2] = length(unique(tmpPD$PATNO))
result[14, 3] = nrow(tmpPD)
result[14, 4] = length(unique(tmpGenPD$PATNO))
result[14, 5] = nrow(tmpGenPD)
result[14, 6] = length(unique(tmpRegPD$PATNO))
result[14, 7] = nrow(tmpRegPD)
result[14, 8] = length(unique(data$PATNO))
result[14, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$a.syn),]
tmpPD = PD[!is.na(PD$a.syn),]
tmpGenPD = GenPD[!is.na(GenPD$a.syn),]
tmpRegPD = RegPD[!is.na(RegPD$a.syn),]
result[15, 2] = length(unique(tmpPD$PATNO))
result[15, 3] = nrow(tmpPD)
result[15, 4] = length(unique(tmpGenPD$PATNO))
result[15, 5] = nrow(tmpGenPD)
result[15, 6] = length(unique(tmpRegPD$PATNO))
result[15, 7] = nrow(tmpRegPD)
result[15, 8] = length(unique(data$PATNO))
result[15, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$abeta),]
tmpPD = PD[!is.na(PD$abeta),]
tmpGenPD = GenPD[!is.na(GenPD$abeta),]
tmpRegPD = RegPD[!is.na(RegPD$abeta),]
result[16, 2] = length(unique(tmpPD$PATNO))
result[16, 3] = nrow(tmpPD)
result[16, 4] = length(unique(tmpGenPD$PATNO))
result[16, 5] = nrow(tmpGenPD)
result[16, 6] = length(unique(tmpRegPD$PATNO))
result[16, 7] = nrow(tmpRegPD)
result[16, 8] = length(unique(data$PATNO))
result[16, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$ptau),]
tmpPD = PD[!is.na(PD$ptau),]
tmpGenPD = GenPD[!is.na(GenPD$ptau),]
tmpRegPD = RegPD[!is.na(RegPD$ptau),]
result[17, 2] = length(unique(tmpPD$PATNO))
result[17, 3] = nrow(tmpPD)
result[17, 4] = length(unique(tmpGenPD$PATNO))
result[17, 5] = nrow(tmpGenPD)
result[17, 6] = length(unique(tmpRegPD$PATNO))
result[17, 7] = nrow(tmpRegPD)
result[17, 8] = length(unique(data$PATNO))
result[17, 9] = nrow(tmpdata)

tmpdata = data[!is.na(data$ttau),]
tmpPD = PD[!is.na(PD$ttau),]
tmpGenPD = GenPD[!is.na(GenPD$ttau),]
tmpRegPD = RegPD[!is.na(RegPD$ttau),]
result[18, 2] = length(unique(tmpPD$PATNO))
result[18, 3] = nrow(tmpPD)
result[18, 4] = length(unique(tmpGenPD$PATNO))
result[18, 5] = nrow(tmpGenPD)
result[18, 6] = length(unique(tmpRegPD$PATNO))
result[18, 7] = nrow(tmpRegPD)
result[18, 8] = length(unique(data$PATNO))
result[18, 9] = nrow(tmpdata)

write.csv(result, "../data/summary.csv", row.names=FALSE)
