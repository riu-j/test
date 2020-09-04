#Library
library("tidyverse")

#Function
StConverter <- function(df){
  target <- grep("ST", df$EVENT_ID)
  for (i in target){
    id <- df[i, "PATNO"]
    replace <- st[st$PATNO==id,2]
    if(length(replace)!=0){df[i,"EVENT_ID"] <- replace}
  }
  del <- grep("ST",df$EVENT_ID)
  if (length(del)!=0){df <- df[-del,]}
  return(df)
}

LoadCsv <- function(df){
  output <- df %>% filter(EVENT_ID %in% target_visit) %>% filter(PATNO %in% target_ids) %>% StConverter %>% na.omit
  return(output)
}

ScoreConverter <- function(df, max){
  output <- df
  output[, 4] <- output[, 3]
  if(max != 0){
    output[, 4] <- (df[, 3] + 0.5) / (max + 1)
    output[, 4] <- output[, 4] / (1 - output[, 4])
  }
  return(output)
}

Renamer <- function(df, name){
  colnames(df) <- c("PATNO", "EVENT_ID", paste(name, "_row", sep=""), name)
  return(df)
}

merge2 <- function(dfs, ...){
  base <- dfs[1]
  lapply(dfs[-1], function(i) base <<- merge(base, i, ...))
  return(base)
}

#set parameters
target_visit <- c("SC", "ST", "BL", "V01", "V02", "V03", "V04", "V05", "V06", "V07", 
                 "V08", "V09", "V10", "V11", "V12", "V13", "V14", "V15", "V16")
target_ids <- read.csv("../data/idlist/target_ids.csv")[, 1]

#St visit
st <- read.csv("../data/csv/ST_CATALOG.csv", na.strings="", stringsAsFactors=FALSE)[,c("PATNO", "STRPLCVS")] %>% na.omit
colnames(st) <- c("PATNO", "STvisit")

#For baseline information and gene_012 output for HWE test
demo <- read.csv("../data/csv/Screening___Demographics.csv")[, c("PATNO", "APPRDX", "GENDER")] %>% 
  filter(PATNO %in% target_ids)
rand <- read.csv("../data/csv/Randomization_table.csv", na.strings=".")[, c("PATNO", "age")] %>% 
  filter(PATNO %in% target_ids)
info <- merge(demo, rand, by="PATNO") %>% na.omit; rm(demo, rand)
colnames(info) <- c("PATNO", "Diagnosis", "Sex", "age")

gene_012 <- read.csv("../data/doc/PPMI_PD_Variants_Genetic_Status_WGS_20180921.csv") %>% 
  filter(PATNO %in% target_ids) %>% select(matches("PATNO|LRRK2|SNCA")) %>% na.omit
colnames(gene_012) <- c("PATNO", "rs356181", "rs3910105", "A53T", "E46K", "A30P", "rs76904798", "R1441G", "R1441C", 
                       "R1628P.H", "Y1699C", "G2019S", "G2385R")
# write.csv(gene_012, "../data/gene.csv", row.names=FALSE)
gene_01 <- replace(gene_012, gene_012==2, 1)

#For DaT scan
dat <- read.csv("../data/csv/DATScan_Analysis.csv", 
                stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mean_striatum")] %>% LoadCsv %>% ScoreConverter(0) %>% Renamer("DaT scan")

#MDS-UPDRS Part1
mds1a <- read.csv("../data/csv/MDS_UPDRS_Part_I.csv", 
                stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs1_sum")] %>% LoadCsv
mds1b <- read.csv("../data/csv/MDS_UPDRS_Part_I__Patient_Questionnaire.csv", 
                 stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs1Q_sum")] %>% LoadCsv
mds1 <- merge(mds1a, mds1b); rm(mds1a, mds1b)
mds1[, 5] <- mds1[, 3] + mds1[, 4]
mds1 <- mds1[,c(1, 2, 5)] %>% ScoreConverter(52) %>% Renamer("MDS-UPDRS Part1")

#MDS=UPDRS Part2
mds2 <- read.csv("../data/csv/MDS_UPDRS_Part_II__Patient_Questionnaire.csv", 
                stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs2_total")] %>% LoadCsv %>% ScoreConverter(52) %>% Renamer("MDS-UPDRS Part2")

#MDS=UPDRS Part3
mds3 <- read.csv('../data/csv/MDS_UPDRS_Part_III.csv', 
                stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "PAG_NAME", "mds.updrs3_total")]
mds3 <- mds3[mds3$PAG_NAME=="NUPDRS3", -3] %>% LoadCsv %>% ScoreConverter(132) %>% Renamer("MDS-UPDRS Part3")

#Hoehn & Yahr
nhy <- read.csv('../data/csv/MDS_UPDRS_Part_III.csv', 
                 stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "PAG_NAME", "NHY")]
nhy <- nhy[nhy$PAG_NAME=="NUPDRS3", -3] %>% LoadCsv %>% ScoreConverter(132) %>% Renamer("NHY")

#SCOAP-AUT
scopa <- read.csv('../data/csv/SCOPA-AUT.csv', 
                 stringsAsFactors=FALSE)[, c('PATNO', 'EVENT_ID', 'scopa.aut_total')] %>% LoadCsv %>% ScoreConverter(69) %>% Renamer("SCOAP-AUT")

#MoCA (Montreal Cognitive Assessment)
moca <- read.csv('../data/csv/Montreal_Cognitive_Assessment__MoCA_.csv', 
               stringsAsFactors=FALSE)[, c('PATNO', 'EVENT_ID', 'MCATOT')] %>% LoadCsv %>% ScoreConverter(30) %>% Renamer("MoCA")

#ESS (Epworth Sleepiness Scale)
ess <- read.csv('../data//csv/Epworth_Sleepiness_Scale.csv', 
               stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "ess_total")] %>% LoadCsv %>% ScoreConverter(24) %>% Renamer("ESS")

#CSF
bio <- read.csv('../data//csv/Current_Biospecimen_Analysis_Results.csv', 
               stringsAsFactors=FALSE)[, c("PATNO", "CLINICAL_EVENT", "TYPE", "TESTNAME", "TESTVALUE", "UNITS")] %>% rename(EVENT_ID=CLINICAL_EVENT)
csf <- bio[bio$TYPE=="Cerebrospinal Fluid", c("PATNO", "EVENT_ID", "TESTNAME", "TESTVALUE")]

asyn <- csf[csf$TESTNAME=="CSF Alpha-synuclein",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="N/A") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0) %>% Renamer("α-synuclein")

abeta <- csf[csf$TESTNAME=="ABeta 1-42",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<200") %>%
  filter(TESTVALUE!="") %>% mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0) %>% Renamer("Amyloid β")

ptau <- csf[csf$TESTNAME=="pTau",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<8") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0) %>% Renamer("pTau")

ttau <- csf[csf$TESTNAME=="tTau",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<80") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0) %>% Renamer("tTau")

prms <- list(dat, mds1, mds2, mds3, nhy, scopa, moca, ess, asyn, abeta, ptau, ttau)
df <- merge2(prms, by=c("PATNO", "EVENT_ID"), all=TRUE, sort=F)

covs <- list(df, info, gene_01)
df <- merge2(covs, all=TRUE)  

times <- c("SC"="-0.12", "BL"="0", "V01"="0.25", "V02"="0.5", "V03"="0.75", "V04"="1", "V05"="1.5", "V06"="2", "V07"="2.5", "V08"="3", "V09"="3.5", "V10"="4", "V11"="4.5", "V12"="5", "V13"="6", "V14"="7", "V15"="8", "V16"="9")
df <- mutate(df, EVENT_ID=str_replace_all(EVENT_ID, pattern=times)) %>% rename(ID=PATNO, TIME=EVENT_ID)

#output
write.csv(df, "../data/dataset_200802.csv", row.names=FALSE, na=".")
