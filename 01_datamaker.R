# Overview
# This code carry on create dataset file. PPMI's csv file was processed and integrated.
# 
# Input
# PPMI's csv files
# 
# Output
# dataset.csv

# Library
library("tidyverse")

# Function
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
  output <- df %>% StConverter %>% na.omit
  return(output)
}

ScoreConverter <- function(df, max, name){
  output <- df
  if(max != 0){
    output[, 4] <- (output[, 3] + 0.5) / (max + 1)
    output[, 4] <- output[, 4] / (1 - output[, 4])
    colnames(output) <- c("PATNO", "EVENT_ID", paste(name, "_row", sep=""), name)
  }
  if(max == 0){
    output[, 4] <- output[, 3]
    colnames(output) <- c("PATNO", "EVENT_ID", paste(name, "_row", sep=""), name)
  }
  return(output)
}

merge2 <- function(dfs, ...){
  base <- dfs[1]
  lapply(dfs[-1], function(i) base <<- merge(base, i, ...))
  return(base)
}

# process time
process_time <-  Sys.time() %>% substring(3, 19) %>% gsub("-", "", .) %>% gsub(":", "", .) %>% gsub(" ", "_", .)
today <- strsplit(process_time, "_")[[1]][1]
time <- strsplit(process_time, "_")[[1]][2]

# ST visit
st <- read.csv("../data/raw/ST_CATALOG.csv", na.strings="", stringsAsFactors=FALSE)[,c("PATNO", "STRPLCVS")] %>% na.omit
colnames(st) <- c("PATNO", "STvisit")

# For baseline information
demo <- read.csv("../data/raw/Screening___Demographics.csv")[, c("PATNO", "APPRDX", "GENDER")]
rand <- read.csv("../data/raw/Randomization_table.csv", na.strings=".")[, c("PATNO", "age")]
info <- merge(demo, rand, by="PATNO") %>% na.omit; rm(demo, rand)
colnames(info) <- c("PATNO", "Diagnosis", "Sex", "age")

# Gene information
gene_012 <- read.csv("../data/raw/PPMI_PD_Variants_Genetic_Status_WGS_20180921.csv") %>% dplyr::select(matches("PATNO|LRRK2|SNCA")) %>% na.omit
colnames(gene_012) <- c("PATNO", "rs356181", "rs3910105", "A53T", "E46K", "A30P", "rs76904798", "R1441G", "R1441C", "R1628P.H", "Y1699C", "G2019S", "G2385R")
gene_01 <- replace(gene_012, gene_012==2, 1)

# DaT scan
dat <- read.csv("../data/raw/DATScan_Analysis.csv", stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mean_striatum")] %>% LoadCsv %>% ScoreConverter(0, "DaT__scan")

# MDS-UPDRS Part1
mds1a <- read.csv("../data/raw/MDS_UPDRS_Part_I.csv", stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs1_sum")] %>% LoadCsv
mds1b <- read.csv("../data/raw/MDS_UPDRS_Part_I__Patient_Questionnaire.csv", stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs1Q_sum")] %>% LoadCsv
mds1 <- merge(mds1a, mds1b); rm(mds1a, mds1b)
mds1[, 5] <- mds1[, 3] + mds1[, 4]
mds1 <- mds1[,c(1, 2, 5)] %>% ScoreConverter(52, "MDS_UPDRS__Part1")

# MDS=UPDRS Part2
mds2 <- read.csv("../data/raw/MDS_UPDRS_Part_II__Patient_Questionnaire.csv", stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mds.updrs2_total")] %>% LoadCsv %>% ScoreConverter(52, "MDS_UPDRS__Part2")

# MDS=UPDRS Part3
mds3 <- read.csv('../data/raw/MDS_UPDRS_Part_III.csv', stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "PAG_NAME", "mds.updrs3_total")]
mds3 <- mds3[mds3$PAG_NAME=="NUPDRS3", -3] %>% LoadCsv %>% ScoreConverter(132, "MDS_UPDRS__Part3")

# Hoehn & Yahr
nhy <- read.csv('../data/raw/MDS_UPDRS_Part_III.csv', stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "PAG_NAME", "NHY")]
nhy <- nhy[nhy$PAG_NAME=="NUPDRS3", -3] %>% LoadCsv %>% ScoreConverter(0, "NHY")

# SCOAP-AUT
scopa <- read.csv('../data/raw/SCOPA-AUT.csv', stringsAsFactors=FALSE)[, c('PATNO', 'EVENT_ID', 'scopa.aut_total')] %>% LoadCsv %>% ScoreConverter(69, "SCOPA_AUT")

# MoCA (Montreal Cognitive Assessment)
moca <- read.csv('../data/raw/Montreal_Cognitive_Assessment__MoCA_.csv', stringsAsFactors=FALSE)[, c('PATNO', 'EVENT_ID', 'MCATOT')] %>% LoadCsv %>% ScoreConverter(30, "MoCA")

# ESS (Epworth Sleepiness Scale)
ess <- read.csv('../data//raw/Epworth_Sleepiness_Scale.csv', stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "ess_total")] %>% LoadCsv %>% ScoreConverter(24, "ESS")

# Biospecimen
bio <- read.csv('../data//raw/Current_Biospecimen_Analysis_Results.csv', stringsAsFactors=FALSE)[, c("PATNO", "CLINICAL_EVENT", "TYPE", "TESTNAME", "TESTVALUE", "UNITS")] %>% rename(EVENT_ID=CLINICAL_EVENT)
# CerebroSpinal Fluid
csf <- bio[bio$TYPE=="Cerebrospinal Fluid", c("PATNO", "EVENT_ID", "TESTNAME", "TESTVALUE")]

# Alpha-synuclein
asyn <- csf[csf$TESTNAME=="CSF Alpha-synuclein",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="N/A") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0, "α_synuclein")

# Amyloid beta
abeta <- csf[csf$TESTNAME=="ABeta 1-42",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<200") %>%
  filter(TESTVALUE!="") %>% mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0, "Amyloid__β")

# p-tau
ptau <- csf[csf$TESTNAME=="pTau",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<8") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0, "pTau")

# t-tau
ttau <- csf[csf$TESTNAME=="tTau",  c("PATNO", "EVENT_ID", "TESTVALUE")] %>% filter(TESTVALUE!="<80") %>%
  mutate_at(vars(-PATNO, -EVENT_ID), as.numeric) %>% with(aggregate(TESTVALUE, list(PATNO, EVENT_ID), mean)) %>% rename(PATNO=Group.1, EVENT_ID=Group.2) %>% LoadCsv %>% ScoreConverter(0, "tTau")

# merge dv items
prms <- list(dat, mds1, mds2, mds3, nhy, scopa, moca, ess, asyn, abeta, ptau, ttau)
df <- merge2(prms, by=c("PATNO", "EVENT_ID"), all=TRUE, sort=F)

# merge covariates
covs <- list(df, info, gene_01)
df <- merge2(covs, by="PATNO")  

# convert event_id to time
times <- c("SC"="-0.12", "BL"="0", "V01"="0.25", "V02"="0.5", "V03"="0.75", "V04"="1", "V05"="1.5", "V06"="2", "V07"="2.5", "V08"="3", "V09"="3.5", "V10"="4", "V11"="4.5", "V12"="5", "V13"="6", "V14"="7", "V15"="8", "V16"="9")
df <- filter(df, !EVENT_ID %in% c("U01", "U02", "PW", "RS1"))
df$TIME <- df$EVENT_ID
df <- mutate(df, TIME=str_replace_all(TIME, pattern=times)) %>% rename(ID=PATNO)

# output
write.csv(gene_012, paste("../data/processed/", process_time, "_gene.csv", sep=""), row.names=FALSE)
write.csv(df, paste("../data/processed/", process_time, "_dataset.csv", sep=""), quote=FALSE, row.names=FALSE, na=".")