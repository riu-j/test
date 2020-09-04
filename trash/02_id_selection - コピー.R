#library
library("dplyr")

#make all patients idlist
data <- read.csv("../data/csv/patient_status.csv", na.strings="")
data <- data[!is.na(data$ENROLL_CAT), c("PATNO", "ENROLL_CAT", "ENROLL_STATUS", "DESCRP_CAT")]
enroll_cat <- unique(data$ENROLL_CAT)

#output of idlist of each enroll category
for (i in enroll_cat){
  tmpid = unique(data[data$ENROLL_CAT==i, 1])
  output = data.frame("PATNO"=tmpid)
  write.csv(output, paste("../data/idlist/ids_", i, ".csv", sep=""), row.names=FALSE)
}

#select target ids
st <- read.csv("../data/csv/ST_CATALOG.csv", na.strings="", stringsAsFactors=FALSE)[,c("PATNO", "STRPLCVS")] %>%
  na.omit
colnames(st) <- c("PATNO", "STvisit")
StConverter <- function(df){
  target <- grep("ST", df$EVENT_ID)
  for (i in target){
    id <- df[i, "PATNO"]
    df[i,"EVENT_ID"] <- st[st$PATNO==id, "STvisit"]
  }
  return(df)
}

target_cat <- c("PD", "GENPD", "REGPD")
target_visit <- c("SC", "V04", "V06", "V10", "ST")
ids <- data.frame(NA)[0, ]
for (i in target_cat){
  tmp <- read.csv(paste("../data/idlist/ids_", i, ".csv", sep=""))
  ids <- rbind(ids, tmp)
}

dat <- read.csv("../data/csv/DATScan_Analysis.csv", stringsAsFactors=FALSE)[, c("PATNO", "EVENT_ID", "mean_striatum")] %>%
  filter(EVENT_ID %in% target_visit)
colnames(dat) <- c("PATNO", "EVENT_ID", "DaT scan")
dat <- StConverter(dat)

demo <- read.csv("../data/csv/Screening___Demographics.csv")[, c("PATNO", "APPRDX", "GENDER")]
rand <- read.csv("../data/csv/Randomization_table.csv", na.strings=".")[, c("PATNO", "age")]
info <- na.omit(merge(demo, rand, by="PATNO")); rm(demo, rand)
colnames(info) <- c("PATNO", "Diagnosis", "Sex", "age")

gene <- read.csv("../data/doc/PPMI_PD_Variants_Genetic_Status_WGS_20180921.csv") %>% 
  select(matches("PATNO|LRRK2|SNCA")) %>% na.omit

target_ids <- data.frame("PATNO" = dat$PATNO %>% intersect(gene$PATNO) %>% intersect(ids$PATNO)
                         %>% intersect(info$PATNO) %>% sort)

write.csv(target_ids, "../data/idlist/target_ids.csv", row.names=FALSE)
