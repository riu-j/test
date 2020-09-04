# OverView
# This code carry on selection by inclusion criteria. Commented out code is for make id list, called only once.
# 
# Input
# dataset.csv
# id_(target categories).csv
# 
# Output
# selected_dataset.csv
# selected_id.csv
# selected_gene_info.csv

# library
library("dplyr")

# # make all patients idlist
# data <- read.csv("../data/raw/patient_status.csv", na.strings="")
# data <- data[!is.na(data$ENROLL_CAT), c("PATNO", "ENROLL_CAT", "ENROLL_STATUS", "DESCRP_CAT")]
# enroll_cat <- unique(data$ENROLL_CAT)
# 
# # output of idlist of each enroll category
# for (i in enroll_cat){
#   tmpid <- unique(data[data$ENROLL_CAT==i, 1])
#   output <- data.frame("PATNO"=tmpid)
#   write.csv(output, paste("../data/processed/id/id_", i, ".csv", sep=""), row.names=FALSE)
# }


process_time <- "200903_235046"

# select target subject
target_cat <- c("PD", "GENPD", "REGPD")
ids <- data.frame(NA)[0, ]
for (i in target_cat){
  tmp <- read.csv(paste("../data/processed/id/id_", i, ".csv", sep=""))
  ids <- rbind(ids, tmp)
}

# load dataset
df <- read.csv(paste("../data/processed/", process_time, "_dataset.csv", sep=""), na.strings=".")

# define selection criteria
hoge <- df[, c("ID", "α_synuclein")] %>% na.omit
hoge_id <- unique(hoge$ID)

# set target id
target_ids <- data.frame("PATNO" = hoge_id %>% intersect(ids$PATNO) %>% sort)

# output target ids
write.csv(target_ids, paste("../data/processed/", process_time, "_id_target.csv", sep=""), quote=FALSE, row.names=FALSE)

# output df selected by target ids
df_target <- filter(df, ID %in% target_ids$PATNO)
write.csv(df_target, paste("../data/processed/", process_time, "_dataset_target.csv", sep=""), quote=FALSE, row.names=FALSE, na=".")

# output gene infomation selected by target ids
gene <- read.csv(paste("../data/processed/", process_time, "_gene.csv", sep=""), na.strings=".") %>% filter(PATNO %in% target_ids$PATNO)
write.csv(gene, paste("../data/processed/", process_time, "_gene_target.csv", sep=""), quote=FALSE, row.names=FALSE, na=".")
