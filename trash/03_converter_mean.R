library("tidyverse")

baseline <- 1.89
GTA <- 5
projectid <- "s131"
dv_names <- c("DaT scan", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "MDS-UPDRS Part3", "SCOPA-AUT")
cov_names <- c("Diagnosis", "Sex", "age", "rs356181", "rs3910105", "rs76904798", "G2019S")

data <- read.csv("../data/dataset.csv", na.strings=".")
cov <- data[!duplicated(data$ID), c("ID", cov_names)] %>% mutate(ID2=row_number())
initial <- data.frame(row.names=c("α", "β", "γ", "ave", "sd", "slope", "intercept", "meanslope"))

for (i in 1:length(dv_names)){
  tmp <- data[, c(1, 2, (2*i+1), (2*i+2))]
  colnames(tmp) <- c("ID", "TIME", "DVrow", "DV")
  tmp$lognorm <- scale(log(tmp$DV))
  tmp <- na.omit(tmp)
  tmp_slope <- tmp %>% group_by(ID) %>% do(model=lm(lognorm ~ TIME, data=.)) %>%
    do(data.frame(ID=.$ID, slope=coef(.$model)[["TIME"]]))
  tmp_x <- aggregate(TIME~ID, data=tmp, FUN=mean)
  tmp_y <- aggregate(lognorm~ID, data=tmp, FUN=mean)
  tmp_c <- aggregate(lognorm~ID, data=tmp, FUN=length)
  
  initial[4, i] <- tmp$DV %>% na.omit %>% log() %>% mean()
  initial[5, i] <- tmp$DV %>% na.omit %>% log() %>% sd()
  tmp_model <- lm(slope~V1, data=data.frame(tmp_slope, tmp_y))
  initial[6, i] <- tmp_model$coefficients[["V1"]]
  initial[7, i] <- tmp_model$coefficients[["(Intercept)"]]
  initial[8, i] <- tmp_slope$slope %>% na.omit %>% mean()
  if(i==6){
    initial[1, i] <- (log(baseline)-initial[4, i])/initial[5, i]
    initial[2, i] <- initial[7, i] + initial[1, i]*initial[6, i]
    initial[3, i] <- initial[6, i]
  }else{
    initial[2, i] <- initial[8, i]/exp(initial[6, i]*GTA)
    initial[3, i] <- initial[6, i]
    initial[1, i] <- (initial[2, i]-initial[7, i])/initial[6, i]
  }
  colnames(tmp_x) <- c("ID", paste("meanx", i, sep=""))
  colnames(tmp_y) <- c("ID", paste("meany", i, sep=""))
  colnames(tmp_c) <- c("ID", paste("count", i, sep=""))
  tmp$BM <- i
  if(i>1){
    meanx <- merge(meanx, tmp_x, by="ID", all=TRUE)
    meany <- merge(meany, tmp_y, by="ID", all=TRUE)
    counts <- merge(counts, tmp_c, by="ID", all=TRUE)
    df <- rbind(df, tmp)
  }else{
    meanx <- tmp_x
    meany <- tmp_y
    counts <- tmp_c
    df <- tmp
  }
}

colnames(initial) <- dv_names

output <- df %>% select(-DV) %>% rename(DV=lognorm) %>% na.omit %>% merge(meanx, on="ID", all = TRUE) %>%
  merge(meany, on="ID", all = TRUE) %>% merge(counts, on="ID", all = TRUE) %>% merge(cov, on="ID", all = TRUE) %>%
  select(ID, ID2, TIME, DV, DVrow, BM, everything())

output <- output[order(output$ID, output$BM, output$TIME), ]

if(!file.exists(paste("../notes/", projectid, "/", sep=""))){
  dir.create(paste("../notes/", projectid, "/", sep=""))
}

write.csv(output, paste("../notes/", projectid, "/", projectid, "_data.csv", sep=""), quote=FALSE, row.names=FALSE, na=".")
write.csv(initial, paste("../notes/", projectid, "/", projectid, "_initialprms.csv", sep=""), na=".")
