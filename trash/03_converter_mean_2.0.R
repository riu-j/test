library("tidyverse")

baseline <- 1600
GTA <- 5
proj <- "s139"
today_time <- "200903_134457"
dv_names <- c("α_synuclein", "DaT__scan", "MDS_UPDRS__Part1", "MDS_UPDRS__Part2", "MDS_UPDRS__Part3", "SCOPA_AUT")
cov_names <- c("NHY", "Diagnosis", "Sex", "age", "rs356181", "rs3910105", "rs76904798", "G2019S")

dataset <- read.csv(paste("../data/processed/", today_time, "_dataset_target.csv", sep=""), na.strings=".")
cov <- dataset[, c("ID", "TIME", cov_names)]
cov$serial <- as.numeric(factor(cov$ID))
df <- dataset[, c("ID", "TIME", paste(dv_names, "_row", sep=""), dv_names)]
numbm <- length(dv_names)

MakeNONMEMSheet <- function(df){
  for (i in 1:numbm){
    tmp <- df[, c(1, 2, (i+2), (i+2+numbm))] %>% setNames(c("ID", "TIME", "DVrow", "DV")) %>% na.omit
    tmp$lognorm <- tmp$DV %>% log() %>% scale()
    tmp_slope <- tmp %>% group_by(ID) %>% do(model=lm(lognorm ~ TIME, data=.)) %>%
      do(data.frame(ID=.$ID, slope=coef(.$model)[["TIME"]]))
    tmp_x <- aggregate(TIME~ID, data=tmp, FUN=mean) %>% setNames(c("ID", paste("meanx", i, sep="")))
    tmp_y <- aggregate(lognorm~ID, data=tmp, FUN=mean) %>% setNames(c("ID", paste("meany", i, sep="")))
    tmp_c <- aggregate(lognorm~ID, data=tmp, FUN=length) %>% setNames(c("ID", paste("count", i, sep="")))
    tmp$BM <- i
    if(i>1){
      meanx <- merge(meanx, tmp_x, by="ID", all=TRUE)
      meany <- merge(meany, tmp_y, by="ID", all=TRUE)
      counts <- merge(counts, tmp_c, by="ID", all=TRUE)
      output <- rbind(output, tmp)
    }else{
      meanx <- tmp_x
      meany <- tmp_y
      counts <- tmp_c
      output <- tmp
    }
  }
  output <- output %>% select(-DV) %>% rename(DV=lognorm) %>% na.omit %>% merge(meanx, on="ID", all = TRUE) %>%
    merge(meany, on="ID", all = TRUE) %>% merge(counts, on="ID", all = TRUE) %>% merge(cov, on=c("ID", "TIME"), all = TRUE) %>%
    select(ID, serial, TIME, DV, DVrow, BM, everything())
  output <- output[order(output$ID, output$BM, output$TIME), ]
  return(output)
}

DefineInitialParameters <- function(df, cov){
  initial <- data.frame(row.names=c("α", "β", "γ", "dY", "dT", "ave", "sd", "slope", "intercept", "meanslope"))
  for (i in 1:numbm){
    if(cov==""){
      tmp <- df %>% select(1, 2, (i+2), (i+2+numbm)) %>% setNames(c("ID", "TIME", "DVrow", "DV")) %>% na.omit
      tmp$lognorm <- tmp$DV %>% log() %>% scale()
    }else{
      tmp <- df %>% select(1, 2, (i+2), (i+2+numbm), all_of(cov)) %>% setNames(c("ID", "TIME", "DVrow", "DV", cov)) %>% na.omit
      tmp$lognorm <- tmp$DV %>% log() %>% scale()
      tmp1 <- tmp %>% filter_(paste(cov, "==", 1))
      tmp <- tmp %>% filter_(paste(cov, "==", 0))
    }
    
    tmp_slope <- tmp %>% group_by(ID) %>% do(model=lm(lognorm ~ TIME, data=.)) %>%
      do(data.frame(ID=.$ID, slope=coef(.$model)[["TIME"]]))
    tmp_y <- aggregate(lognorm~ID, data=tmp, FUN=mean)

    initial["ave", i] <- tmp$DV %>% na.omit %>% log() %>% mean()
    initial["sd", i] <- tmp$DV %>% na.omit %>% log() %>% sd()
    tmp_model <- lm(slope~V1, data=data.frame(tmp_slope, tmp_y))
    initial["slope", i] <- tmp_model$coefficients[["V1"]]
    initial["intercept", i] <- tmp_model$coefficients[["(Intercept)"]]
    initial["meanslope", i] <- tmp_slope$slope %>% na.omit %>% mean()
    
    if(i==1){
      initial["α", i] <- (log(baseline)-initial["ave", i])/initial["sd", i]
      initial["β", i] <- initial["intercept", i] + initial["α", i]*initial["slope", i]
      initial["γ", i] <- initial["slope", i]
    }else{
      initial["β", i] <- initial["meanslope", i]/exp(initial["slope", i]*GTA)
      initial["γ", i] <- initial["slope", i]
      initial["α", i] <- (initial["β", i]-initial["intercept", i])/initial["slope", i]
    }
    
    if(cov!=""){
      print("calc cov prms")
      tmp1_slope <- tmp1 %>% group_by(ID) %>% do(model=lm(lognorm ~ TIME, data=.)) %>%
        do(data.frame(ID=.$ID, slope=coef(.$model)[["TIME"]]))
      tmp1_y <- aggregate(lognorm~ID, data=tmp1, FUN=mean)
      tmp1_model <- lm(slope~V1, data=data.frame(tmp1_slope, tmp1_y))
      initial["dT", i] <- tmp1_model$coefficients[["V1"]] / initial["γ", i] - 1
      initial["dY", i] <-
        (initial["β", i] * initial["dT", i] - tmp1_model$coefficients[["(Intercept)"]]) / tmp1_model$coefficients[["V1"]] - initial["α", i]
    }
  }
  return(initial)
}


data <- MakeNONMEMSheet(df)
initial <- DefineInitialParameters(df, "") %>% setNames(dv_names)

if(!file.exists(paste("../notes/", proj, sep=""))){
  dir.create(paste("../notes/", proj, sep=""))
}

write.csv(data, paste("../notes/", proj, "/", proj, "_data.csv", sep=""), quote=FALSE, row.names=FALSE, na=".")
write.csv(initial, paste("../notes/", proj, "/", proj, "_initialprms.csv", sep=""), na=".")
