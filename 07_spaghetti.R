#input data and make dataframe
{
  proj <- "s139"
  folder_path <- paste("../notes/", proj, "/", proj,  "_", sep="")
  data <- read.csv(paste(folder_path, "data.csv", sep=""), na.strings=".")
  init <- read.csv(paste(folder_path, "initialprms.csv", sep=""), row.names=1, na.strings=".")
  dv_name <- colnames(init) %>% gsub("__", " ", .) %>% gsub("_", "-", .)
  numbm <- length(dv_name)
  cex_size <- 0.6
  lwd_size <- 0.6
  row <- sqrt(numbm) %>% trunc()
  col <- (numbm/row) %>% ceiling()
}

{
  svg(paste(folder_path, "spaghetti.svg", sep=""), width=col*3, height=row*4, bg="transparent")
  par(mfrow=c(row, col), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    df <- data[data$BM==j, c("ID", "TIME", "DVrow")] %>% na.omit
    xrange <- c(min(df$TIME), max(df$TIME))
    yrange <- c(min(df$DVrow), max(df$DVrow))
    plot(df$DVrow, df$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         lwd=1, type="n", cex.lab=2, cex.main=2.5, cex.axis=2, tcl=0.5, las=1, lab=c(5, 5, 3))

    #draw individual plots
    for (i in unique(df$ID)) {
      dsub <- df[df$ID==i,]
      lines(DVrow~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o")
    }
  }
  dev.off()
}

{
  svg(paste(folder_path, "spaghetti.svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(row, col), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    list = c("id", "TIME", paste("BM", j, sep=""), target_cov)
    data = as.data.frame(lapply(df[, list], as.numeric))
    data = na.omit(data)
    colnames(data) = c("id", "TIME", "DV", target_cov)
    if(xmax[j]==-1){
      data$DV = data$DV/1000
    }
    if(xmax[j]>0){
      data$DV = round(data$DV/(data$DV+1)*xmax[j]-0.5)
    }
    yrange = c(min(data$DV), max(data$DV))
    plot(data$DV, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         lwd=1, type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    if(xmax[j]==-1){
      mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    }
    
    #draw individual plots
    for (i in unique(data$id)) {
      dsub = data[data$id==i,]
      HY_score = mean(na.omit(df[df$id==i,]$NHY))
      if(HY_score<dev1){
        lines(DV~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col="blue")
      }
    }
    
    for (i in unique(data$id)) {
      dsub = data[data$id==i,]
      HY_score = mean(na.omit(df[df$id==i,]$NHY))
      if(HY_score>=dev1 &&HY_score<dev2){
        lines(DV~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col="green")
      }
    }
    
    for (i in unique(data$id)) {
      dsub = data[data$id==i,]
      HY_score = mean(na.omit(df[df$id==i,]$NHY))
      if(HY_score>=dev2){
        lines(DV~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col="red")
      }
    }
    
  }
  
  plot(0, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  legend("bottomright"
         , legend=c(paste("meanH&Y<",dev1), paste(dev1, "=<meanH$Y<", dev2), paste(dev2, "=<meanH&Y"))
         , col=color, box.lwd = 1, lty=1, cex=2, lwd=1.5, pch=16)
  
  dev.off()
}