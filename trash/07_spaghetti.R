#input data and make dataframe
{
  proj = "s115"
  file_name = "before_long"
  dv_name = c("DaT scan (SBR)", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "SCOPA-AUT", "CSF α-synuclein")
  cov_name = c("Sex", "LRRK2", "SNCA1", "SNCA2", "NHY", "Diag")
  target_cov = "NHY"
  xmax = c(0, 53, 53, 70, -1)
                          
  folder_path = paste("../note/", proj, "/", proj,  "_", sep="")
  data_before = 
    read.csv(paste(folder_path, "data_before.csv", sep=""), stringsAsFactors=FALSE, na.strings =".")
  offsetT = read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  numbm = length(dv_name)
  colnames(data_before)=c("id", "TIME", paste("BM", seq(from=1, to=numbm), sep=""), cov_name)
  
  df = merge(data_before, offsetT, by="id", all=T)
  df$TIME2 = df$TIME+df$offsetT

  xrange = c(min(df$TIME2), max(df$TIME2))
  color = c("blue", "green", "red")
  dev1 = 1.8
  dev2 = 2.05
  cex_size = 0.6
  lwd_size = 0.6
}

{
  svg(paste(folder_path, file_name, ".svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
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