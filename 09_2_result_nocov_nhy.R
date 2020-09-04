#input data and make dataframe
{
  proj = "s132"
  folder_name = "base_3045.046/"
  file_name = "plots"
  dv_name = c("DaT scan (SBR)", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "SCOPA-AUT", "CSF α-synuclein")#, "MDS-UPDRS Part3", "Aβ1-42"
  # cov_name = c("G2019S", "Diag", "Sex", "age")
  cov_name = c("NHY", "rs769", "G2019S", "rs356", "rs391", "Diag", "Sex", "age")
  xmax = c(0, 52, 52, 70, -1)#, 132, -1

  folder_path = paste("../notes/", proj, "/", folder_name, proj,  "_", sep="")
  data_before = 
    read.csv(paste(folder_path, "data_before.csv", sep=""), stringsAsFactors=FALSE, na.strings =".")
  offsetT = read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  info = read.csv(paste(folder_path, "info.csv", sep=""), stringsAsFactors=FALSE)
  iteration = read.csv(paste(folder_path, "iteration.csv", sep=""), sep = "", header = TRUE
                       , skip = 1, stringsAsFactors=FALSE)
  prms = iteration[iteration$ITERATION==-1000000000,]
  numbm = length(dv_name)
  colnames(data_before) = c("id", "TIME", paste("BM", seq(from=1, to=numbm), sep=""), cov_name)
  
  df = merge(data_before, offsetT, by="id", all=T)
  df$TIME = df$TIME + df$offsetT
  xrange = c(min(df$TIME), max(df$TIME))
  alpha = as.vector(prms[,2:(1+numbm)], mode="numeric")
  beta  = as.vector(prms[,(2+numbm):(1+2*numbm)], mode="numeric")
  gamma = as.vector(prms[,(2+2*numbm):(1+3*numbm)], mode="numeric")
  ave = info$ave
  sd  = info$sd
  c = gamma
  b = beta/gamma*sd
  a = (exp(alpha-beta/gamma)^sd)*exp(ave)
  # dt = (1 + as.numeric(prms[17]))
  dt = (1 + as.numeric(prms[17])) * (1 + as.numeric(prms[18])) * (1 + as.numeric(prms[19]))
  dy = rep(0, numbm)
  dev1 = 1.8
  dev2 = 2.05
  cex_size = 0.6
  lwd_size = 0.6
}

#define functions
{
  #for main biomarker
  graph_main = function(x){return(a[j]*exp(b[j]*exp(c[j]*dt*x))/exp(sd[j]*dy[j]*(exp(c[j]*dt*x)-1)))}
  graph_main_non = function(x){return(a[j]*exp(b[j]*exp(c[j]*x))/exp(sd[j]*0*(exp(c[j]*x)-1)))}
  
  #for other biomarkers
  graph_sub = function(x){return(a[j]*exp(b[j]*exp(c[j]*dt*x))*exp(sd[j]*dy[j]))}
  graph_sub_non = function(x){return(a[j]*exp(b[j]*exp(c[j]*1*x))*exp(sd[j]*0))}
  
  #for a-synuclein
  graph_sub_2 = function(x){return(a[j]*exp(b[j]*exp(c[j]*dt*x))*exp(sd[j]*dy[j])/1000)}
  graph_sub_2_non = function(x){return(a[j]*exp(b[j]*exp(c[j]*1*x))*exp(sd[j]*0)/1000)}
  
  #other score biomarker's function with covariate 
  graph_sub_score = function(x){return(graph_sub(x)/(1+graph_sub(x))*xmax[j]-0.5)}
  graph_sub_score_non = function(x){return(graph_sub_non(x)/(1+graph_sub_non(x))*xmax[j]-0.5)}
}

##no covariate model using one color
{
  svg(paste(folder_path, file_name, ".svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    list = c("id", "TIME", paste("BM", j, sep=""))
    data = as.data.frame(lapply(df[, list], as.numeric))
    data = na.omit(data)
    colnames(data) = c("id", "TIME", "DV")
    if(xmax[j]==-1){
      data$DV = data$DV/1000
    }
    if(xmax[j]>0){
      data$DV = round(data$DV/(data$DV+1)*xmax[j]-0.5)
    }
    yrange = c(min(data$DV), max(data$DV))
    plot(data$DV, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    if(xmax[j]==-1){
      mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    }
    #各観測値を描く
    for (i in unique(data$id)) {
      dsub = data[data$id==i,]
      lines(DV~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col="black")
    }
    
    if (xmax[j]==0){
      plot(graph_main, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if(xmax[j]==-1){
      plot(graph_sub_2, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if (xmax[j]>0){
      plot(graph_sub_score, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
  }
  dev.off()
}

##no covariate model using Hoehn&Yahr
{
  svg(paste(folder_path, file_name, ".svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    list = c("id", "TIME", paste("BM", j, sep=""))
    data = as.data.frame(lapply(df[, list], as.numeric))
    data = na.omit(data)
    colnames(data) = c("id", "TIME", "DV")
    if(xmax[j]==-1){
      data$DV = data$DV/1000
    }
    if(xmax[j]>0){
      data$DV = round(data$DV/(data$DV+1)*xmax[j]-0.5)
    }
    yrange = c(min(data$DV), max(data$DV))
    plot(data$DV, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    if(xmax[j]==-1){
      mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    }
    #各観測値を描く
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

    if (xmax[j]==0){
      plot(graph_main, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if(xmax[j]==-1){
      plot(graph_sub_2, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if (xmax[j]>0){
      plot(graph_sub_score, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
  }
  
  #add legend
  plot(graph_main, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  legend("bottomright"
         , legend=c(paste("meanH&Y<",dev1), paste(dev1, "=<meanH&Y<", dev2), paste(dev2, "=<meanH&Y"))
         , col=color, box.lwd = 1, lty=1, cex=2, lwd=1.5, pch=16)
  
  dev.off()
}

##covariate model
{
  svg(paste(folder_path, file_name, "_women.svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    list = c("id", "TIME", paste("BM", j, sep=""), target_cov, "Sex")
    data = as.data.frame(lapply(df[, list], as.numeric))
    data = na.omit(data)
    colnames(data) = c("id", "TIME", "DV", target_cov, "Sex")
    if(xmax[j]==-1){
      data$DV = data$DV/1000
    }
    if(xmax[j]>0){
      data$DV = round(data$DV/(data$DV+1)*xmax[j]-0.5)
    }
    yrange = c(min(data$DV), max(data$DV))
    data = data[data$Sex==1,]
    #draw empty plot
    plot(data$DV, data$TIME, xlim=xrange, ylim=yrange, xlab="Time (year)", ylab="", main =dv_name[j],
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    if(xmax[j]==-1){
      mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    }

    #draw each plot
    for (i in unique(data$id)) {
      dsub = data[data$id==i,]
      lines(DV~TIME, dsub,
            cex=cex_size, lwd=lwd_size, pch=16, type="o",
            col=switch(dsub[target_cov][1,]+1, color[1], color[2]))
    }

    #draw line
    if(j==1){
      plot(graph_main, lwd=4, xlim=xrange, ylim=yrange, col=color[2], add=TRUE)
      plot(graph_main_non, lwd=4, xlim=xrange, ylim=yrange, col=color[1], add=TRUE)
    }else if(xmax[j]==-1){
      plot(graph_sub_2, lwd=4, xlim=xrange, ylim=yrange, col=color[2], add=TRUE)
      plot(graph_sub_2_non, lwd=4, xlim=xrange, ylim=yrange, col=color[1], add=TRUE)
    }else{
      plot(graph_sub_score, lwd=4, xlim=xrange, ylim=yrange, col=color[2], add=TRUE)
      plot(graph_sub_score_non, lwd=4, xlim=xrange, ylim=yrange, col=color[1], add=TRUE)
    }
  }

  #add legend
  plot(0, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  legend("bottomright", legend=leg, col=color, box.lwd = 1, lty=1, cex=2, lwd=1.5, pch=16)

  dev.off()
}
