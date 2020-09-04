library("tidyverse")

#input data and make dataframe
{
  projectid <- "s132"
  folder_name <- "G2019S_dT_numericdif_e8/" #ex) "base/"
  file_name <- "plots"
  target_cov <- "G2019S"
  xmax = c(0, 52, 52, 132, 70, -1)
  
  folder_path <- paste("../notes/", projectid, "/", folder_name, projectid, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""))
  initialprms <- read.csv(paste(folder_path, "initialprms.csv", sep=""), stringsAsFactors=FALSE, row.names=1, na.strings =".")
  offsetT <- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  iteration <- read.csv(paste(folder_path, "iteration.csv", sep=""), stringsAsFactors=FALSE, sep="", header=TRUE, skip=1)
  dv_name <- colnames(initialprms)
  cov_name <- table %>% select(-(1:7)) %>% colnames()
  prms <- iteration[iteration$ITERATION==-1000000000, ] %>% as.numeric()
  numbm <- length(dv_name)

  df <- merge(table, offsetT, by="ID", all=T)
  df$TIME <- df$TIME + df$offsetT
  xrange <- c(min(df$TIME), max(df$TIME))
  alpha <- as.vector(prms[2:(1 + numbm)], mode="numeric")
  beta  <- as.vector(prms[(2 + numbm):(1 + 2 * numbm)], mode="numeric")
  gamma <- as.vector(prms[(2 + 2 * numbm):(1 + 3 * numbm)], mode="numeric")
  ave <- initialprms["ave", ] %>% as.numeric()
  sd  <- initialprms["sd", ] %>% as.numeric()
  c <- gamma
  b <- beta / gamma * sd
  a <- (exp(alpha - beta / gamma) ^ sd) * exp(ave)
  dt <- 1 + as.numeric(prms[20])
  dy = numeric(numbm)
  cex_size <- 0.6
  lwd_size <- 0.6
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
  svg(paste(folder_path, file_name, ".svg", sep=""), height=8, width=12, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    data <- df %>% filter(BM==j)
    # if(xmax[j]==-1){
    #   data$DV <- data$DV/1000
    # }
    yrange <- c(min(data$row), max(data$row))
    # if (j==6){
    #   yrange[2] <- 5000
    # }
    plot(data$row, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    # if(xmax[j]==-1){
    #   mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    # }
    #各観測値を描く
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=switch(dsub[target_cov][1, ]+1, "blue", "red"))
    }
    # ggplot(data, aes(x=TIME, y=row, group=ID))+
    #   geom_line(alpha=0.3) +
    #   geom_point(size=0.8, alpha=0.3)
    
    #draw line
    if(j==6){
      plot(graph_main, lwd=4, xlim=xrange, ylim=yrange, col="red", add=TRUE)
      plot(graph_main_non, lwd=4, xlim=xrange, ylim=yrange, col="blue", add=TRUE)
    }else if(xmax[j]==0){
      plot(graph_sub, lwd=4, xlim=xrange, ylim=yrange, col="red", add=TRUE)
      plot(graph_sub_non, lwd=4, xlim=xrange, ylim=yrange, col="blue", add=TRUE)
    }else{
      plot(graph_sub_score, lwd=4, xlim=xrange, ylim=yrange, col="red", add=TRUE)
      plot(graph_sub_score_non, lwd=4, xlim=xrange, ylim=yrange, col="blue", add=TRUE)
    }

  }
  dev.off()
}
