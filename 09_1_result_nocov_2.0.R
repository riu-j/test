library("tidyverse")

#input data and make dataframe
{
  projectid <- "s132"
  folder_name <- "base_3045.046/" #ex) "base/"
  file_name <- "plots"
  xmax = c(0, 52, 52, 132, 70, -1)
  
  folder_path <- paste("../notes/", projectid, "/", folder_name, projectid, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""), na.strings=".")
  initialprms <- read.csv(paste(folder_path, "initialprms.csv", sep=""), stringsAsFactors=FALSE, row.names=1, na.strings =".", check.names=FALSE)
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
  cex_size <- 0.6
  lwd_size <- 0.6
}

#define functions
{
  #for main biomarker
  graph_main <- function(x){return(a[j]*exp(b[j]*exp(c[j]*x))/exp(sd[j]*0*(exp(c[j]*x)-1)))}
  
  #for other biomarkers
  graph_sub <- function(x){return(a[j]*exp(b[j]*exp(c[j]*1*x))*exp(sd[j]*0))}
  
  #for a-synuclein
  graph_sub_2 <- function(x){return(a[j]*exp(b[j]*exp(c[j]*1*x))*exp(sd[j]*0)/1000)}
  
  #other score biomarker's function with covariate 
  graph_sub_score <- function(x){return(graph_sub(x)/(1+graph_sub(x))*xmax[j]-0.5)}
}

##no covariate model using one color
{
  svg(paste(folder_path, file_name, ".svg", sep=""), height=8, width=12, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    data <- df %>% filter(BM==j)
    yrange <- c(min(data$row), max(data$row))
    plot(data$row, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    #各観測値を描く
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=rgb(0, 0, 0, alpha=0.7))
    }

    if (j==6){
      plot(graph_main, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if (j==1){
      plot(graph_sub, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if (j!=1 && j!=6){
      plot(graph_sub_score, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
  }
  dev.off()
}

##no covariate model using Hoehn&Yahr
{
  svg(paste(folder_path, file_name, "_nhy_short.svg", sep=""), width=12, height = 8, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  dev1 <- 1.8
  dev2 <- 2.05
  nhy_list <- df[!duplicated(paste(df$ID, df$TIME, sep=",")),] %>% select(ID, TIME, NHY)
  for (j in 1:numbm){
    data <- df %>% filter(BM==j)
    yrange <- c(min(data$row), max(data$row))
    plot(data$row, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    #各観測値を描く
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      HY_score = nhy_list[nhy_list$ID==i,"NHY"] %>% na.omit() %>%mean()
      if(HY_score<dev1){
        lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=rgb(0, 0, 1, alpha=0.7))
      }
    }
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      HY_score = nhy_list[nhy_list$ID==i,"NHY"] %>% na.omit() %>%mean()
      if(HY_score>=dev1 &&HY_score<dev2){
        lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=rgb(0, 1, 0, alpha=0.7))
      }
    }
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      HY_score = nhy_list[nhy_list$ID==i,"NHY"] %>% na.omit() %>%mean()
      if(HY_score>=dev2){
        lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=rgb(1, 0, 0, alpha=0.7))
      }
    }    

    # if (j==6){
    #   plot(graph_main, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    # }
    # if (j==1){
    #   plot(graph_sub, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    # }
    # if (j!=1 && j!=6){
    #   plot(graph_sub_score, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    # }
  }
  
  #add legend
  # plot(graph_main, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  # legend("bottomright"
  #        , legend=c(paste("meanH&Y<",dev1), paste(dev1, "=<meanH&Y<", dev2), paste(dev2, "=<meanH&Y"))
  #        , col=color, box.lwd = 1, lty=1, cex=2, lwd=1.5, pch=16)
  
  dev.off()
}
