library("tidyverse")

#input data and make dataframe
{
  projectid <- "s132"
  folder_name <- "base_3045.046/" #ex) "base/"
  file_name <- "plots"

  folder_path <- paste("../notes/", projectid, "/", folder_name, projectid, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT <- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  iteration <- read.csv(paste(folder_path, "iteration.csv", sep=""), stringsAsFactors=FALSE, sep="", header=TRUE, skip=1)
  prms <- iteration[iteration$ITERATION==-1000000000, ] %>% as.numeric()
  numbm <- 3

  df <- merge(table, offsetT, by="ID", all=T)
  # df$TIME <- df$TIME + df$offsetT
  xrange <- c(min(df$TIME), max(df$TIME))
  alpha <- as.vector(prms[2:(1 + numbm)], mode="numeric")
  beta  <- as.vector(prms[(2 + numbm):(1 + 2 * numbm)], mode="numeric")
  gamma <- as.vector(prms[(2 + 2 * numbm):(1 + 3 * numbm)], mode="numeric")
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
  svg(paste(folder_path, file_name, "short.svg", sep=""), height=8, width=12, bg="transparent")
  par(mfrow=c(2, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  for (j in 1:numbm){
    data <- df %>% filter(BM==j)
    yrange <- c(min(data$DV), max(data$DV))
    # if (j==6){
    #   yrange[2] <- 5000
    # }
    # plot(data$row, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange, 
    #      type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    # if(xmax[j]==-1){
    #   mtext(expression(paste("(", x10^3, "pg/ml", ")")), adj=0.7, cex=0.8, at=xrange[1])
    # }
    #各観測値を描く
    # for (i in unique(data$ID)) {
    #   dsub <- data[data$ID==i,]
    #   lines(row~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=rgb(0, 0, 0, alpha=0.5))
    # }
    for(i in 1:3){
      df$BM <- replace(df$BM, which(df$BM==i), paste("Biomarker", i, sep=""))
    }
    ggplot(df, aes(x=TIME, y=DV, group=ID, color=BM))+
      geom_line(alpha=0.5) +
      geom_point(size=0.8, alpha=0.5) +
      guides(color=FALSE) +
      facet_wrap(~ BM, scales="free")
    myfunc <- function(x){
      y <- alpha + beta/gamma(exp(gamma*x)-1)
      return(y)
    }
    # if (j==6){
    #   plot(graph_main, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    # }
    # if(xmax[j]==-1){
    #   plot(graph_sub_2, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    # }
    # if (xmax[j]>0)
    if (j==1){
      plot(graph_sub, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
    if (j!=1 && j!=6){
      plot(graph_sub_score, xlim=xrange, ylim=yrange, col="black", lwd=4, add=TRUE)
    }
  }
  dev.off()
}
