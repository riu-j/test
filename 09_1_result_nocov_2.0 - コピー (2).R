library("tidyverse")
library("RColorBrewer")

#input data and make dataframe
{
  projectid <- "s136"
  folder_name <- "analytical/" #ex) "base/"
  file_name <- "plots"
  
  ans <- read.csv(paste("../notes/", projectid, "/", "build/", projectid, "_", "answer.csv", sep=""), stringsAsFactors=FALSE)
  colnames(ans) <- c("ID", "answer")

  folder_path <- paste("../notes/", projectid, "/", folder_name, projectid, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT <- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  iteration <- read.csv(paste(folder_path, "iteration.csv", sep=""), stringsAsFactors=FALSE, sep="", header=TRUE, skip=1)
  dv_name <- c("Biomarker 1", "Biomarker 2", "Biomarker 3")
  prms <- iteration[iteration$ITERATION==-1000000000, ] %>% as.numeric()
  numbm <- length(dv_name)

  df <- merge(table, offsetT, by="ID", all=T)
  df <- merge(df, ans, by="ID", all=T)
  df$TIME2 <- df$TIME + df$answer
  df$TIME <- df$TIME + df$offsetT
  xrange <- c(min(df$TIME), max(df$TIME))

  cex_size <- 0.6
  lwd_size <- 0.6
  
  cols <- brewer.pal(3, "Set1")
}

#define functions
{
  #for main biomarker
  graph_main <- function(x){return(alpha[j]+beta[j]/gamma[j]*(exp(gamma[j]*x)-1))}
}

##no covariate model using one color
{
  svg(paste(folder_path, file_name, "_before.svg", sep=""), height=4, width=12, bg="transparent")
  par(mfrow=c(1, 3), mar=c(5.1, 5.1, 4.1, 2.1))
  alpha <- c(1.5, -2, 2)
  beta <- c(-0.2, 0.05, -0.1)
  gamma <- c(-0.1, 0.1, 0.1)
  for (j in 1:numbm){
    data <- df %>% filter(BM==j)
    yrange <- c(min(data$DV), max(data$DV))
    plot(data$DV, data$TIME2, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange,
         type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    for (i in unique(data$ID)) {
      dsub <- data[data$ID==i,]
      lines(DV~TIME2, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=cols[j])
    }
    plot(graph_main, xlim=xrange, ylim=yrange, col=cols[j], lwd=4, add=TRUE)
  }
  # alpha <- as.vector(prms[2:(1 + numbm)], mode="numeric")
  # beta  <- as.vector(prms[(2 + numbm):(1 + 2 * numbm)], mode="numeric")
  # gamma <- as.vector(prms[(2 + 2 * numbm):(1 + 3 * numbm)], mode="numeric")
  # for (j in 1:numbm){
  #   data <- df %>% filter(BM==j)
  #   yrange <- c(min(data$DV), max(data$DV))
  #   plot(data$DV, data$TIME, main =dv_name[j], xlab="Time (year)", ylab="", xlim=xrange, ylim=yrange,
  #        type="n", cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
  #   for (i in unique(data$ID)) {
  #     dsub <- data[data$ID==i,]
  #     lines(DV~TIME, dsub, cex=cex_size, lwd=lwd_size, pch=16, type="o", col=cols[j])
  #   }
  #   plot(graph_main, xlim=xrange, ylim=yrange, col=cols[j], lwd=4, add=TRUE)
  # }
  dev.off()
}
