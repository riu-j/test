library("tidyverse")

#Set parameters
{
  proj <- "s139"
  folder_path <- paste("../notes/", proj, "/", proj,  "_", sep="")
  data <- read.csv(paste(folder_path, "data.csv", sep=""), na.strings=".")
  init <- read.csv(paste(folder_path, "initialprms.csv", sep=""), row.names=1, na.strings=".")
  dv_name <- colnames(init) %>% gsub("__", " ", .) %>% gsub("_", "-", .)
  numbm <- length(dv_name)
  aa <- 0.05/10
}

MakeDateFrame <- function(data){
  df <- data.frame(ID=NA, TIME=NA) %>% na.omit
  for (i in 1:numbm){
    tmp_data <- data[data$BM==i, c("ID", "TIME", "DV")] %>% na.omit
    colnames(tmp_data) <- c("ID", "TIME", dv_name[i])
    df <- merge(df, tmp_data, by=c("ID", "TIME"), all=TRUE)
  }
  return(df)
}
df <- MakeDateFrame(data)

#Difine functions for pairs plot
{
  ##Upper panel function
  upperf <- function(x, y){
    usr <- par("usr")
    p <- cor.test(x, y)$p.value
    corr <- round(cor.test(x, y)$estimate, 3)
    
    if((p<=aa) && (abs(corr)>=0.9)){
      rect(usr[1], usr[3], usr[2], usr[4], col="coral")
      text((usr[1]+usr[2])/2, (usr[3]+usr[4])/2, lab=corr, cex = 3.0)
    }else{
      text((usr[1]+usr[2])/2, (usr[3]+usr[4])/2, lab=corr, cex = 3.0)
    } 
  }
  
  ##Diagonal panel function
  diagf <- function(x){ 
    usr <- par("usr")
    on.exit(par(usr))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$density
    par(usr = c(usr[1:2], 0, max(y)*1.5))
    rect(breaks[-nB], 0, breaks[-1], y)
    lines(density(na.omit(x)), col = "red", lwd = 3)
  }
  
  ##Lower panel function
  lowerf <- function(x, y){
    points(x, y, pch=16, cex=1)
    abline(lm(y ~ x), col="red", lwd=3)
  }
}

#Make svg file (plot)
{
  svg(paste(folder_path, "scatter_matrix.svg", sep=""), width=2*numbm, height=2*numbm, bg="transparent")
  pairs(df[,-c(1, 2)],  diag.panel=diagf, lower.panel = lowerf, upper.panel = upperf)
  dev.off()
}
