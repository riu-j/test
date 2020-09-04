#set parameters
{
  proj = "s128"
  folder_name = "G2019S_dydx/"
  dv_name = c("DaT scan (SBR)", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "MDS-UPDRS Part3", "SCOPA-AUT", "α-synuclein", "Aβ1-42")

  folder_path = paste("../notes/", proj, "/", folder_name, proj,  "_", sep="")
  table = read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT = read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)

  df = merge(table, offsetT, by.x="ID", by.y="id", all=T)
  df$TIME = df$TIME + df$offsetT
  df$onset = df$age-df$offsetT
  fun = function(x){x}
  fun2 = function(x){x*0}
  numbm = length(dv_name)
}

#DV vs CIPRED
{
  file_name = "vpc_CIPRED"
  svg(paste(folder_path, file_name, ".svg", sep=""), width=16, height=8, bg="transparent")
  par(mfrow=c(2, 4), mar=c(5.1, 6.1, 4.1, 2.1))
  for (i in 1:numbm){
    df_a = df[df$BM==i, ]
    # df_a = df_a[df_a$LRRK2==0, ]
    min_val = min(min(df_a$DV), min(df_a$CIPRED))
    max_val = max(max(df_a$DV), max(df_a$CIPRED))
    
    plot(df_a$DV, df_a$CIPRED, main=dv_name[i], xlab="DV", ylab="CIPRED", xlim=c(min_val, max_val)
         , ylim=c(min_val, max_val), cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    plot(fun, min_val, max_val, add=TRUE, lwd=3, lty="dotted")
    lm_cipred = lm(CIPRED~DV, data=df_a)
    abline(lm_cipred, col="red", lwd=3, lty="dotted")
  }
  dev.off()
}

#TIME vs CWRES
{
  file_name = "vpc_CWRES"
  svg(paste(folder_path, file_name, ".svg", sep=""), width=16, height=8, bg="transparent")
  par(mfrow=c(2, 4), mar=c(5.1, 6.1, 4.1, 2.1))
  for (i in 1:numbm){
    df_a = df[df$BM==i, ]
    # df_a = df_a[df_a$LRRK2==0, ]
    
    plot(df_a$TIME, df_a$CWRES, main=dv_name[i], xlab="TIME (year)", ylab="CWRES"
         , cex.lab=2.5, cex.main=3, cex.axis=2.5, tcl=0.5, las=1, lab=c(5, 5, 3))
    plot(fun2, min(df_a$TIME), max(df_a$TIME), add=TRUE, lwd=3, lty="dotted")
    lm_cwres = lm(CWRES~TIME, data=df_a)
    abline(lm_cwres, col="red", lwd=3, lty="dotted")
  }
  dev.off()
}

#offsetT vs onset
{
  file_name = "vpc_onset"
  svg(paste(folder_path, file_name, ".svg", sep=""), width=7, height=5, bg="transparent")
  par(mar=c(5.1, 5.1, 2.1, 2.1), mgp=c(3.5, 1, 0))
  # df_a = df[df$LRRK2==0, ]
  df_a = df[!duplicated(df$ID),]
  plot(df_a$offsetT, df_a$onset,xlab="offsetT", ylab="onset age", cex=0.8,
       cex.lab=2, cex.axis=2, tcl=0.5, las=1, lab=c(5, 5, 3))
  lm_onset = lm(onset~offsetT, data=df)
  abline(lm_onset, col="red", lwd=3, lty="dotted")
  dev.off()
}

#offsetT vs onset (contour)
{
  library(MASS)
  file_name="vpc_contour_onset" 
  svg(paste(folder_path, file_name, ".svg", sep=""), width=7, height=5, bg="transparent")
  par(mar=c(5.1, 5.1, 2.1, 2.1), mgp=c(3.5, 1, 0))
  # df_a = df[df$LRRK2==1, ]
  df_a = df[!duplicated(df$id),]
  image(kde2d(df_a$offsetT, df_a$onset,n=80),xlab="offsetT", ylab="onset age", cex=0.8,
        cex.lab=2, cex.axis=2, tcl=0.5, las=1, lab=c(5, 5, 3))
  dev.off()
}

