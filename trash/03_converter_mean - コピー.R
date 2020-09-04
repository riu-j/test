library("dplyr")

data = read.csv("../data/data_trim.csv", na.strings=".")
cov = data[, c(1:2, 14:21)]

for (i in 1:9){
  tmp = data[,c(1, 2, (i+2))]
  colnames(tmp) = c("PATNO", "EVENT_ID", "DV")
  tmp$DV = scale(log(tmp$DV))
  # tmp = na.omit(tmp)
  tmp_x = aggregate(EVENT_ID~PATNO, data=na.omit(tmp), FUN=mean)
  colnames(tmp_x) = c("PATNO", paste("meanx", i, sep=""))
  if(i>1) meanx = merge(meanx, tmp_x, by="PATNO", all=TRUE)
  if(i==1) meanx = tmp_x

  tmp_y = aggregate(DV~PATNO, data=tmp, FUN=mean)
  colnames(tmp_y) = c("PATNO", paste("meany", i, sep=""))
  if(i>1) meany = merge(meany, tmp_y, by="PATNO", all=TRUE)
  if(i==1) meany = tmp_y
  
  tmp_c = aggregate(DV~PATNO, data=tmp, FUN=length)
  colnames(tmp_c) = c("PATNO", paste("count", i, sep=""))
  if(i>1) counts = merge(counts, tmp_c, by="PATNO", all=TRUE)
  if(i==1) counts = tmp_c
  
  tmp$BM = i
  if(i>1) df = rbind(df, tmp)
  if(i==1) df = tmp
}

df %>% na.omit %>% merge(meanx, on="PATNO") %>% merge(meany, on="PATNO") %>% merge(counts, on="PATNO") -> output

output = output[order(output$PATNO, output$BM, output$EVENT_ID),]

write.csv(output, "data.csv", row.names=FALSE, na=".")


write.csv(df, "df.csv", row.names=FALSE, na=".")
write.csv(aa[,c(1:4)], "df2.csv", row.names=FALSE, na=".")
