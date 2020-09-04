for(i in 1:6){
table2 <- table[table$BM==i,]
rows <- nrow(table2)
nrows2 <- nrow(table2[!duplicated(table2$ID),])
print(paste("number of observations BM", i, ": ", rows, sep=""))
print(paste("number of subjects BM", i, ": ", nrows2, sep=""))
}

attach(table)                       # 使うデータを固定
fc <- factor(ID)                    # 要素をグループ化
levels(fc)                               # グループ化されているかを確認

test <- tapply(TIME, ID, max)                 # breaks の平均をグループごとに計算
