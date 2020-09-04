#import
library("genetics")

#input
data <- read.csv("../data/processed/200903_235046_gene_target.csv")
process_time <- "200903_235046"

#analysis
gene_name <- colnames(data)[-1]
result <- data.frame(gene=gene_name, major=rep(NA, length(gene_name)), hetero=rep(NA, length(gene_name)), 
                    minor=rep(NA, length(gene_name)), maf=rep(NA, length(gene_name)), 
                    x2_p=rep(NA, length(gene_name)), exact_p=rep(NA, length(gene_name)))
for (i in 1:length(gene_name)){
  gene <- gene_name[i]
  tmp_gene <- data[gene]
  major <- max(sum(tmp_gene == 0), sum(tmp_gene == 2))
  hetero <- sum(tmp_gene == 1)
  minor <- min(sum(tmp_gene == 0), sum(tmp_gene == 2))
  maf <- min(c((major*2+hetero)/(major+hetero+minor)/2, (minor*2+hetero)/(major+hetero+minor)/2))
  result[i, 2] <- major
  result[i, 3] <- hetero
  result[i, 4] <- minor
  result[i, 5] <- maf %>% signif(digits=4)
  if (maf < 0.05){
    result[i, 6] <- "."
    result[i, 7] <- "."
  }else{
    three.data <- c(rep("A/A", major),
                   rep("A/C", hetero),
                   rep("C/C", minor))
    g3 <- genotype(three.data)
    hwe_chi <- HWE.chisq(g3)
    hwe_exa <- HWE.exact(g3)
    result[i, 6] <- hwe_chi$p.value %>% signif(digits=4)
    result[i, 7] <- hwe_exa$p.value %>% signif(digits=4)
  }
}

#output
write.csv(result, paste("../data/processed/", process_time, "_hwe_result.csv", sep=""), quote=FALSE, row.names=FALSE)
