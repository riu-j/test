library("tidyverse")
library("plotROC")
library("pROC")

data <- read.csv("../data/keep/data2.csv")
threshold <- coords(roc(data$diagnosis, data$dv), "best", transpose=TRUE)['threshold']

## draw plots
rocplot <- ggplot(data, aes(d = diagnosis, m = dv)) +
  geom_roc(n.cuts = 0, lwd=1.3) +
  style_roc(xlab = "1-Specificity", ylab = "Sensitivity") +
  ggtitle("ROC plots") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot)["AUC"], 3))) +
  annotate("text", x = .75, y = .15, label = paste("threshold =", round(threshold, 2)))

rocplot

ggsave("../data/DaT_ROC.pdf", rocplot)
