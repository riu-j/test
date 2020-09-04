#set parameters
{
  proj <- "s129"
  subfolder <- ""
  dv_name <- c("DaT scan (SBR)", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "MDS-UPDRS Part3")#, "SCOPA-AUT", "α-synuclein", "Aβ1-42")

  folder_path <- paste("../notes/", proj, "/", subfolder , proj, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT <- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  
  df <- merge(table, offsetT, by.x="ID", by="id", all=T)
  rm(table, offsetT)
  df$TIME <- df$TIME + df$offsetT
  df$onset <- df$age-df$offsetT
  
  for(i in 1:length(dv_name)){
    df$BM <- replace(df$BM, which(df$BM==i), dv_name[i])
  }
  df$Sex <- replace(df$Sex, which(df$Sex==1), "female")
  df$Sex <- replace(df$Sex, which(df$Sex==0), "male")
  df$G2019s <- replace(df$G2019s, which(df$G2019s==1), "+")
  df$G2019s <- replace(df$G2019s, which(df$G2019s==0), "-")
}

#libraries
library("tidyverse")

#CIPRED vs DV
cipred <- ggplot(df, aes(x=DV, y=CIPRED)) +
  geom_point(alpha=0.3) + 
  stat_smooth(method="lm") + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  facet_wrap(~ factor(BM), scales="free") +
  ggtitle("CIPRED vs TIME")
ggsave(paste(folder_path, "vpc_cipred.pdf", sep=""), cipred)

#CWRES vs TIME
cwres <- ggplot(df, aes(x=TIME, y=CWRES)) +
  geom_point(alpha=0.3) + 
  stat_smooth(method="lm") + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  facet_wrap(~ factor(BM), scales="free") +
  ggtitle("CWRES vs TIME")
ggsave(paste(folder_path, "vpc_cwres.pdf", sep=""), cwres)

#onset vs offsetT
onset <- ggplot(df[!duplicated(df$ID),], aes(x=offsetT, y=onset)) +
  geom_point(shape=16, size=1.5) +
  stat_smooth(method="lm") +
  ggtitle("onset vs offsetT")
ggsave(paste(folder_path, "vpc_onset.pdf", sep=""), onset, width=4, height=4)

# {svg(paste(folder_path, "vpc_cipred.svg", sep=""), width=16, height=12, bg="transparent")
#   df %>% 
#     ggplot(aes(x=DV, y=CIPRED, colour=interaction(factor(Sex), factor(G2019s), sep=""))) +
#     geom_point(shape=16, size=2, alpha=0.5) + 
#     geom_smooth(method="lm", linetype=2) + 
#     geom_abline(slope=1, intercept=0, linetype=2) +
#     theme_bw(base_size=11, base_family="") +
#     scale_color_hue(name="") +
#     facet_wrap(~ BM, scales="free") %>%
#     print()
#   dev.off()
# }
#offsetT vs onset
# {svg(paste(folder_path, "vpc_onset.svg", sep=""), width=5, height=4, bg="transparent")
# df[!duplicated(df$ID),] %>%
#     ggplot(aes(x=offsetT, y=onset, colour=interaction(factor(Sex), factor(G2019s), sep=""))) +
#     geom_point(shape=16, size=1, alpha=0.5) + 
#     geom_smooth(method="lm", linetype=2) + 
#     scale_color_hue(name="") +
#     theme_bw(base_size=11, base_family="") %>%
#   print()
# dev.off()
# }

ggplot(df, aes(x=TIME, y=DV, colour=interaction(factor(Sex), factor(G2019s), sep=""), group=ID))+
  geom_line(alpha=0.3) +
  geom_point(size=0.8, alpha=0.3) +
  scale_color_hue(name="group") +
  facet_wrap(~ factor(BM), scales="free")
