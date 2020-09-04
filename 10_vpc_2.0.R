#libraries
library("tidyverse")
library("ggExtra")

#set parameters
{
  projectid <- "s132"
  folder_name <- "base_3045.046/" #ex) "base/"
  dv_name <- c("DaT scan (SBR)", "MDS-UPDRS Part I", "MDS-UPDRS Part II", "MDS-UPDRS Part III", "SCOPA-AUT", "α-synuclein")#, "Aβ1-42")

  folder_path <- paste("../notes/", projectid, "/", folder_name , projectid, "_", sep="")
  table <- read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT <- read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  
  df <- merge(table, offsetT, by="ID", all=T)
  # rm(table, offsetT)
  df$TIME <- df$TIME + df$offsetT
  df$onset <- df$age - df$offsetT
  
  for(i in 1:length(dv_name)){
    df$BM <- replace(df$BM, which(df$BM==i), dv_name[i])
  }
  # df$Sex <- replace(df$Sex, which(df$Sex>0), "female")
  # df$Sex <- replace(df$Sex, which(df$Sex==0), "male")
  # df$G2019s <- replace(df$G2019s, which(df$G2019s==1), "+")
  # df$G2019s <- replace(df$G2019s, which(df$G2019s==0), "-")
}

#CIPRED vs DV
cipred <- ggplot(df, aes(x=DV, y=CIPRED)) +
  geom_point(size=0.5, alpha=0.7) + 
  stat_smooth(method="lm") + 
  geom_abline(slope=1, intercept=0, linetype=2) +
  facet_wrap(~ BM, scales="free", ncol=6) +
  ggtitle("CIPRED vs DV")
# ggsave(paste(folder_path, "vpc_cipred.pdf", sep=""), cipred)
ggsave(paste(folder_path, "vpc_cipred2.png", sep=""), cipred, width=18, height=3.3)

#CWRES vs TIME
cwres <- ggplot(df, aes(x=TIME, y=CWRES)) +
  geom_point(size=0.5, alpha=0.7) + 
  stat_smooth(method="lm") + 
  geom_abline(slope=0, intercept=0, linetype=2) +
  facet_wrap(~ BM, scales="free", ncol=6) +
  ggtitle("CWRES vs TIME")
# ggsave(paste(folder_path, "vpc_cwres.pdf", sep=""), cwres)
ggsave(paste(folder_path, "vpc_cwres2.png", sep=""), cwres, width=18, height=3.3)

#onset vs offsetT
onset <- ggplot(df[!duplicated(df$ID),], aes(x=offsetT, y=onset)) +
  geom_point(shape=16, size=1.5) +
  stat_smooth(method="lm") +
  ggtitle("onset vs offsetT") +
  guides(colour=FALSE)

onset <- ggMarginal(
  onset, 
  type = "density",
  margins = "both",
  size = 5,
  groupColour = FALSE,
  groupFill = FALSE
)
# ggsave(paste(folder_path, "vpc_onset.pdf", sep=""), onset, width=4, height=4)
ggsave(paste(folder_path, "vpc_onset2.png", sep=""), onset, width=4, height=4)

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
