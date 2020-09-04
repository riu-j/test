#set parameters
{
  proj = "s128"
  folder_name = "G2019S_dydx/"
  dv_name = c("DaT scan (SBR)", "MDS-UPDRS Part1", "MDS-UPDRS Part2", "MDS-UPDRS Part3", "SCOPA-AUT", "α-synuclein", "Aβ1-42")

  folder_path = paste("../notes/", proj, "/", folder_name, proj,  "_", sep="")
  table = read.csv(paste(folder_path, "table.csv", sep=""))
  offsetT = read.csv(paste(folder_path, "offsetT.csv", sep=""), stringsAsFactors=FALSE)
  
  df = merge(table, offsetT, by.x="ID", by="id", all=T)
  df$TIME = df$TIME + df$offsetT
  df$onset = df$age-df$offsetT
  
  for(i in 1:length(dv_name)){
    df$BM = replace(df$BM, which(df$BM==i), dv_name[i])
  }
  df$Sex = replace(df$Sex, which(df$Sex==1), "female")
  df$Sex = replace(df$Sex, which(df$Sex==0), "male")
  df$G2019s = replace(df$G2019s, which(df$G2019s==1), "+")
  df$G2019s = replace(df$G2019s, which(df$G2019s==0), "-")
}


#define function
{
  MainBiomarker <- function(x, bm){
    alpha[bm] + beta[bm]/gamma[bm]*(exp(gamma[bm]*x*dt)-1)
  }
  
  OtherBiomarker <- function(x, bm){
    alpha[bm] + beta[bm]/gamma[bm]*(exp(gamma[bm]*x*dt)-1)
  }
}

#libraries
{
  library("ggplot2")
  library("dplyr")
}

#DV vs CIPRED
{svg(paste(folder_path, "vpc_cipred2.svg", sep=""), width=16, height = 12, bg="transparent")
(df %>% 
    ggplot(aes(x = DV, y = CIPRED, colour = interaction(factor(Sex), factor(G2019s), sep=""))) +
    geom_point(shape=16, size=2, alpha=0.5) + 
    geom_smooth(method = "lm", linetype=2) + 
    geom_abline(slope=1, intercept=0, linetype=2) +
    theme_bw(base_size = 11, base_family = "") +
    scale_color_hue(name="") +
    facet_wrap(~ BM, scales = "free") -> fig_cipred)%>%
  print()
dev.off()
# plotly::ggplotly()
# ggsave("cipred.pdf", fig_cipred)
}

#TIME vs CWRES
{svg(paste(folder_path, "vpc_cwres2.svg", sep=""), width=16, height = 12, bg="transparent")
(df %>% 
    ggplot(aes(x = TIME, y = CWRES, colour = interaction(factor(Sex), factor(G2019s), sep=""))) +
    geom_point(shape=16, size=2, alpha=0.5) + 
    geom_smooth(method = "lm", linetype=2) + 
    geom_abline(slope=0, intercept=0, linetype=2) +
    theme_bw(base_size = 11, base_family = "") +
    scale_color_hue(name="") +
    facet_wrap(~ factor(BM), scales = "free") -> fig_cwres) %>%
  print()
dev.off()
}

#offsetT vs onset
{svg(paste(folder_path, "vpc_onset2.svg", sep=""), width=5, height=4, bg="transparent")
(df[!duplicated(df$ID),] %>% 
    ggplot(aes(x = offsetT, y = onset, colour = interaction(factor(Sex), factor(G2019s), sep=""))) +
    geom_point(shape=16, size=1, alpha=0.5) + 
    geom_smooth(method = "lm", linetype=2) + 
    scale_color_hue(name="") +
    theme_bw(base_size = 11, base_family = "") -> fig_onset) %>%
  print()
dev.off()
# plotly::ggplotly()
# ggsave("onset.pdf", fig_onset)
}

# df[!duplicated(df$ID),] %>% 
#     ggplot(aes(x = offsetT, y = onset, colour = interaction(factor(Sex), factor(G2019s), sep=""))) +
#     geom_point(shape=16, size=1, alpha=0.5) + 
#     geom_smooth(method = "lm", linetype=2) + 
#     scale_color_hue(name="") +
#     theme_bw(base_size = 11, base_family = "") %>%
#   print()

#plot
# df_a = df[is.na(df$BM1)==F,]
(df %>% 
    ggplot(aes(x = TIME, y = DV, colour = interaction(factor(Sex), factor(G2019S), sep=""), group=id)) +
    geom_line(size=0.3, alpha=0.5) + 
    geom_point(size=0.3, alpha=0.5) +
    # stat_function(fun=graph_main_non) +
    scale_color_hue(name="") +
    theme_bw(base_size = 11, base_family = "") +
    facet_wrap(~ factor(BM), scales = "free")
    -> fig_cwres) %>%
  plotly::ggplotly()
ggsave("cwres.pdf", fig_cwres)

