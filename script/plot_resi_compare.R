#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(viridis)
library(qualpalr)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_E6_vs_Calu3 <- function(df, graphname){
  print (paste('correlation for:', graphname, cor(df$`mean_fit_P1-Calu3`, df$`mean_fit_P1-E6`)))
  textsize <- 7
  p <- ggplot(df,aes(x=`mean_fit_P1-Calu3`, y=`mean_fit_P1-E6`)) +
    geom_point(size=0.6, alpha=0.5, color='grey30', pch=16) +
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold"),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    labs(x=bquote(bold('mutation tolerance (Calu3)')),y=bquote(bold('mutation tolerance (Vero E6)'))) +
    xlim(-0.5,1.7) +
    ylim(-0.5,1.7)
  ggsave(graphname, p, height=2, width=2, dpi=600)
  }

df <- read_tsv('result/FP_DMS_fit_by_resi.tsv')
plot_E6_vs_Calu3(df, 'graph/compare_E6_vs_Calu3.png')
