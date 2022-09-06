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

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plot_hist <- function(t, main_title){
  textsize <- 7
  p <- ggplot(t,aes(x=fit)) +
    geom_histogram(binwidth=0.1) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize,face="bold", hjust=0.5),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    ggtitle(main_title) +
    scale_fill_manual(values=c('black'),drop=FALSE) +
    labs(y=expression(bold('count')),x=expression(bold('fitness'))) +
    coord_cartesian(xlim=c(-2.5,1))
  return (p)
}

plot_replicate_cor <- function(df, graphname, param){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  textsize <- 7
  df$density <- get_density(df$rep1, df$rep2, n = 100)
  p <- ggplot(df,aes(x=rep1, y=rep2))+#, color=density)) +
    #geom_hex(bins = 70) +
    #scale_fill_continuous(type = "viridis") +
    geom_point(size = 0.3, alpha=0.5, color='grey30', pch=16)+
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
    labs(x=bquote(bold(paste(.(param),' (replicate 1)'))),y=bquote(bold(paste(.(param),' (replicate 2)'))))
  ggsave(graphname, p, height=2, width=2.5, dpi=3000)
}

mut_classification <- function(mut_class, resi){
  if (mut_class=='missense' & resi=='W64'){return ('W64X')}
  else(return (mut_class))
}

plot_by_class <- function(df, graphname, ylab){
  df <- df %>%
    filter(mut_class != 'WT')
  textsize <- 7
  p <- ggplot(df,aes(x=mut_class, y=score, group=mut_class)) +
    geom_violin(width=1, color="black") +
    geom_sina(pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='black', alpha=0.2) +
    geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0) + 
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold"),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    ylab(ylab)
  ggsave(graphname, p, height=2, width=2,dpi=3000)
}

t_test <- function(df, class_1, class_2){
  p_value <- t.test(filter(df_exp, mut_class==class_1)$score, filter(df_exp, mut_class==class_2)$score)$p.value
  print (paste("p-value of diff between", class_1, 'vs', class_2, ':', p_value))
}

df <- read_tsv('result/FP_DMS_fit.tsv') %>%
  filter(avg_ipt_freq >= 0) %>%
  filter(mut != "WT")
print (nrow(df))
df_exp <- df %>%
  rename(rep1=fit_Lib1) %>%
  rename(rep2=fit_Lib2) %>%
  rename(score=fit)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit.png', "fitness")
plot_by_class(df_exp, 'graph/QC_fit_by_class.png', 'fitness')
t_test(df, 'silent', 'nonsense')
t_test(df, 'silent', 'missense')
t_test(df, 'missense', 'nonsense')
#write.table(select(filter(df_exp, mut_class=='missense'), mut, Input_freq, score), 'result/FP_DMS_expression_score.tsv', quote=FALSE, sep="\t", row.names=FALSE)
