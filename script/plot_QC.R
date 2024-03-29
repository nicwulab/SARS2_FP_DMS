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

plot_replicate_cor <- function(df, graphname, param, title){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  textsize <- 7
  df$density <- get_density(df$rep1, df$rep2, n = 100)
  p <- ggplot(df,aes(x=rep1, y=rep2)) +
    geom_point(size=0.6, alpha=0.5, color='grey30', pch=16) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize,face="bold",family="Arial",hjust=0.5),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold",family="Arial"),
          axis.text=element_text(size=textsize,face="bold",family="Arial"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold",family="Arial"),
          legend.text=element_text(size=textsize,face="bold",family="Arial"),
          legend.position='right') +
    ggtitle(title) +
    labs(x=bquote(bold(paste(.(param),' (replicate 1)'))),y=bquote(bold(paste(.(param),' (replicate 2)'))))
  ggsave(graphname, p, height=1.5, width=1.5, dpi=600)
}

mut_classification <- function(mut_class, resi){
  if (mut_class=='missense' & resi=='W64'){return ('W64X')}
  else(return (mut_class))
}

plot_by_class <- function(df, graphname, ylab, title){
  t_test(df, 'silent', 'nonsense')
  t_test(df, 'silent', 'missense')
  t_test(df, 'missense', 'nonsense')

  df <- df %>%
    filter(mut_class != 'WT') %>%
    mutate(mut_class = str_replace(mut_class, "silent", "sil")) %>%
    mutate(mut_class = str_replace(mut_class, "missense", "mis")) %>%
    mutate(mut_class = str_replace(mut_class, "nonsense", "non"))
  textsize <- 7
  p <- ggplot(df,aes(x=mut_class, y=score, group=mut_class)) +
    geom_violin(width=1, color="black") +
    geom_sina(pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='black', alpha=0.2) +
    geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0) + 
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize,face="bold",family="Arial", hjust = 0.5),
          plot.background = element_rect(fill = "white"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=textsize,face="bold",family="Arial"),
          axis.text.x=element_text(size=textsize,face="bold",family="Arial",angle=90,hjust=1,vjust=0.5),
          axis.text.y=element_text(size=textsize,face="bold",family="Arial"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold",family="Arial"),
          legend.text=element_text(size=textsize,face="bold",family="Arial"),
          legend.position='right') +
    ggtitle(title) +
    ylab(ylab) +
    xlab("") #+
    #ylim(0,4.6)
  ggsave(graphname, p, height=1.5, width=1.5,dpi=600)
}

t_test <- function(df_exp, class_1, class_2){
  p_value <- t.test(filter(df_exp, mut_class==class_1)$score, filter(df_exp, mut_class==class_2)$score)$p.value
  print (paste("p-value of diff between", class_1, 'vs', class_2, ':', p_value))
}

df <- read_tsv('result/FP_DMS_fit.tsv') %>%
  filter(avg_ipt_freq >= 0.0001) %>%
  filter(mut != "WT")
print (nrow(df))

df_exp <- df %>%
  rename(rep1=`fit_P0_Rep1`) %>%
  rename(rep2=`fit_P0_Rep2`) %>%
  rename(score=`fit_P0`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_P0.png', "fitness", "P0")
plot_by_class(df_exp, 'graph/QC_fit_by_class_P0.png', 'fitness', "P0")

df_exp <- df %>%
  rename(rep1=`fit_Calu3_noAb_Rep1`) %>%
  rename(rep2=`fit_Calu3_noAb_Rep2`) %>%
  rename(score=`fit_Calu3_noAb`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_Calu3_noAb.png', "fitness", "Calu-3 (no Ab)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_Calu3_noAb.png', 'fitness', "Calu-3 (no Ab)")

df_exp <- df %>%
  rename(rep1=`fit_E6_noAb_Rep1`) %>%
  rename(rep2=`fit_E6_noAb_Rep2`) %>%
  rename(score=`fit_E6_noAb`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_E6_noAb.png', "fitness", "Vero (no Ab)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_E6_noAb.png', 'fitness', "Vero (no Ab)")

df_exp <- df %>%
  rename(rep1=`fit_Calu3_CoV44-62_Rep1`) %>%
  rename(rep2=`fit_Calu3_CoV44-62_Rep2`) %>%
  rename(score=`fit_Calu3_CoV44-62`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_Calu3_CoV44-62.png', "fitness", "Calu-3 (COV44-62)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_Calu3_CoV44-62.png', 'fitness', "Calu-3 (COV44-62)")

df_exp <- df %>%
  rename(rep1=`fit_E6_CoV44-62_Rep1`) %>%
  rename(rep2=`fit_E6_CoV44-62_Rep2`) %>%
  rename(score=`fit_E6_CoV44-62`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_E6_CoV44-62.png', "fitness", "Vero (COV44-62)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_E6_CoV44-62.png', 'fitness', "Vero (COV44-62)")

df_exp <- df %>%
  rename(rep1=`fit_Calu3_CoV44-79_Rep1`) %>%
  rename(rep2=`fit_Calu3_CoV44-79_Rep2`) %>%
  rename(score=`fit_Calu3_CoV44-79`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_Calu3_CoV44-79.png', "fitness", "Calu-3 (COV44-79)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_Calu3_CoV44-79.png', 'fitness', "Calu-3 (COV44-79)")

df_exp <- df %>%
  rename(rep1=`fit_E6_CoV44-79_Rep1`) %>%
  rename(rep2=`fit_E6_CoV44-79_Rep2`) %>%
  rename(score=`fit_E6_CoV44-79`)
plot_replicate_cor(df_exp, 'graph/QC_replicate_fit_E6_CoV44-79.png', "fitness", "Vero (COV44-79)")
plot_by_class(df_exp, 'graph/QC_fit_by_class_E6_CoV44-79.png', 'fitness', "Vero (COV44-79)")

