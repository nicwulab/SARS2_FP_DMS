#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(qualpalr)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(data.table)
library(gridExtra)
library(stringr)
require(cowplot)

plot_mean_fit <- function(df, graphname){
  textsize <- 7
  palette <- brewer.pal(3,"Accent")
  p <- ggplot(df,aes(x=pos, y=value, color=variable)) +
    geom_point(size=1, alpha=1, pch=16) +
    geom_line() +
    scale_color_manual(values=palette,drop=FALSE,
                       labels = c("Calu-3", "Vero E6")) +
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.05, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize,face="bold"),
          legend.justification = "center",
          legend.position='top') +
    labs(x=bquote(bold('residue position')),y=bquote(bold('mutation tolerance')))
  ggsave(graphname, p, height=2, width=3, dpi=600)
  }

df <- read_tsv('result/FP_DMS_fit_by_resi.tsv') %>%
        select(resi, pos, mean_fit_Calu3_noAb, mean_fit_E6_noAb) %>%
        data.table(.) %>%
        melt(., id=c('resi', 'pos'))
print (df)
plot_mean_fit(df, 'graph/mean_fit.png')


