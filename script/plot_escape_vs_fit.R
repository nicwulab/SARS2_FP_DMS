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

plot_escape_vs_fit <- function(df, graphname, cell_line, Ab){
  textsize <- 7
  p <- ggplot(df,aes(x=fit, y=escape)) +
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
  labs(x=bquote(bold(paste('fitness (',.(cell_line),')',sep=''))),y=bquote(bold(paste(.(Ab),' escape (',.(cell_line),')',sep=''))))
  ggsave(graphname, p, height=2, width=2.5, dpi=3000)
  }

df <- read_tsv('result/FP_DMS_fit.tsv') %>%
  filter(avg_ipt_freq >= 0.0001) %>%
  filter(mut != "WT")

df_esp <- df %>%
  rename(fit=`fit_P1-Calu3_noAb`) %>%
  mutate(escape=`fit_P1-Calu3_CoV44-62` - fit)
plot_escape_vs_fit(df_esp, 'graph/fit_vs_escape_Calu3_CoV44-62.png', 'Calu3', 'CoV44-62')

df_esp <- df %>%
  rename(fit=`fit_P1-Calu3_noAb`) %>%
  mutate(escape=`fit_P1-Calu3_CoV44-79` - fit)
plot_escape_vs_fit(df_esp, 'graph/fit_vs_escape_Calu3_CoV44-79.png', 'Calu3', 'CoV44-79')

df_esp <- df %>%
  rename(fit=`fit_P1-E6_noAb`) %>%
  mutate(escape=`fit_P1-E6_CoV44-62` - fit)
plot_escape_vs_fit(df_esp, 'graph/fit_vs_escape_E6_CoV44-62.png', 'Vero E6', 'CoV44-62')

df_esp <- df %>%
  rename(fit=`fit_P1-E6_noAb`) %>%
  mutate(escape=`fit_P1-E6_CoV44-79` - fit)
plot_escape_vs_fit(df_esp, 'graph/fit_vs_escape_E6_CoV44-79.png', 'Vero E6', 'CoV44-79')
