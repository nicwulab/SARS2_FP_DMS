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

plot_mut_count <- function(df, graphname, h, w){
  textsize <- 8
  palette <- brewer.pal(3,"Set1")
  p <-  ggplot(df,aes(x=num_mut, y=read_count, fill=rep)) +
	  geom_bar(stat="identity",position=position_dodge(),width=0.8) +
          scale_fill_manual(values=palette,drop=FALSE) +
	  theme_cowplot(12) +
	  theme(plot.title=element_text(size=textsize,face='bold',hjust = 0.5),,
		plot.background = element_rect(fill = "white"),
		axis.title=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		legend.key.size=unit(0.1,'in'),
		legend.spacing.x=unit(0.05, 'in'),
		legend.title=element_blank(),
		legend.text=element_text(size=textsize,face="bold"),
		legend.justification = "center",
		legend.position='right') +
          ggtitle('BAC mutant library') +
	  labs(x='# of amino acid mutations',y=bquote(bold('% of reads'))) +
	  scale_y_continuous(limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,20,40,60,80,100)) +
	  scale_x_continuous(limits=c(-0.5,3.5), breaks=seq(0,3),labels=c('0','1','2','≥3'))

  ggsave(graphname, p, height=h, width=w, dpi=600)
  }

df_rep1 <- read_tsv('result/Lib1_mut_count.tsv') %>%
             mutate(read_count=read_count/sum(read_count)) %>%
             mutate(rep='Replicate 1')
df_rep2 <- read_tsv('result/Lib2_mut_count.tsv') %>%
             mutate(read_count=read_count/sum(read_count)) %>%
             mutate(rep='Replicate 2')
df <- rbind(df_rep1, df_rep2)
print (filter(df, num_mut<4))
plot_mut_count(df, 'graph/Mut_count_input_lib.png', 1.5, 3)
