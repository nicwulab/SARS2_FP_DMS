#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(stringr)
require(cowplot)

plot_score_heatmap <- function(fitness_table, WTresibox, start_resi, end_resi, legend_title){
  textsize <- 5.5
  fitness_table <- fitness_table %>%
    filter(Pos >= start_resi & Pos <= end_resi)
  WTresibox     <- WTresibox %>%
    filter(Pos >= start_resi & Pos <= end_resi) %>%
    mutate(x=x-min(x)+1)
  p <-  ggplot() +
    geom_tile(data=fitness_table,aes(x=resi,y=aa,fill=parameter),color='black') +
    scale_fill_gradientn(colours=c("green","green","white","white","purple","purple"),
                         limits=c(-2.6,1.6),
                         values=rescale(c(-2.6, -1, -0.2, 0.2, 1, 1.5)),
                         breaks=c(-2,-1,0,1),
                         labels=c('-2','-1','0','1'),
                         guide="colorbar",
                         na.value="grey50") +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title=element_text(size=7,face="bold"),
          axis.line = element_line(colour = 'black', linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                 label.theme=element_text(size=7,face="bold",colour='black'),
                                 frame.colour="black",
                                 frame.linewidth = 0.5,
                                 ticks = TRUE,
                                 ticks.colour = "black",
                                 barwidth = 0.5, barheight = 6, title=legend_title)) +
    geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.2) +
    xlab("") +
    ylab("amino acid")
}

wrapper <- function(df_plot, graphname, legend_title){
  df_plot <- df_plot %>%
    mutate(parameter=case_when(str_sub(resi,1,1)==aa ~ 0, TRUE ~ parameter)) %>% #Set WT parameter (usually 0 or 1)
    mutate(Mutation=paste(resi,aa,sep='')) %>%
    select(Mutation, avg_ipt_freq, resi, Pos, aa, parameter) #Variable name for input freq may need to chance
  print (range(df_plot$parameter,na.rm=T))
  p1 <- plot_score_heatmap(df_plot, WTresibox, 808, 855, legend_title)
  p <- grid.arrange(p1, nrow=1)
  ggsave(graphname,p,width=4, height=2, dpi=300) #Output file may need to adjust
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','*'))

start_pos <- 808
end_pos   <- 855

df <- read_tsv('result/FP_DMS_fit.tsv') %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  mutate(resi=str_sub(mut,1,-2))
residues <- unique(df$resi)

df <- df %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  filter(avg_ipt_freq>= 0.0001) %>%
  mutate(resi=str_sub(mut,1,-2)) %>%
  mutate(aa=str_sub(mut,-1,-1)) %>%
  mutate(aa=ifelse(aa=='_', '*', aa)) %>%
  filter(aa %in% aa_level) %>%
  mutate(aa=factor(aa,levels=aa_level)) %>%
  mutate(resi=factor(resi,levels=unique(residues))) %>%
  complete(resi, aa) %>%
  mutate(Pos=str_sub(resi,2,-1)) %>%
  mutate(Pos=as.numeric(as.character(Pos))) %>%
  arrange(Pos)

WTresibox  <- df %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(1,end_pos-start_pos+1)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,Pos,x, y)

legend_title <- "escape"
fit_cutoff <- 0.75

df_plot <- df %>%
             mutate(parameter=`fit_Calu3_CoV44-62` - `fit_Calu3_noAb`) %>%
             mutate(parameter=ifelse(`fit_Calu3_noAb`>fit_cutoff, parameter, NA)) 
wrapper(df_plot, 'graph/FP_escape_Calu3_CoV44-62.png', legend_title)

df_plot <- df %>%
             mutate(parameter=`fit_Calu3_CoV44-79` - `fit_Calu3_noAb`) %>%
             mutate(parameter=ifelse(`fit_Calu3_noAb`>fit_cutoff, parameter, NA)) 
wrapper(df_plot, 'graph/FP_escape_Calu3_CoV44-79.png', legend_title)

df_plot <- df %>%
             mutate(parameter=`fit_E6_CoV44-62` - `fit_E6_noAb`) %>%
             mutate(parameter=ifelse(`fit_E6_noAb`>fit_cutoff, parameter, NA)) 
wrapper(df_plot, 'graph/FP_escape_E6_CoV44-62.png', legend_title)

df_plot <- df %>%
             mutate(parameter=`fit_E6_CoV44-79` - `fit_E6_noAb`) %>%
             mutate(parameter=ifelse(`fit_E6_noAb`>fit_cutoff, parameter, NA)) 
wrapper(df_plot, 'graph/FP_escape_E6_CoV44-79.png', legend_title)
