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

plot_score_heatmap <- function(fitness_table, start_resi, end_resi, legend_title, title){
  textsize <- 4.5
  fitness_table <- fitness_table %>%
    filter(pos >= start_resi & pos <= end_resi)
  p <-  ggplot() +
    geom_tile(data=fitness_table,aes(x=pos,y=codon,fill=log10(parameter)),color='black') +
    scale_fill_gradientn(colours=c("white","white","orange"),
                         limits=c(-6.2,-0.36),
                         values=rescale(c(-6.2, -5, 0.36)),
                         breaks=c(-6,-5,-4,-3,-2,-1),
                         labels=c(expression(bold('10'^'-6')), expression(bold('10'^'-5')), expression(bold('10'^'-4')),
                                  expression(bold('10'^'-3')), expression(bold('10'^'-2')), expression(bold('10'^'-1'))),
                         guide="colorbar",
                         na.value="grey50") +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          plot.title=element_text(size=textsize+2,face="bold",hjust=0.5),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title=element_text(size=textsize+1,face="bold"),
          axis.line = element_line(colour = 'black', linewidth = 0),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    guides(fill = guide_colorbar(title.theme=element_text(size=textsize+1,face="bold",colour='black',hjust=0.5),
                                 label.theme=element_text(size=textsize,face="bold",colour='black'),
                                 frame.colour="black",
                                 frame.linewidth = 0.5,
                                 ticks = TRUE,
                                 ticks.colour = "black",
                                 barwidth = 0.5, barheight = 6, title=legend_title)) +
    ggtitle(title) +
    xlab("") +
    ylab("codon")
}

wrapper <- function(df_plot, graphname, legend_title, title){
  df_plot <- df_plot %>%
    select(pos, codon, parameter) #Variable name for input freq may need to chance
  print (range(df_plot$parameter,na.rm=T))
  p1 <- plot_score_heatmap(df_plot, 808, 855, legend_title, title)
  p <- grid.arrange(p1, nrow=1)
  ggsave(graphname,p,width=3.3, height=2.3, dpi=300) #Output file may need to adjust
  }

start_pos <- 808
end_pos   <- 855

df <- read_tsv('result/FP_DMS_codon_freq.tsv') %>%
        filter(!grepl('WT',pos)) %>%
        mutate(codon=factor(codon,levels=sort(unique(codon))))

legend_title <- "freq"

df_plot <- df %>% mutate(parameter=`ipt_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_ipt_heatmap.png', legend_title, 'BAC mutant library (input)')

df_plot <- df %>% mutate(parameter=`Calu3_noAb_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_Calu3_noAb_heatmap.png', legend_title, 'Calu-3 (no Ab)')

df_plot <- df %>% mutate(parameter=`E6_noAb_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_E6_noAb_heatmap.png', legend_title, 'Vero (no Ab)')

df_plot <- df %>% mutate(parameter=`Calu3_CoV44-62_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_Calu3_CoV44-62_heatmap.png', legend_title, 'Calu-3 (COV44-62)')

df_plot <- df %>% mutate(parameter=`E6_CoV44-62_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_E6_CoV44-62_heatmap.png', legend_title, 'Vero (COV44-62)')

df_plot <- df %>% mutate(parameter=`Calu3_CoV44-79_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_Calu3_CoV44-79_heatmap.png', legend_title, 'Calu-3 (COV44-79)')

df_plot <- df %>% mutate(parameter=`E6_CoV44-79_freq`) 
wrapper(df_plot, 'graph/FP_codon_freq_E6_CoV44-79_heatmap.png', legend_title, 'Vero (COV44-79)')
