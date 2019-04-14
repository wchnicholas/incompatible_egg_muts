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
require(cowplot)

plot_freq <- function(t, sample_ID, mut_levels, freq_cutoff){
  palette      <- c(brewer.pal(9,"Set1"))
  palette      <- c('grey30',palette[1],palette[3:4],palette[7],palette[2],palette[9],'yellow3',palette[5])
  P_levels <- c('0','1','2','3','4','5')
  sample_t <- t %>%  
		filter(sample==sample_ID) %>%
		select(-sample) %>%
                filter(mut!='WT') %>%
		cast(mut~passage, value.var=mut_freq) %>%
		mutate(max_freq=apply(., 1, max)) %>%
		filter(max_freq > freq_cutoff) %>%
		melt() %>%
                mutate(variable=unlist(lapply(variable, function(x) str_replace(x,'P','')))) %>%
                mutate(variable=unlist(lapply(variable, function(x) str_replace(x,'DNA','0')))) %>%
		filter(variable %in% P_levels) %>%
		mutate(variable=factor(variable, levels=P_levels)) %>%
                mutate(mut=factor(mut, levels=mut_levels))
  textsize <- 7
  p <-  ggplot(sample_t, aes(x=variable,y=value,group=mut,color=mut)) + 
	  geom_line() + 
	  scale_color_manual(values=palette,drop=FALSE) +
	  theme(legend.title=element_blank(),
		legend.key.size=unit(3,"mm"),
                legend.position='none',
                #legend.position='bottom',
		legend.text=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
		axis.title=element_text(size=textsize,face="bold")) +
	  ylab(bquote(bold('Frequency'))) +
	  xlab(bquote(bold('Passage'))) +
          scale_y_continuous(breaks=c(0,0.5,1),labels=c('0.0','0.5','1.0'), lim=c(0,1)) +
          guides(color=guide_legend(nrow=3, override.aes = list(size=0.5))) 
  ggsave(paste('graph/Freq_',sample_ID,'.png',sep=''),p,height=1,width=1.8)
  #ggsave(paste('graph/Freq_',sample_ID,'.png',sep=''),p,height=1,width=3)
  }

t <- read_tsv('result/mut_freq.tsv')
freq_cutoff <- 0.1
sample_IDs <- c('G186V_rep1','G186V_rep2','G186V_rep3','L194P_rep1')
mut_levels <- t %>%
                filter(mut_freq > freq_cutoff) %>%
                filter(sample %in% sample_IDs) %>%
                filter(mut != 'WT') %>%
                .$mut %>%
                unique()
print (mut_levels)
for (sample_ID in sample_IDs){
  plot_freq(t, sample_ID, mut_levels, freq_cutoff)
  }
