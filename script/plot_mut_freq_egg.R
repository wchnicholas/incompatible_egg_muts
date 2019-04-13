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
library(data.table)
require(cowplot)

mutconverter <- function(m){
  if (m=='V186'){return ('G186V')}
  else if (m=='Y219'){return ('S219Y')}
  else if (m=='F219'){return ('S219F')}
  else if (m=='P194'){return ('L194P')}
  else if (m=='Q156'){return ('H156Q')}
  else if (m=='H246'){return ('N246H')}
  else if (m=='S246'){return ('N246S')}
  else if (m=='L183'){return ('H183L')}
  else if (m=='R156'){return ('H156R')}
  else {return ('NA')}
  }

plot_eggmutfreq <- function(EggFreqTable){
  palette      <- c(brewer.pal(9,"Set1"))
  textsize     <- 7
  p <-  ggplot(EggFreqTable,aes(x=year,y=freq,color=mut)) + 
	  geom_line() +
	  scale_color_manual(values=palette,drop=FALSE) +
	  theme(plot.title=element_text(size=textsize,face="bold"),
                legend.title=element_blank(),
		legend.key.size=unit(3,"mm"),
		legend.text=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=textsize,face="bold")) +
	  ylab(bquote(bold('Frequency'))) +
	  xlab(bquote(bold('Year'))) +
          ggtitle('Major H3N2 egg-adaptive mutations') +
	  scale_y_continuous(lim=c(0,0.75),breaks=c(0,0.25,0.5,0.75),labels=c('0%','25%','50%','75%')) +
	  scale_x_continuous(breaks=c(2003,2005,2007,2009,2011,2013,2015,2017),
			     labels=c('2003','2005','2007','2009','2011','2013','2015','2017'))
  return (p)
  }

plot_seqcount <- function(CountSeqTable){
  textsize <- 7
  p <-  ggplot(filter(CountSeqTable,year>=2003),aes(x=year,y=value,color=variable,group=variable,fill=variable)) + 
	  geom_line() + 
	  geom_point(size=0.7) +
	  geom_area(alpha=0.5,position = "identity") +
	  theme(legend.title=element_text(size=textsize,face="bold"),
		legend.key.size=unit(3,"mm"),
		legend.text=element_text(size=textsize,face="bold"),
		axis.text=element_text(size=textsize,face="bold"),
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=textsize,face="bold")) +
	  ylab(bquote(bold('Total count'))) +
	  xlab(bquote(bold('Year'))) +
          ggtitle(bquote(bold('Major H3N2 egg-adaptive mutations')))+
	  scale_y_continuous(expand=c(0,0),lim=c(0,300)) +
	  scale_x_continuous(breaks=c(2003,2005,2007,2009,2011,2013,2015,2017),
			     labels=c('2003','2005','2007','2009','2011','2013','2015','2017'))
  return (p)
  }

MutLevels    <- c('Q156','R156','L183','V186','P194','Y219','F219','H246','S246')
MutLevels    <- mapply(mutconverter, MutLevels)
EggFreqTable <- read_tsv('table/mut_freq_egg_filtered.tsv') %>% 
                  mutate(mut=mapply(mutconverter,mut)) %>%
                  mutate(mut=factor(mut,levels=MutLevels))
CountSeqTable <- read_tsv('table/seq_count.tsv') %>%
                   mutate(year=as.character(year))
setDT(CountSeqTable)
CountSeqTable <- melt(CountSeqTable,id='year') %>%
                   mutate(year=as.numeric(year))
p1 <- plot_eggmutfreq(EggFreqTable)
ggsave('graph/MutFreq_Egg.png',p1,height=1.5,width=2.8,dpi=600)
