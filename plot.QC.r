data<-read.table('./QC.stats')
colnames(data)=c('Category','Number','Sample')
data$Sample=factor(data$Sample,levels=unique(data[,3]))
data$Category=factor(data$Category,levels=unique(data[,1]))
library(ggplot2)
library(cowplot)
library(extrafont)
library(RColorBrewer)

pdf('QC_HiC_pairs.pdf',width=5,height=4,family='Arial')
ggplot(data, aes (x=Category, y = Number,fill=Sample)) + geom_bar(stat='identity',position=position_dodge()) + ylab("Reads") +theme_cowplot()+theme(axis.text.x = element_text(angle =90, hjust = 1))+xlab('')+scale_y_continuous(limits=c(-1,1200000000),expand=c(0,0))+theme(axis.text=element_text(size=9),axis.title=element_text(size=9,face="plain"))+theme(legend.title=element_text(size=9),legend.key.size=unit(0.2,'cm'),legend.text=element_text(size=8))

#+ scale_fill_manual(values=rev(brewer.pal(10,'Set3'))[c(2,1,3)])
dev.off()
