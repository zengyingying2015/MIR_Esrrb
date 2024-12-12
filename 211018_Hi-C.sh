
###E14_2i_HiC-Arima-R1-filtered.mcool    E14_RA_HiC-Arima-R1-filtered.mcool
#E14_EPSC_HiC-Arima-R1-filtered.mcool

#5kb
pyBHFDR -O E14_2i_HiC-Arima-BHFDR-loops.5kb.txt -p E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/5000  --pw 4  --ww 7 --nproc 20
pyHICCUPS -O E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt -p E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/5000 --pw 4 --ww  7 --only-anchors --nproc 20


pyBHFDR -O E14_EPSC_HiC-Arima-BHFDR-loops.5kb.txt -p E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/5000  --pw 4  --ww 7 --nproc 20
pyHICCUPS -O E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt -p E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/5000 --pw  4 --ww  7 --only-anchors --nproc 20

pyBHFDR -O E14_RA_HiC-Arima-BHFDR-loops.5kb.txt -p E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/5000  --pw 4  --ww 7 --nproc 20
pyHICCUPS -O E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt -p E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/5000 --pw  4 --ww  7 --only-anchors --nproc 20


#10kb
pyBHFDR -O E14_2i_HiC-Arima-BHFDR-loops.10kb.txt -p E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000  --pw 2  --ww 5 --nproc 20
pyHICCUPS -O E14_2i_HiC-Arima-HICCUPS-loops.10kb.txt -p E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 --pw 2 --ww  5 --only-anchors --nproc 20
pyBHFDR -O E14_EPSC_HiC-Arima-BHFDR-loops.10kb.txt -p E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000  --pw 2  --ww 5 --nproc 20
pyHICCUPS -O E14_EPSC_HiC-Arima-HICCUPS-loops.10kb.txt -p E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 --pw  2 --ww  5 --only-anchors --nproc 20
pyBHFDR -O E14_RA_HiC-Arima-BHFDR-loops.10kb.txt -p E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000  --pw 2  --ww 5 --nproc 20
pyHICCUPS -O E14_RA_HiC-Arima-HICCUPS-loops.10kb.txt -p E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 --pw  2 --ww  5 --only-anchors --nproc 20



################
#detect single level TAD
domaincaller --uri E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -O E14_2i_HiC-Arima.tad.bed --DI-output E14_2i_HiC-Arima.DIs.bedGraph --removeCache -p 20 --logFile E14_2i_HiC-Arima_domaincaller.log

domaincaller --uri E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -O E14_EPSC-Arima.tad.bed --DI-output E14_EPSC-Arima.DIs.bedGraph --removeCache -p 20 --logFile E14_EPSC-Arima_domaincaller.log

domaincaller --uri E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -O E14_RA_HiC-Arima.tad.bed --DI-output E14_RA_HiC-Arima.DIs.bedGraph --removeCache -p 20 --logFile E14_RA_HiC-Arima_domaincaller.log


##Aggregate Peak Analysis

apa-analysis -O E14_2i_HiC_E14_2i_HiC-loops-Arima-apa.png -p E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/5000 -I E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt -U
apa-analysis -O E14_EPSC_HiC_E14_EPSC_HiC-loops-Arima-apa.png -p E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/5000 -I E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt -U
apa-analysis -O E14_RA_HiC_E14_RA_HiC-loops-Arima-apa.png -p E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/5000 -I E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt -U


##########TAD

######E14_2i_HiC-Arima-R1-filtered.mcool    E14_RA_HiC-Arima-R1-filtered.mcool
#E14_EPSC_HiC-Arima-R1-filtered.mcool

coolpup.py HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 HiC_mod-Arima.tad.bed --n_proc 20 --outname E14_ESC_HiC_ESC_tad_10kb --local --rescale

coolpup.py E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_2i_HiC-Arima.tad.bed --n_proc 20 --outname E14_2i_HiC_2i_tad_10kb   --local --rescale

coolpup.py E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC-Arima.tad.bed  --n_proc 20 --outname E14_EPSC_HiC_EPSC_tad_10kb --local --rescale

coolpup.py E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_RA_HiC-Arima.tad.bed --n_proc 20 --outname E14_RA_HiC_RA_tad_10kb --local --rescale

plotpup.py E14_ESC_HiC_ESC_tad_10kb E14_2i_HiC_2i_tad_10kb E14_EPSC_HiC_EPSC_tad_10kb E14_RA_HiC_RA_tad_10kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 0 --output ESC_2i_EPSC_RA_tad_10kb.png



coolpup.py ../../HiC_mod-Arima-R1-filtered.mcool::/resolutions/5000 HiC_mod-Arima-HICCUPS-loops.txt --n_proc 20 --outname E14_ESC_HiC_ESC_loops_5kb

coolpup.py ../../E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/5000 E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt --n_proc 20 --outname E14_2i_HiC_2i_loops_5kb

coolpup.py ../../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/5000 E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt  --n_proc 20 --outname E14_EPSC_HiC_EPSC_loops_5kb

coolpup.py ../../E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/5000 E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt --n_proc 20 --outname E14_RA_HiC_RA_loops_5kb

plotpup.py E14_ESC_HiC_ESC_loops_5kb E14_2i_HiC_2i_loops_5kb E14_EPSC_HiC_EPSC_loops_5kb E14_RA_HiC_RA_loops_5kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 3 --output ESC_2i_EPSC_RA_loops_5kb.png

###Produce the WashU format.
E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt
E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt
E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt



data<-read.table('E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt')
data1=data
ma=data.frame(data1[,c(1,2,3)],paste(data1[,4],':',data1[,5],'-',data1[,6],',',data1[,8],sep=''),1:dim(data1)[1],data1[,9])
colnames(ma)=c('chr','start','end','interacting','ID','strand')
ma1=data.frame(data1[,c(4,5,6)],paste(data1[,1],':',data1[,2],'-',data1[,3],',',data1[,8],sep=''),(1+dim(data1)[1]):(dim(data1)[1]+dim(data1)[1]),data1[,10])
colnames(ma1)=c('chr','start','end','interacting','ID','strand')

mas=rbind(ma,ma1)


write.table(mas,'E14_2i_HiC-Arima-HICCUPS-loops.5kb.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


sort  -k1,1 -k2,2n  E14_2i_HiC-Arima-HICCUPS-loops.5kb.bed > E14_2i_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
bgzip E14_2i_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
tabix -p bed E14_2i_HiC-Arima-HICCUPS-loops.5kb.sorted.txt.gz



data<-read.table('E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt')

data1=data

ma=data.frame(data1[,c(1,2,3)],paste(data1[,4],':',data1[,5],'-',data1[,6],',',data1[,8],sep=''),1:dim(data1)[1],data1[,9])
colnames(ma)=c('chr','start','end','interacting','ID','strand')
ma1=data.frame(data1[,c(4,5,6)],paste(data1[,1],':',data1[,2],'-',data1[,3],',',data1[,8],sep=''),(1+dim(data1)[1]):(dim(data1)[1]+dim(data1)[1]),data1[,10])
colnames(ma1)=c('chr','start','end','interacting','ID','strand')

mas=rbind(ma,ma1)


write.table(mas,'E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


sort  -k1,1 -k2,2n  E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.bed > E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
bgzip E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
tabix -p bed E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt.gz



data<-read.table('E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt')

data1=data

ma=data.frame(data1[,c(1,2,3)],paste(data1[,4],':',data1[,5],'-',data1[,6],',',data1[,8],sep=''),1:dim(data1)[1],data1[,9])
colnames(ma)=c('chr','start','end','interacting','ID','strand')
ma1=data.frame(data1[,c(4,5,6)],paste(data1[,1],':',data1[,2],'-',data1[,3],',',data1[,8],sep=''),(1+dim(data1)[1]):(dim(data1)[1]+dim(data1)[1]),data1[,10])
colnames(ma1)=c('chr','start','end','interacting','ID','strand')

mas=rbind(ma,ma1)


write.table(mas,'E14_RA_HiC-Arima-HICCUPS-loops.5kb_loops.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


sort  -k1,1 -k2,2n  E14_RA_HiC-Arima-HICCUPS-loops.5kb_loops.bed > E14_RA_HiC-Arima-HICCUPS-loops.5kb_loops.sorted.txt
bgzip E14_RA_HiC-Arima-HICCUPS-loops.5kb_loops.sorted.txt
tabix -p bed E14_RA_HiC-Arima-HICCUPS-loops.5kb_loops.sorted.txt.gz



data<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/HiC_mod-Arima-HICCUPS-loops.txt')
data1=data
#data1[,1]=c1
#data1[,4]=c4
ma=data.frame(data1[,c(1,2,3)],paste(data1[,4],':',data1[,5],'-',data1[,6],',',data1[,8],sep=''),1:dim(data1)[1],data1[,9])
colnames(ma)=c('chr','start','end','interacting','ID','strand')
ma1=data.frame(data1[,c(4,5,6)],paste(data1[,1],':',data1[,2],'-',data1[,3],',',data1[,8],sep=''),(1+dim(data1)[1]):(dim(data1)[1]+dim(data1)[1]),data1[,10])
colnames(ma1)=c('chr','start','end','interacting','ID','strand')

mas=rbind(ma,ma1)


write.table(mas,'E14_ESC_HiC-Arima-HICCUPS-loops.5kb.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


sort  -k1,1 -k2,2n  E14_ESC_HiC-Arima-HICCUPS-loops.5kb.bed > E14_ESC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
bgzip E14_ESC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt
tabix -p bed E14_ESC_HiC-Arima-HICCUPS-loops.5kb.sorted.txt.gz

####
#find the differential loops between 2i and EPSC.
E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt
E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt
E14_RA_HiC-Arima-HICCUPS-loops.5kb.txt

bedtools pairtopair -a E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt -b E14_2i_HiC-Arima-HICCUPS-loops.5kb.txt -type neither > EPSCVS2i_EPSC-specific-loops.bedpe;



data<-read.table('EPSCVS2i_EPSC-specific-loops.bedpe')

data1=data
data1[,1]=c1
data1[,4]=c4
ma=data.frame(data1[,c(1,2,3)],paste(data1[,4],':',data1[,5],'-',data1[,6],',',data1[,12],sep=''),1:dim(data1)[1],data1[,9])
colnames(ma)=c('chr','start','end','interacting','ID','strand')
ma1=data.frame(data1[,c(4,5,6)],paste(data1[,1],':',data1[,2],'-',data1[,3],',',data1[,12],sep=''),(1+dim(data1)[1]):(dim(data1)[1]+dim(data1)[1]),data1[,10])
colnames(ma1)=c('chr','start','end','interacting','ID','strand')

mas=rbind(ma,ma1)


write.table(mas,'EPSCVS2i_EPSC-specific-loops.bedpe.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


sort  -k1,1 -k2,2n  EPSCVS2i_EPSC-specific-loops.bedpe.bed > EPSCVS2i_EPSC-specific-loops.bedpe.sorted.txt
bgzip EPSCVS2i_EPSC-specific-loops.bedpe.sorted.txt
tabix -p bed EPSCVS2i_EPSC-specific-loops.bedpe.sorted.txt.gz

bedtools intersect -wa -wb -a EPSCVS2i_EPSC-specific-loops.bedpe.bed -b E14_EPSC_Esrrb-HiChIP_peaks.narrowPeak > EPSCVS2i_EPSC-specific-loops_Esrrb.bed
bedtools intersect -wa -wb -a E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.bed -b E14_EPSC_Esrrb-HiChIP_peaks.narrowPeak > EPSC-loops_Esrrb.bed

#4052 EPSCVS2i_EPSC-specific-loops.bedpe

####find  the ctcf motifs
findMotifsGenome.pl EPSCVS2i_EPSC-specific-loops.bedpe.bed mm9 MotifOutputDirectory/  -p 30 -find ctcf.motif > EPSCVS2i_EPSC-specific-loops_ctcfoutputfile.txt
findMotifsGenome.pl E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.bed mm9 MotifOutputDirectory_1/  -p 30 -find ctcf.motif > E14_EPSC_HiC-loops_ctcfoutputfile.txt


####Whether Loops anchors in EPSC were bound by Esrrb.


computeMatrix reference-point  --referencePoint center  -S ./E14_2i_Esrrb.HiChIP.mapped.PT.bw ./E14_EPSC_Esrrb.HiChIP.mapped.PT.bw  -R  EPSCVS2i_EPSC-specific-loops_Esrrb_peakpos.bed   --sortUsingSamples 1 --beforeRegionStartLength 5000  --afterRegionStartLength 5000  -o EPSC-Esrrb_matrix2.mat.gz -p 28  --samplesLabel 2i-Esrrb  EPSC-Esrrb  --missingDataAsZero --blackListFileName /home/jlab/Genome/mm9-blacklist.bed
#H3K27ac H3K4me1 H3K4me3 H3K27me3 ATAC-seq

plotHeatmap -m EPSC-Esrrb_matrix2.mat.gz --plotFileFormat pdf -out EPSC-Esrrb_matrix2_heatmap1.pdf   --heatmapHeight 6   --heatmapWidth 1 --colorList "white,forestgreen"

plotProfile -m  EPSC-Esrrb_matrix2.mat.gz -out EPSC-Esrrb_matrix2_plotprofile.pdf --perGroup  --plotHeight  5  --plotWidth  8


####aggregate peak analysis of Esrrb peaks called based on HiChIP using Hi-C data.


coolpup.py E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC_Esrrb-HiChIP_peaks.pos3.bed --n_proc 20 --outname E14_2i_HiC_EsrrbPeak_10kb   --local

coolpup.py E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC_Esrrb-HiChIP_peaks.pos3.bed  --n_proc 20 --outname E14_EPSC_HiC_EsrrbPeak_10kb --local

coolpup.py E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC_Esrrb-HiChIP_peaks.pos3.bed --n_proc 20 --outname E14_RA_HiC_EsrrbPeak_10kb --local

plotpup.py  E14_EPSC_HiC_EsrrbPeak_10kb E14_2i_HiC_EsrrbPeak_10kb E14_RA_HiC_EsrrbPeak_10kb --col_names   EPSC,2i,RA --n_cols 3 --enrichment 3 --output EPSC_2i_RA_2i_EsrrbPeak_10kb.png


####plot compartments.
fanc compartments E14_EPSC_HiC-Arima-R1-filtered.mcool@100000 E14_EPSC_HiC_fanc_100kb.ab
fanc compartments E14_2i_HiC-Arima-R1-filtered.mcool@100000 E14_2i_HiC_fanc_100kb.ab
fanc compartments HiC_mod-Arima-R1-filtered.mcool@100000 E14_ESC_HiC_fanc_100kb.ab
fanc compartments E14_RA_HiC-Arima-R1-filtered.mcool@100000 E14_RA_HiC_fanc_100kb.ab

fanc compartments -g  ~/Genome/STAR/mm9_fa/mm9.fa -v E14_EPSC_HiC_fanc_100kb.ab.ev_gc.txt  E14_EPSC_HiC_fanc_100kb.ab
fanc compartments -g  ~/Genome/STAR/mm9_fa/mm9.fa -v E14_2i_HiC_fanc_100kb.ab.ev_gc.txt  E14_2i_HiC_fanc_100kb.ab
fanc compartments -g  ~/Genome/STAR/mm9_fa/mm9.fa -v E14_RA_HiC_fanc_100kb.ab.ev_gc.txt  E14_RA_HiC_fanc_100kb.ab
fanc compartments -g  ~/Genome/STAR/mm9_fa/mm9.fa -v E14_ESC_HiC_fanc_100kb.ab.ev_gc.txt  E14_ESC_HiC_fanc_100kb.ab



fancplot -o fanc_example_100kb.EPSC.ab_chr17.png chr17 \
     -p square E14_EPSC_HiC_fanc_100kb.ab \
     -vmin -0.75 -vmax 0.75 -c RdBu_r


fancplot -o fanc_example_100kb.ab_chr17.ESC.png chr17 \
          -p square E14_ESC_HiC_fanc_100kb.ab \
          -vmin -0.75 -vmax 0.75 -c RdBu_r

fancplot -o fanc_example_100kb.ab_chr17.2i.png chr17 \
          -p square E14_2i_HiC_fanc_100kb.ab \
          -vmin -0.75 -vmax 0.75 -c RdBu_r

fancplot -o fanc_example_100kb.ab_chr17.RA.png chr17 \
          -p square E14_RA_HiC_fanc_100kb.ab \
          -vmin -0.75 -vmax 0.75 -c RdBu_r

fanc compartments -g ~/Genome/STAR/mm9_fa/mm9.fa \
                            -e E14_RA_HiC.100kb.ab_profile.png \
                            E14_RA_HiC-Arima-R1-filtered.mcool@100000 \
                            E14_RA_HiC_fanc_100kb.ab

###produce the insulation score.
fanc insulation E14_EPSC_HiC-Arima-R1-filtered.mcool@100000 \
                insulation/EPSC_100kb.insulation \
                -w 500000 1000000 --output-format bigwig


###produce the insulation score.
fanc insulation E14_2i_HiC-Arima-R1-filtered.mcool@100000 \
                insulation/2i_100kb.insulation \
                -w 500000 1000000 --output-format bigwig
###produce the insulation score.
fanc insulation E14_RA_HiC-Arima-R1-filtered.mcool@100000 \
                insulation/RA_100kb.insulation \
                -w 500000 1000000 --output-format bigwig
    ###produce the insulation score.
fanc insulation HiC_mod-Arima-R1-filtered.mcool@100000 \
                insulation/ESC_100kb.insulation \
                -w 500000 1000000 --output-format bigwig

###extend 20kb.
data<-read.table('E14_EPSC-Arima.tad.bed',sep='\t',stringsAsFactors=FALSE)
data1<-data.frame(chr=data[,1],start=data[,2]-20000,end=data[,2]+20000)
write.table(data1,'EPSC-Boundaries.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


data<-read.table('E14_RA_HiC-Arima.tad.bed',sep='\t',stringsAsFactors=FALSE)
data1<-data.frame(chr=data[,1],start=data[,2]-20000,end=data[,2]+20000)
write.table(data1,'RA-Boundaries.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

data<-read.table('E14_2i_HiC-Arima.tad.bed',sep='\t',stringsAsFactors=FALSE)
data1<-data.frame(chr=data[,1],start=data[,2]-20000,end=data[,2]+20000)
write.table(data1,'2i-Boundaries.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

data<-read.table('HiC_mod-Arima.tad.bed',sep='\t',stringsAsFactors=FALSE)
data1<-data.frame(chr=data[,1],start=data[,2]-20000,end=data[,2]+20000)
write.table(data1,'ESC-Boundaries.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


###identify the differential boundaries for each cell type.
mkdir 1;mergePeaks ESC-Boundaries.bed 2i-Boundaries.bed EPSC-Boundaries.bed  RA-Boundaries.bed -prefix 1/1 -venn 1/1.venn.txt

####APA pileup analysis of the cell-type-specific boundaries.
files<-c('1_ESC-Boundaries.bed','1_2i-Boundaries.bed','1_EPSC-Boundaries.bed','1_RA-Boundaries.bed')

for(i in 1:length(files)){
  data<-read.table(files[i],sep='\t',stringsAsFactors=FALSE)
  outname=strsplit(files[i],split='_')[[1]][2]
  outname=strsplit(outname,split='-')[[1]][1]

  data1<-data[,c(2:4)]
  write.table(data1,paste('../',outname,'_specific_boundaries.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


}

data1<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/EPSC_specific_boundaries.txt',sep='\t')

data2<-data.frame(data1[,1],data1[,2]+20000-1)
data2_names=paste(data2[,1],data2[,2],sep='|')
options(scipen = 100)
EPSC<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/E14_EPSC-Arima.tad.bed',sep='\t')
EPSC_names<-paste(EPSC[,1],EPSC[,2],sep='|')
n=match(data2_names,EPSC_names)
data3=EPSC[n,]
write.table(data3,'EPSC-specific-TADs.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

data1<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/ESC_specific_boundaries.txt',sep='\t')

data2<-data.frame(data1[,1],data1[,2]+20000-1)
data2_names=paste(data2[,1],data2[,2],sep='|')
options(scipen = 100)
EPSC<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/HiC_mod-Arima.tad.bed',sep='\t')
EPSC_names<-paste(EPSC[,1],EPSC[,2],sep='|')
n=match(data2_names,EPSC_names)
data3=EPSC[n,]
write.table(data3,'ESC-specific-TADs.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

data1<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/2i_specific_boundaries.txt',sep='\t')

data2<-data.frame(data1[,1],data1[,2]+20000-1)
data2_names=paste(data2[,1],data2[,2],sep='|')
options(scipen = 100)
EPSC<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/E14_2i_HiC-Arima.tad.bed',sep='\t')
EPSC_names<-paste(EPSC[,1],EPSC[,2],sep='|')
n=match(data2_names,EPSC_names)
data3=EPSC[n,]
write.table(data3,'2i-specific-TADs.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


data1<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/RA_specific_boundaries.txt',sep='\t')

data2<-data.frame(data1[,1],data1[,2]+20000-1)
data2_names=paste(data2[,1],data2[,2],sep='|')
options(scipen = 100)
EPSC<-read.table('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/E14_RA_HiC-Arima.tad.bed',sep='\t')
EPSC_names<-paste(EPSC[,1],EPSC[,2],sep='|')
n=match(data2_names,EPSC_names)
data3=EPSC[n,]
write.table(data3,'RA-specific-TADs.bed',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


coolpup.py ../HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 EPSC-specific-TADs.bed --n_proc 20 --outname E14_ESC_HiC_EPSC_specific_tad_10kb --local --rescale

coolpup.py ../E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 EPSC-specific-TADs.bed --n_proc 20 --outname E14_2i_HiC_EPSC_specific_tad_10kb   --local --rescale

coolpup.py ../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 EPSC-specific-TADs.bed --n_proc 20 --outname E14_EPSC_HiC_EPSC_specific_tad_10kb --local --rescale

coolpup.py ../E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 EPSC-specific-TADs.bed --n_proc 20 --outname E14_RA_HiC_EPSC_specific_tad_10kb --local --rescale

plotpup.py E14_ESC_HiC_EPSC_specific_tad_10kb E14_2i_HiC_EPSC_specific_tad_10kb E14_EPSC_HiC_EPSC_specific_tad_10kb E14_RA_HiC_EPSC_specific_tad_10kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 3 --output ESC_2i_EPSC_RA_EPSC_specific_tad_10kb.png


coolpup.py HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 HiC_mod-Arima.tad.bed --n_proc 20 --outname E14_ESC_HiC_ESC_tad_10kb --local --rescale

coolpup.py E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_2i_HiC-Arima.tad.bed --n_proc 20 --outname E14_2i_HiC_2i_tad_10kb   --local --rescale

coolpup.py E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC-Arima.tad.bed  --n_proc 20 --outname E14_EPSC_HiC_EPSC_tad_10kb --local --rescale

coolpup.py E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_RA_HiC-Arima.tad.bed --n_proc 20 --outname E14_RA_HiC_RA_tad_10kb --local --rescale

plotpup.py E14_ESC_HiC_ESC_tad_10kb E14_2i_HiC_2i_tad_10kb E14_EPSC_HiC_EPSC_tad_10kb E14_RA_HiC_RA_tad_10kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 0 --output ESC_2i_EPSC_RA_tad_10kb.png

coolpup.py HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 HiC_mod-Arima.tad.bed --n_proc 20 --outname E14_ESC_HiC_ESC_tad_10kb --local --rescale

coolpup.py E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_2i_HiC-Arima.tad.bed --n_proc 20 --outname E14_2i_HiC_2i_tad_10kb   --local --rescale

coolpup.py E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC-Arima.tad.bed  --n_proc 20 --outname E14_EPSC_HiC_EPSC_tad_10kb --local --rescale

coolpup.py E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_RA_HiC-Arima.tad.bed --n_proc 20 --outname E14_RA_HiC_RA_tad_10kb --local --rescale

plotpup.py E14_ESC_HiC_ESC_tad_10kb E14_2i_HiC_2i_tad_10kb E14_EPSC_HiC_EPSC_tad_10kb E14_RA_HiC_RA_tad_10kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 0 --output ESC_2i_EPSC_RA_tad_10kb.png


coolpup.py HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 HiC_mod-Arima.tad.bed --n_proc 20 --outname E14_ESC_HiC_ESC_tad_10kb --local --rescale

coolpup.py E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_2i_HiC-Arima.tad.bed --n_proc 20 --outname E14_2i_HiC_2i_tad_10kb   --local --rescale

coolpup.py E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_EPSC-Arima.tad.bed  --n_proc 20 --outname E14_EPSC_HiC_EPSC_tad_10kb --local --rescale

coolpup.py E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 E14_RA_HiC-Arima.tad.bed --n_proc 20 --outname E14_RA_HiC_RA_tad_10kb --local --rescale

plotpup.py E14_ESC_HiC_ESC_tad_10kb E14_2i_HiC_2i_tad_10kb E14_EPSC_HiC_EPSC_tad_10kb E14_RA_HiC_RA_tad_10kb --col_names   ESC,2i,EPSC,RA --n_cols 4 --enrichment 0 --output ESC_2i_EPSC_RA_tad_10kb.png


###########/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/insulation/








computeMatrix  reference-point  --referencePoint TSS  -S /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/insulation/EPSC_100kb.insulation_1mb.bigwig /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/insulation/ESC_100kb.insulation_1mb.bigwig /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/insulation/2i_100kb.insulation_1mb.bigwig /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/insulation/RA_100kb.insulation_1mb.bigwig  -R   EPSC-specific-TADs.bed --sortUsingSamples 1 --beforeRegionStartLength 1000000  --afterRegionStartLength 1000000 --skipZeros -o EPSC_specific_boundaries_matrix.mat.gz -p 28  --samplesLabel EPSC ESC 2i RA --blackListFileName /home/jlab/Genome/mm9-blacklist.bed --maxThreshold 200;

plotHeatmap -m EPSC_specific_boundaries_matrix.mat.gz --plotFileFormat pdf -out EPSC_specific_boundaries_matrix_heatmap1.pdf   --heatmapHeight 6   --heatmapWidth 1 --colorMap  colorwarm

plotProfile -m  EPSC_specific_boundaries_matrix.mat.gz -out EPSC_specific_boundaries_matrix_plotprofile.pdf --perGroup  --plotHeight  5  --plotWidth  8


###use cooltools to call compartments.
#first produce gc bedgraph as reference.
cooltools genome binnify ~/Genome/STAR/mm9/chrNameLength.txt 100000 > mm9.bin
cooltools genome gc mm9.bin  ~/Genome/STAR/mm9_fa/mm9.fa >gc.bedgraph
cooltools call-compartments --bigwig   --reference-track gc.bedgraph  -o E14_EPSC.compartments E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/100000
cooltools call-compartments --bigwig --reference-track gc.bedgraph  -o E14_2i.compartments E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/100000
cooltools call-compartments --bigwig --reference-track gc.bedgraph -o E14_RA.compartments E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/100000
cooltools call-compartments --bigwig --reference-track gc.bedgraph -o E14_ESC.compartments HiC_mod-Arima-R1-filtered.mcool::/resolutions/100000


############################33
cooler dump  ../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -b -H -m -o EPSC_HiC-Arima.matrix.gz
cooler dump  ../HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 -b -H -m -o ESC_HiC-Arima.matrix.gz
cooler dump ../E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -b -H -m -o 2i_HiC-Arima.matrix.gz
cooler dump  ../E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -b -H -m -o RA_HiC-Arima.matrix.gz

cooler dump --join ../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 > EPSC.10kb.txt
cooler dump --join ../HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000 > ESC.10kb.txt
cooler dump --join ../E14_2i_HiC-Arima-R1-filtered.mcool::/resolutions/10000 > 2i.10kb.txt
cooler dump --join ../E14_RA_HiC-Arima-R1-filtered.mcool::/resolutions/10000 > RA.10kb.txt

cooler dump  ../../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000 -b -m -o EPSC_HiC-Arima.matrix


#insulation score.
perl ~/Softwares/crane-nature-2015/scripts/matrix2insulation.pl -i EPSC_HiC-Arima.matrix.gz  -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v
perl ~/Softwares/crane-nature-2015/scripts/matrix2insulation.pl -i  ESC_HiC-Arima.matrix.gz  -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v
perl ~/Softwares/crane-nature-2015/scripts/matrix2insulation.pl -i  ../2i_HiC-Arima.matrix.gz  -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v
perl ~/Softwares/crane-nature-2015/scripts/matrix2insulation.pl -i  ../RA_HiC-Arima.matrix.gz  -is 500000 -ids 200000 -im mean -bmoe 3 -nt 0.1 -v


####use insulation2tad.pl;
perl ~/Softwares/cworld-dekker/scripts/perl/insulation2tads.pl -i insulationfile -b boundaryfile -o prefix

perl ~/Softwares/cworld-dekker/scripts/perl/matrix2insulation.pl
#####
tadtool plot chr7_145-155Mb.matrix.txt examples/chr12_20-35Mb_regions.bed chr12:31000000-33000000


###use arrowhead to call the TADs.
/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/filtered-mm9/
E14_2i_HiC-Arima-R1.hic    E14_ESC_HiC-Arima-R1.hic
E14_EPSC_HiC-Arima-R1.hic  E14_RA_HiC-Arima-R1.hic

java -Xmx40g  -jar  ~/Softwares/juicer_tools_1.13.02.jar  arrowhead -m 2000 -r 10000 -k KR --threads 20   /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/filtered-mm9/E14_EPSC_HiC-Arima-R1.hic EPSC.contact_domains.bedpe
java -Xmx40g  -jar  ~/Softwares/juicer_tools_1.13.02.jar  arrowhead -m 2000 -r 10000 -k KR --threads 20   /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/filtered-mm9/E14_2i_HiC-Arima-R1.hic 2i.contact_domains.bedpe

java -Xmx40g  -jar  ~/Softwares/juicer_tools_1.13.02.jar  arrowhead -m 2000 -r 10000 -k KR --threads 20   /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/filtered-mm9/E14_ESC_HiC-Arima-R1.hic ESC.contact_domains.bedpe

java -Xmx40g  -jar  ~/Softwares/juicer_tools_1.13.02.jar  arrowhead -m 2000 -r 10000 -k KR --threads 20   /media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/filtered-mm9/E14_RA_HiC-Arima-R1.hic RA.contact_domains.bedpe



############chr14:63200000-64200000
from tadlib.visualize.heatmaps import *
vis = Triangle('../../E14_EPSC_HiC-Arima-R1-filtered.mcool::/resolutions/10000', 'chr14', 63200000, 64200000)
vis.matrix_plot(vmin=0, vmax=0.03)
vis.plot_loops('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt',marker_size=10)
vis.plot_TAD('EPSC.tad.bed', linewidth=1.5)
vis.show()


####
from tadlib.visualize.heatmaps import *
vis = Triangle('../../HiC_mod-Arima-R1-filtered.mcool::/resolutions/10000', 'chr14', 63200000, 64200000)
vis.matrix_plot(vmin=0, vmax=0.03)
vis.plot_loops('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/HiC_mod-Arima-HICCUPS-loops.txt',marker_size=10)
vis.plot_TAD('ESC.tad.bed', linewidth=1.5)
vis.show()
