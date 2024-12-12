library(GENOVA)

#/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9
###E14_2i_HiC-Arima-R1-filtered.mcool    HiC_mod-Arima-R1-filtered.mcool
#E14_EPSC_HiC-Arima-R1-filtered.mcool  HiC_std-Arima-R1-filtered.mcool
#E14_RA_HiC-Arima-R1-filtered.mcool

#####> brewer.pal(9,'Set1')
[1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628"
[8] "#F781BF" "#999999"



EPSC <- load_contacts(signal_path = '/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/E14_EPSC_HiC-Arima-R1-filtered.mcool',resolution=10000,sample_name = "EPSC",balancing = T,colour = "#E41A1C")
ESC <- load_contacts(signal_path = '/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/HiC_mod-Arima-R1-filtered.mcool',resolution=10000,sample_name = "ESC",balancing = T,colour = "black")
RA <- load_contacts(signal_path = '/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/E14_RA_HiC-Arima-R1-filtered.mcool',resolution=10000,sample_name = "RA",balancing = T,colour = "#377EB8")
E2i <- load_contacts(signal_path = '/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/E14_2i_HiC-Arima-R1-filtered.mcool',resolution=10000,sample_name = "2i",balancing = T,colour = "#4DAF4A")

# 'chr14', 63200000, 64200000
ID <- insulation_domainogram(
EPSC,
chrom = 'chr14',
start = 63200000,
end =  64200000,
window_range = c(1, 101),
step = 2
)

library(extrafont)
pdf('EPSC_chr14_boundaries_Gata4.pdf',family='Arial',width=5,height=1.5)

visualise(ID)
dev.off()


ID1 <- insulation_domainogram(
ESC,
chrom = 'chr14',
start = 63200000,
end =  64200000,
window_range = c(1, 101),
step = 2
)
library(extrafont)
pdf('ESC_chr14_boundaries_Gata4.pdf',family='Arial',width=5,height=1.5)

visualise(ID1)
dev.off()


###show the insulation score difference.
EPSC_ESC_10kb_insulation <- insulation_score(
list(ESC, EPSC),
window = 25
)

library(ggplot2)
library(RColorBrewer)

pdf('EPSCVSESC_chr14_boundaries_Gata4_IS_difference.pdf',family='Arial',width=6,height=1.5)
visualise(EPSC_ESC_10kb_insulation,
chrom = 'chr14',
start = 63200000,
end =  64200000,
contrast = 1)
dev.off()



boundaries=read.delim('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/IS/EPSC_10kb.insulation.boundaries',h=F)
ESC_boundaries=read.delim('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/IS/ESC_10kb.insulation.boundaries',h=F)


cell_10kb_insulation <- insulation_score(
list(ESC,E2i, EPSC,RA),
window = 50
)

tornado_insulation(cell_10kb_insulation, ESC_TADs, bed_pos = 'start')
pdf('ESC-TADs_ID_tornado_insulation.pdf',family='Arial',width=6,height=6)
tornado_insulation(cell_10kb_insulation, ESC_TADs, bed_pos = 'start')
dev.off()

pdf('ESC-TADs_ID_tornado_insulation_profile_window50.pdf',family='Arial',width=6,height=3)
tornado_insulation(cell_10kb_insulation, ESC_TADs, bed_pos = 'start',mode = 'profile')
dev.off()


ESC_TADs = read.delim('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/HiC_mod-Arima.tad.bed', h = F)

EPSC_TADs = read.delim('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/diffBoundaries/E14_EPSC-Arima.tad.bed', h = F)
EPSC_Loops = read.delim('/media/jlab/JL23/ZYY/Esrrb_Project/Hi-C/coolers-mm9/E14_EPSC_HiC-Arima-HICCUPS-loops.5kb.txt', h = F, skip = 1)



hic_matrixplot(exp1 =EPSC,
chrom = 'chr14',
start = 63200000,
end =  64200000,
loops =EPSC_Loops,
tads = EPSC_TADs, # see ATA
tads.type = 'upper', # only plot in lower triangle
tads.colour = '#91cf60', # green TAD-borders
cut.off = 50, # upper limit of contacts
#skipAnn = T) # skip the outside annotation
plot(ID, minimalist = TRUE)
