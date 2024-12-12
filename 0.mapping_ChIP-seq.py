import pandas as pd
import os,glob
from multiprocessing import Pool

def run(cmd):
    print(cmd)
    os.system(cmd)

samples=glob.glob('*_1*.fastq.gz')
for f in samples:
    a=f.split('_1')
    f2=a[0]+'_2'+a[1]
    fout=f.split('_1')[0]
    cmd='STAR --genomeDir /home/jlab/Genome/STAR/STAR_mm9/ --readFilesIn '+f+' '+f2+' --readFilesCommand zcat --runThreadN 10 --outFileNamePrefix ./Mapped/'+fout+' --outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 3 --outFilterMultimapNmax 500  --alignIntronMax 1 --alignEndsType EndToEnd '
    run(cmd)
