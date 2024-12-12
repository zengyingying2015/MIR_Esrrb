import os,glob
from multiprocessing import Pool


files=glob.glob('*.bam')

def HTSeq_count(f):
    f_out=f.split('Aligned.sortedByCoord.out.bam')[0]
    cmd='htseq-count '+f+' mm9_ready.gtf -f bam -r pos -t exon -s no  > HTSeq/'+f_out
    print(cmd)
    os.system(cmd)



p = Pool(len(files))
p.map(HTSeq_count,files)
