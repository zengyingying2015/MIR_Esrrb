import os,glob
from multiprocessing import Pool
def run(cmd):
    print(cmd)
    os.system(cmd)

run('mkdir bigwig')
run('samtools index *.bam')
bams=glob.glob('*.bam')

for bam in bams:
    fout=bam.split('Align')[0]
    cmd='bamCoverage -b  '+bam+' -o bigwig/'+fout+'.bw --normalizeUsing RPKM  --effectiveGenomeSize 2620345972 --numberOfProcessors 20  --ignoreDuplicates      '
    print(cmd)
    os.system(cmd)
