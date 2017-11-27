# Gandalf SNPs workflow

## Download VCFtools and Install where you want

https://urgi.versailles.inra.fr/download/gandalf/VCFtools-1.2.tar.gz

(set PATH in GandalfWorkflow.sh)

## Index your genome

bwa index data/genome.fasta

## Compute low complexity regions

python mdust_wrapper.py -i data/genome.fasta -v 28 -w 3 -m N -f bed -o data/lowcomplex.bed

## For each couple of paired end reads

./GandalfWorkflow.sh -i TEST -s fungi -f data/R1.fastq.gz -r data/R2.fastq.gz -g data/genome.fasta -t data/TE.bed -l data/lowcomplex.bed -v 2 

## Remove intermediate file (use prefix Index) and compress vcf

sh TEST_job.clean.sh
