#!/usr/bin/env bash

PICARD=/usr/local/bin/picard.jar
TRIMMOMATIC=/usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTERS=/usr/local/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
VCFFILTERING=/home/nlapalu/VCFtools-1.2/bin/VCFFiltering.py
export PYTHONPATH=/home/nlapalu/VCFtools-1.2
PWD=`pwd`


function show_usage {
	echo "
This script launches the Gandalf SNP workflow\

usage: GandalfWorkflow.sh [options]

    -h                help, print usage 
    -v 0|1|2|3        log level 0=off, 1=error, 2=info, 3=debug
    -i [ID]  Identifier (mandatory)
    -s [species]  Species (mandatory)
    -f [forward reads]  R1, in gz or not (mandatory)
    -r [reverse reads]  R2, in gz or not (mandatory)
    -g [genome] reference genome (mandatory)
    -t [TE bed file] TE bed file (mandatory)
    -l [low complexity bed file] low complexity bed file (mandatory)
	"
}

function log {
	case "$1" in
		error)
			level=1
			;;
		info)	
			level=2
			;;
		debug)
			level=3
			;;
	esac

	if [ $level -le $log_level ]
	then
		echo "###" "level: ${1}" "---" `date` "###"
		echo "###" "message: ${2}" "###"
	fi
}


## set options
while getopts "hv:i:s:f:r:g:t:l:" opt
do
	case "$opt" in
		h|\?)
			show_usage
			exit 0
			;;
		v)
			log_level=$OPTARG
			;;
		i)
                        # ID used as identifier
			ID=$OPTARG
			;;
		s)
                        # Species
			SPECIES=$OPTARG
			;;
		f)
                        # PE - R1 full path file (in gz or not)
			R1=$OPTARG
			;;
		r)
                        # PE - R2 full path file (in gz or not)
			R2=$OPTARG
			;;
		g)
                        # Reference genome
			REF=$OPTARG
			;;
		t)
                        # TE file
			TE=$OPTARG
			;;
		l)
                        # Low Complexity 
			LOW=$OPTARG
			;;
	esac
done
shift $((OPTIND-1))

## check options
if [ -z $log_level ]
then
	log_level=0
fi

if [ -z $ID ]
then
	echo "missing ID, set this option with -i"
	exit 0
fi

if [ -z $SPECIES ]
then
	echo "missing SPECIES, set this option with -s"
	exit 0
fi

if [[ -z $R1 || ! -f $R1 ]]
then
	echo "missing R1 (reads forward), set this option with -f"
	exit 0
fi

if [[ -z $R2 || ! -f $R2 ]]
then
	echo "missing R2 (reads forward), set this option with -r"
	exit 0
fi

if [[ -z $REF || ! -f $REF ]]
then
	echo "missing Reference (genome in fasta), set this option with -g"
	exit 0
fi

if [[ -z $TE || ! -f $TE ]]
then
	echo "missing TE file (in bed), set this option with -t"
	exit 0
fi

if [[ -z $LOW || ! -f $LOW ]]
then
	echo "missing Low Complexity (in bed), set this option with -l"
	exit 0
fi

if [[ -z $PICARD || ! -f $PICARD ]]
then
	echo "missing PICARD path or incorrect, set this option in file"
	exit 0
fi

if [[ -z $TRIMMOMATIC || ! -f $TRIMMOMATIC ]]
then
	echo "missing TRIMMOMATIC path or incorrect, set this option in file"
	exit 0
fi

if [[ -z $ADAPTERS || ! -f $ADAPTERS ]]
then
	echo "missing ADAPTERS path or incorrect, set this option in file"
	exit 0
fi

if [[ -z $VCFFILTERING || ! -f $VCFFILTERING ]]
then
	echo "missing VCFFILTERING path or incorrect, set this option in file"
	exit 0
fi




#if [[ -z "$dir_files" || ! -d "$dir_files" ]]
#then
#	echo "files directory not set (-f) or does not exist: ${dir_files}"
#	exit 0
#fi



log "info" "Workflow started !!!!!"

timestamp=`date "+%F"`
#dump=$dir_backup_dumps/$dbname-$timestamp.gz


# steps 
## clean reads trimmomatic
cmd1="java -jar $TRIMMOMATIC PE -threads 1 -phred33 $R1 $R2 $ID.R1.paired.gz $ID.R1.unpaired.gz $ID.R2.paired.gz $ID.R2.unpaired.gz ILLUMINACLIP:$ADAPTERS:2:30:12:1:true SLIDINGWINDOW:4:20 MINLEN:16"
## mapping
cmd2="bwa mem -t 1 -M $REF $ID.R1.paired.gz $ID.R2.paired.gz | samtools view -bS - | samtools sort -m 1G -o $ID.sorted.bam -"
## add RG
cmd3="java -Xms1024m -Xmx1024m -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT RGLB=$ID RGPL=illumina RGPU=gandalf I=$ID.sorted.bam O=RG.$ID.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE TMP_DIR=$ID RGSM=${ID}_${SPECIE} ID=${ID}_${SPECIES}"
## remove duplicates
cmd4="samtools rmdup RG.${ID}.sorted.bam dup.${ID}.sorted.bam"
## clean mapping
cmd5="samtools view -b -F 0x0100 -f 0x0002 -o min30.${ID}.sorted.bam -q 30 dup.${ID}.sorted.bam"
## extract list of reads
cmd6="samtools view min30.${ID}.sorted.bam | cut -f 1 | sort -T $PWD | uniq -c | grep ' 1 ' | cut -f8 -d' ' > min30.${ID}.lst"
## remove unpaired reads
cmd7="java -Xms1024m -Xmx1024m -jar $PICARD FilterSamReads I=min30.${ID}.sorted.bam FILTER=excludeReadList RLF=min30.${ID}.lst O=paired.min30.${ID}.sorted.bam"
## launch freebayes
cmd8="freebayes --report-monomorphic --ploidy 2 -X -u -f $REF paired.min30.${ID}.sorted.bam > min30.mono.${ID}.vcf"
## vcf filtering
cmd9="$VCFFILTERING -f min30.mono.${ID}.vcf -b ${TE} -b ${LOW} -o ${ID}.filtering.vcf"

cat > ${ID}_job.sh<<EOF
#$ -S /bin/bash
#$ -V
#$ -N ${ID}_job
#$ -cwd
#$ -o ${ID}_job.stdout
#$ -e ${ID}_job.stderr
#$ -pe parallel_smp 1
#$ -l mem_free=1G


echo "#### Trimming ####"
$cmd1
echo "#### Mapping ####"
$cmd2
echo "#### Read Group ####"
$cmd3
echo "#### Duplicates ####"
$cmd4
echo "#### Mapping filtering ####"
$cmd5
echo "#### Extracting list of unpaired reads ####"
$cmd6
echo "#### Remove unpaired reads  ####"
$cmd7
echo "#### Freebayes ####"
$cmd8
echo "#### VCFFiltering ####"
$cmd9
EOF

cat > ${ID}_job.clean.sh<<EOF
echo "removing raw mapping"
rm $ID.sorted.bam

echo "removing RG mapping"
rm RG.${ID}.sorted.bam
rm RG.${ID}.sorted.bai

echo "removing tmp dir"
rmdir ${ID}

echo "removing duplicate mapping"
rm dup.${ID}.sorted.bam

echo "removing unpaired bam"
rm min30.${ID}.lst
rm min30.${ID}.sorted.reads
rm min30.${ID}.sorted.bam
rm paired.min30.${ID}.sorted.reads

echo "removing non filtered vcf"
rm  min30.mono.${ID}.vcf

echo "compressing vcf file bgzip"
bgzip ${ID}.filtering.vcf

EOF


qsub ${ID}_job.sh

echo "At the end of the workflow launch ${ID}_job.clean.sh"
