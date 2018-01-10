#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh

P=8
SEQDATA=$BASEDIR/data/fastq

TMP=/tmp/$LSB_JOBID

mkdir -p $TMP

# Build bowtie indexes
mkdir -p $BASEDIR/data/bowtie_indexes
# Make symbolic link to genome fasta file
# (tophat looks for it in the folder where the indexes are)

if [ ! -f "$BASEDIR/data/bowtie_indexes/hg38.fa" ]; then
	ln -s $BASEDIR/data/hg38.fa $BASEDIR/data/bowtie_indexes/hg38.fa
fi
if [ ! -f "$BASEDIR/data/bowtie_indexes/mm10.fa" ]; then
	ln -s $BASEDIR/data/mm10.fa $BASEDIR/data/bowtie_indexes/mm10.fa
fi

if [ ! -f "$BASEDIR/data/bowtie_indexes/hg38.1.bt2" ]; then
	$BOWTIE2/bowtie2-build $BASEDIR/data/hg38.fa $BASEDIR/data/bowtie_indexes/hg38
fi

if [ ! -f "$BASEDIR/data/bowtie_indexes/mm10.1.bt2" ]; then
	$BOWTIE2/bowtie2-build $BASEDIR/data/mm10.fa $BASEDIR/data/bowtie_indexes/mm10
fi

# Create transcriptome-index
if [ ! -d "$BASEDIR/transcriptome_index/GencodeV21" ]; then
	mkdir -p $BASEDIR/transcriptome_index/GencodeV21
	$TOPHAT -p $P -G $BASEDIR/data/GencodeV21.gtf --transcriptome-index $BASEDIR/transcriptome_index/GencodeV21 $BASEDIR/data/bowtie_indexes/hg38
fi

if [ ! -d "$BASEDIR/transcriptome_index/GencodeM4" ]; then
	mkdir -p $BASEDIR/transcriptome_index/GencodeM4
	$TOPHAT -p $P -G $BASEDIR/data/GencodeM4.gtf --transcriptome-index $BASEDIR/transcriptome_index/GencodeM4 $BASEDIR/data/bowtie_indexes/mm10
fi


mkdir -p $BASEDIR/transcriptomes/tophat
mkdir -p $BASEDIR/transcriptomes/tophat/logs

# Define stranded datasets
STRANDED_MOUSE="SRR496249$|SRR496250$|SRR496221$|SRR496222$|SRR530639$|SRR530640$"
STRANDED_HUMAN="SRR307919$|SRR317040$|SRR307925$|SRR317039$|SRR307911$|SRR307912$|SRR307923$|SRR307924$"
STRANDED=$(echo "$STRANDED_MOUSE|$STRANDED_HUMAN")
DRY=""
ECHO=""
# Run tophat for mouse and human single end data
for i in $SEQDATA/single/hsa/*fastq.gz; do
	SAMP_NAME=$(basename $i .fastq.gz)
	if [[ $SAMP_NAME =~ $STRANDED ]]; then 
		LIBTYPE="fr-firststrand";
	else
		LIBTYPE="fr-unstranded";
	fi
	# Some samples require massive amounts of tmp disk space.
	# To avoid running out of space in /tmp we save the results directly 
	# to NFS (i.e., lsfTophat with -n option)
	if [[ $(du $i | awk '{print $1}') -gt 5243000 ]]; then
		NOTMP="-n"
	else
		NOTMP=""
	fi
	if [[ ! -d $BASEDIR/transcriptomes/tophat/hsa/$SAMP_NAME ]]; then
		$ECHO bsub -q research-rh6 -J $SAMP_NAME -oo $BASEDIR/transcriptomes/tophat/logs/$SAMP_NAME.log -M 50000 -n $P -R 'rusage[mem=50000]' $BIN/lsfTophat.sh $DRY $NOTMP -t $TOPHAT -p "--library-type $LIBTYPE -p $P --transcriptome-index $BASEDIR/transcriptome_index/GencodeV21/GencodeV21 --b2-sensitive --zpacker pigz" -o $BASEDIR/transcriptomes/tophat/hsa/$SAMP_NAME -i "$BASEDIR/data/bowtie_indexes/hg38" -f "$SEQDATA/single/hsa/$SAMP_NAME.fastq.gz"
	fi
done;

for i in $SEQDATA/single/mmu/*fastq.gz; do
	SAMP_NAME=$(basename $i .fastq.gz)
	if [[ $SAMP_NAME =~ $STRANDED ]]; then 
		LIBTYPE="fr-firststrand";
	else
		LIBTYPE="fr-unstranded";
	fi
	
	# SRR549335 and SRR549339 get stuck at the stage 
	# "Searching for junctions via segment mapping"
	# The option --no-coverage-search should solve
	# the problem.
	if [[ $SAMP_NAME =~ SRR549335$|SRR549339$ ]]; then
		NOCOVSEARCH="--no-coverage-search";
	else
		NOCOVSEARCH=""
	fi
	# Some samples require massive amounts of tmp disk space.
	# To avoid running out of space in /tmp we save the results directly 
	# to NFS (i.e., lsfTophat with -n option)
	if [[ $(du $i | awk '{print $1}') -gt 5243000 ]]; then
		NOTMP="-n"
	else
		NOTMP=""
	fi

	if [[ ! -d $BASEDIR/transcriptomes/tophat/mmu/$SAMP_NAME ]]; then
		$ECHO bsub -q research-rh6 -J $SAMP_NAME -oo $BASEDIR/transcriptomes/tophat/logs/$SAMP_NAME.log -M 50000 -n $P -R 'rusage[mem=50000]' $BIN/lsfTophat.sh $DRY $NOTMP -t $TOPHAT -p "--library-type $LIBTYPE -p $P $NOCOVSEARCH --transcriptome-index $BASEDIR/transcriptome_index/GencodeM4/GencodeM4 --b2-sensitive --zpacker pigz" -o $BASEDIR/transcriptomes/tophat/mmu/$SAMP_NAME -i "$BASEDIR/data/bowtie_indexes/mm10" -f "$SEQDATA/single/mmu/$SAMP_NAME.fastq.gz"
	fi
done;


# Run tophat for paired end data.
# The jobs require at least 100GB of free space in /tmp

for i in $SEQDATA/paired/hsa/*_1.fastq.gz; do
	SAMP_NAME=$(basename $i _1.fastq.gz)
	if [[ $SAMP_NAME =~ $STRANDED ]]; then 
		LIBTYPE="fr-firststrand";
	else
		LIBTYPE="fr-unstranded";
	fi
	# Some samples require massive amounts of tmp disk space.
	# To avoid running out of space in /tmp we save the results directly 
	# to NFS (i.e., lsfTophat with -n option)
	if [[ $(du $i | awk '{print $1}') -gt 5243000 ]]; then
		NOTMP="-n"
	else
		NOTMP=""
	fi
	if [[ ! -d $BASEDIR/transcriptomes/tophat/hsa/$SAMP_NAME ]]; then
		$ECHO bsub -q research-rh6 -J $SAMP_NAME -oo $BASEDIR/transcriptomes/tophat/logs/$SAMP_NAME.log -M 50000 -n $P -R 'rusage[mem=50000, tmp=50000]' $BIN/lsfTophat.sh $DRY $NOTMP -t $TOPHAT -p "--library-type $LIBTYPE -p $P --transcriptome-index $BASEDIR/transcriptome_index/GencodeV21/GencodeV21 --b2-sensitive --zpacker pigz" -o $BASEDIR/transcriptomes/tophat/hsa/$SAMP_NAME -i "$BASEDIR/data/bowtie_indexes/hg38" -f "$SEQDATA/paired/hsa/${SAMP_NAME}_1.fastq.gz $SEQDATA/paired/hsa/${SAMP_NAME}_2.fastq.gz"
	fi
done;

for i in $SEQDATA/paired/mmu/*_1.fastq.gz; do
	SAMP_NAME=$(basename $i _1.fastq.gz)
	if [[ $SAMP_NAME =~ $STRANDED ]]; then 
		LIBTYPE="fr-firststrand";
	else
		LIBTYPE="fr-unstranded";
	fi
	# SRR530639 and SRR530640 require massive amounts of tmp disk space.
	# To avoid running out of space in /tmp we save the results directly 
	# to NFS (i.e., lsfTophat with -n option)
	if [[ $(du $i | awk '{print $1}') -gt 5243000 ]]; then
		NOTMP="-n"
	else
		NOTMP=""
	fi

	if [[ ! -d $BASEDIR/transcriptomes/tophat/mmu/$SAMP_NAME ]]; then
		$ECHO bsub -q research-rh6 -J $SAMP_NAME -oo $BASEDIR/transcriptomes/tophat/logs/$SAMP_NAME.log -M 50000 -n $P -R 'rusage[mem=50000, tmp=50000]' $BIN/lsfTophat.sh $DRY $NOTMP -t $TOPHAT -p "--library-type $LIBTYPE -p $P --transcriptome-index $BASEDIR/transcriptome_index/GencodeM4/GencodeM4 --b2-sensitive --zpacker pigz" -o $BASEDIR/transcriptomes/tophat/mmu/$SAMP_NAME -i "$BASEDIR/data/bowtie_indexes/mm10" -f "$SEQDATA/paired/mmu/${SAMP_NAME}_1.fastq.gz $SEQDATA/paired/mmu/${SAMP_NAME}_2.fastq.gz"
	fi
done;

if [[ -z $ECHO ]]; then
	touch $BASEDIR/.map
fi
