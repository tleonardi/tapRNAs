#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


mkdir -p $BASEDIR/data

# Download chrominfo for mm10
if [ ! -f $BASEDIR/data/chromInfo_mm10.txt ]; then
    wget -qO - http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz | zcat > $BASEDIR/data/chromInfo_mm10.txt
fi

# Download chrominfo for hg19
if [ ! -f $BASEDIR/data/chromInfo_hg38.txt ]; then
	wget -qO - hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz | zcat > $BASEDIR/data/chromInfo_hg38.txt
fi

# Download hg38 to mm10 liftover chains
if [ ! -f $BASEDIR/data/hg38ToMm10.over.chain ]; then
	wget -qO - http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz | zcat > $BASEDIR/data/hg38ToMm10.over.chain
fi

# Download hg19 to hg38 liftover chains
if [ ! -f $BASEDIR/data/hg19ToHg38.over.chain ]; then
	wget -qO - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | zcat > $BASEDIR/data/hg19ToHg38.over.chain
fi

# Download mouse genome in 2bit and convert to fasta
if [ ! -f $BASEDIR/data/mm10.2bit ]; then
	wget -O - ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit > $BASEDIR/data/mm10.2bit
fi

if [ ! -f $BASEDIR/data/mm10.fa ]; then
	$twoBitToFa $BASEDIR/data/mm10.2bit $BASEDIR/data/mm10.fa
fi

# Download human genome in 2bit and convert to fasta
if [ ! -f $BASEDIR/data/hg38.2bit ]; then
	wget -O - ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit > $BASEDIR/data/hg38.2bit
fi

if [ ! -f $BASEDIR/data/hg38.fa ]; then
	$twoBitToFa $BASEDIR/data/hg38.2bit $BASEDIR/data/hg38.fa
fi

# Download Gencode V21 BED file
if [ ! -f $BASEDIR/data/GencodeV21.bed ]; then
	wget -O - "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz" | gunzip  > $BASEDIR/data/GencodeV21.gtf
	$gtfToGenePred -ignoreGroupsWithoutExons -infoOut=$BASEDIR/data/GencodeV21_info.txt $BASEDIR/data/GencodeV21.gtf $BASEDIR/data/GencodeV21.genePred
	$genePredToBed $BASEDIR/data/GencodeV21.genePred $BASEDIR/data/GencodeV21_unsorted.bed
	awk 'BEGIN{OFS=FS="\t"}$1!=$2{sub("\\.[0-9]+$", "", $1); sub("\\.[0-9]+$", "", $2); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10'} $BASEDIR/data/GencodeV21_info.txt > $BASEDIR/data/GencodeV21_info_tr.txt
	awk 'BEGIN{OFS=FS="\t"}{sub("\\.[0-9]+$", "", $4); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $BASEDIR/data/GencodeV21_unsorted.bed > $BASEDIR/data/GencodeV21_unsorted_shortID.bed
	sort -k1,1 -k2,2n $BASEDIR/data/GencodeV21_unsorted_shortID.bed >$BASEDIR/data/GencodeV21.bed
	ln -s $BASEDIR/data/GencodeV21.bed $BASEDIR/data/Gencode-current_hsa.bed
fi

# Download Gencode M4
if [ ! -f $BASEDIR/data/GencodeM4.bed ]; then

	wget -O - "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz" | gunzip > $BASEDIR/data/GencodeM4.gtf
	$gtfToGenePred -ignoreGroupsWithoutExons -infoOut=$BASEDIR/data/GencodeM4_info.txt $BASEDIR/data/GencodeM4.gtf $BASEDIR/data/GencodeM4.genePred
	$genePredToBed $BASEDIR/data/GencodeM4.genePred $BASEDIR/data/GencodeM4_unsorted.bed
	awk 'BEGIN{OFS=FS="\t"}$1!=$2{sub("\\.[0-9]+$", "", $1); sub("\\.[0-9]+$", "", $2); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10'} $BASEDIR/data/GencodeM4_info.txt > $BASEDIR/data/GencodeM4_info_tr.txt
	awk 'BEGIN{OFS=FS="\t"}{sub("\\.[0-9]+$", "", $4); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $BASEDIR/data/GencodeM4_unsorted.bed > $BASEDIR/data/GencodeM4_unsorted_shortID.bed
	sort -k1,1 -k2,2n $BASEDIR/data/GencodeM4_unsorted_shortID.bed > $BASEDIR/data/GencodeM4.bed
	ln -s $BASEDIR/data/GencodeM4.bed $BASEDIR/data/Gencode-current_mmu.bed
fi

# Download Gencode lincRNA annotation (human)
if [ ! -f $BASEDIR/data/Gencode.v21.lincRNAs.bed ]; then
	wget -O - "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.long_noncoding_RNAs.gtf.gz" | gunzip > $BASEDIR/data/Gencode.v21.lincRNAs.gtf
	$gtfToGenePred -ignoreGroupsWithoutExons $BASEDIR/data/Gencode.v21.lincRNAs.gtf $BASEDIR/data/Gencode.v21.lincRNAs.genePred
	$genePredToBed $BASEDIR/data/Gencode.v21.lincRNAs.genePred $BASEDIR/data/Gencode.v21.lincRNAs_unsorted.bed
	awk 'BEGIN{OFS=FS="\t"}{sub("\\.[0-9]+$", "", $4); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' $BASEDIR/data/Gencode.v21.lincRNAs_unsorted.bed > $BASEDIR/data/Gencode.v21.lincRNAs_unsorted_shortID.bed
	sort -k1,1 -k2,2n $BASEDIR/data/Gencode.v21.lincRNAs_unsorted_shortID.bed > $BASEDIR/data/Gencode.v21.lincRNAs.bed
fi


# Download human repeat masker
if [[ ! -f $BASEDIR/data/repeatMasker.txt ]]; then
	wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz | gunzip > $BASEDIR/data/repeatMasker.txt
fi

# Download Ensembl orthology annotation
if [ ! -f $BASEDIR/data/hg38ToMm10_orthologs_mart.txt ]; then
	wget -O $BASEDIR/data/hg38ToMm10_orthologs_mart.txt "http://www.ensembl.org/biomart/martservice?query=$(cat $DATA/biomart_query_orthology.xml | tr '\n' ' ')"
fi
if [ ! -f $BASEDIR/data/hg38ToZz10_orthologs_mart.txt ]; then
	wget -O $BASEDIR/data/hg38ToZz10_orthologs_mart.txt "http://www.ensembl.org/biomart/martservice?query=$(cat $DATA/biomart_query_orthology_zebrafish.xml| tr '\n' ' ')"
fi



touch $BASEDIR/.download
