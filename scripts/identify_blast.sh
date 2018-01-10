#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


P=16
TMP=/tmp/$LSB_JOBID

# Prepare tmp folder and download data
mkdir -p $TMP


# Get sequence for 500bp promoter of each human ncRNA
awk 'BEGIN{OFS="\t"}{if ($6=="+") {print $1,$2-500,$2,$4,$5,$6} else {print $1,$3,$3+500,$4,$5,$6}}' $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed | sort -k1,1 -k2,2n> $BASEDIR/nc/GCode_nc_splice_5pEnd_prom.bed
$PSTR getDna $BASEDIR/data/hg38.fa $BASEDIR/nc/GCode_nc_splice_5pEnd_prom.bed > $BASEDIR/nc/GCode_nc_splice_5pEnd_prom.fa


# This block checks the conservation of pcRNAs in mouse and human.
# the file produced ($BASEDIR/nc/ncRNAcmb_mouseCons.bedx) will be
# used by annotate.sh for the final annotation.

# Get sequence for ncRNAs to see if they are conserved
$PSTR getDna $BASEDIR/data/hg38.fa $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.fa

# Make a mask file based on the soft-masking of the genome
# and then create a BLAST DB
$CONVERT2BLASTMASK -in $BASEDIR/data/mm10.fa -masking_algorithm repeat -masking_options "repeatmasker and tandem repeats from UCSC" -outfmt maskinfo_asn1_bin -out $BASEDIR/data/mm10_mask.asnb 
$MAKEBLASTDB -in $BASEDIR/data/mm10.fa -mask_data $BASEDIR/data/mm10_mask.asnb -dbtype nucl


$PSTR blasttranscript -tempDir $TMP -useBlastPlus -p $P -o $BASEDIR/coding/hsaNcIncNovVsMm10.blast $BASEDIR/data/mm10.fa $BASEDIR/mouse/transcripts/CombinedRefAnddeNovo.bed $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.fa >$BASEDIR/nc/Gencode-currentncPlusNovelInMouse.bedx

# Same for allmRNA
$PSTR blasttranscript -tempDir $TMP -useBlastPlus -p $P -o $BASEDIR/coding/hsaNcIncNovVsMm10.blast $BASEDIR/data/mm10.fa $BASEDIR/data/mm10_allmrna_canonical.bed $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.fa >$BASEDIR/nc/Gencode-currentPlusNovelInMouse_allmrna.bedx
#Create list of combined transcripts from working set and allmrna
sort -k4,4 -k13,13 $BASEDIR/nc/Gencode-currentncPlusNovelInMouse.bedx | awk 'BEGIN{OFS="\t"}{print $4,$13,$14*$15/100,$19-$18}' | groupBy -g 1,2 -c 3,4 -o sum,sum | awk 'BEGIN{OFS="\t"}{pc=sprintf("%.2f",$3/$4);print $1,$2,pc,$4}' >  $BASEDIR/nc/Gencode-currentncPlusNovelInMouse.txt
sort -k4,4 -k13,13 $BASEDIR/nc/Gencode-currentPlusNovelInMouse_allmrna.bedx | awk 'BEGIN{OFS="\t"}{print $4,$13,$14*$15/100,$19-$18}' | groupBy -g 1,2 -c 3,4 -o sum,sum | awk 'BEGIN{OFS="\t"}{pc=sprintf("%.2f",$3/$4);print $1,$2,pc,$4}' >  $BASEDIR/nc/Gencode-currentPlusNovelInMouse_allmrna.txt
cat $BASEDIR/nc/Gencode-currentncPlusNovelInMouse.txt $BASEDIR/nc/Gencode-currentPlusNovelInMouse_allmrna.txt | sort -k1,1 -k4,4nr | sort -k1,1 -u | awk '$4>50' > $BASEDIR/nc/Gencode-currentncPlusNovelInMouse_cmb.txt
cut -f1-12 $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed | pstr join -x -idCol 1 -columns 2,3,4 $BASEDIR/nc/Gencode-currentncPlusNovelInMouse_cmb.txt stdin | awk 'BEGIN{OFS="\t"}{LEN=0;split($11,SZ,",");for (i in SZ){LEN+=SZ[i]};pc=sprintf("%.2f",$15/LEN); print $0,pc}' > $BASEDIR/nc/ncRNAcmb_mouseCons.bedx

 


### Feature association - ncRNA to closest coding gene in human
# For each ncRNA find the closest upstream and downstream (and overlapping) protein coding genes
# Field 13 indicates the name of the ncRNA
$PSTR closest -p $P $BASEDIR/coding/Gencode-current_hsa_coding.bed $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_closestCoding.bedx

# Join the closest coding transcripts to their gene ids
$PSTR join -idCol 1 --columns 2 $BASEDIR/data/GencodeV21_info_tr.txt $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_closestCoding.bedx > $BASEDIR/nc/hg38ncRNAToClosestCodingGene.bedx

# Join the non-coding BED line
# Awk swaps coding and nc coords. Final out is:
# 1-5: Coding coords
# 6-10: nc coords
# 11-15: distances between coding and nc (pstr closest)
# 16: Gene id of coding
join -1 4 -2 13 -t $'\t' <(cut -f1-6 $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed | sort -k4,4) <(sort -k13,13 $BASEDIR/nc/hg38ncRNAToClosestCodingGene.bedx) | awk 'BEGIN{OFS=FS="\t"}{print $7,$8,$9,$10,$12,$2,$3,$4,$1,$6,$19,$20,$21,$22,$23,$24}'> $BASEDIR/nc/hg38ncRNAToClosestCodingGene_withNcAnnot.bedx

# Do positional annotation of nc

function posAnnot() {
	awk ' BEGIN{OFS="\t";BDDIST=2000}
		{if($2<$8 && $3>$7){ # Transcripts overlap
			if ($5==$10) { # Same strand -> OLAP
				if ($5=="+") print $9,"OLAP",$2-$7,$10,$11,$12,$13,$14,$15,$16;
				else print $9,"OLAP",$8-$3,$10,$11,$12,$13,$14,$15,$16;
			} 
			else {		# Different strand -> Asense
				if ($5=="+") print $9,"Asense",$8-$2,$10,$11,$12,$13,$14,$15,$16;
				else print $9,"Asense",$3-$7,$10,$11,$12,$13,$14,$15,$16;
			}
		}
		else if($5!=$10){
			if ($5=="+") {
				if ($8<$2){
					if (($2-$8)<=BDDIST) print $9,"BIDIR",$2-$8,$10,$11,$12,$13,$14,$15,$16;
					else print $9,"US-AS",$8-$2,$10,$11,$12,$13,$14,$15,$16;
				} 
				else if ($3<=$7) print $9,"DS-AS",$2-$8,$10,$11,$12,$13,$14,$15,$16;
			};
			if ($5=="-" ){
				if($7>$3){
					if (($7-$3)<=BDDIST) print $9,"BIDIR",$7-$3,$10,$11,$12,$13,$14,$15,$16;
					else print $9,"US-AS",$7-$3,$10,$11,$12,$13,$14,$15,$16;
				} 
				else if ($8<=$2) print $9,"DS-AS",$3-$7,$10,$11,$12,$13,$14,$15,$16;
			} 
		}
		else if($5==$10){
			if ($5=="+" ) {
				if ($8<$2) print $9,"US-S",$2-$7,$10,$11,$12,$13,$14,$15,$16;
				else if ($7>$3) print $9,"DS-S",$7-$2,$10,$11,$12,$13,$14,$15,$16;
			};
			if ($5=="-") {
				if ($7>$3)  print $9,"US-S",$8-$3,$10,$11,$12,$13,$14,$15,$16;
				else if ($8<$2) print $9,"DS-S",$3-$8,$10,$11,$12,$13,$14,$15,$16;
			}
		}
	}' $1 | awk 'function abs(value) {return (value<0?-value:value)}BEGIN{OFS="\t"}{print $1,$2,$3,abs($3),$10}'
}
	
posAnnot $BASEDIR/nc/hg38ncRNAToClosestCodingGene_withNcAnnot.bedx > $BASEDIR/nc/hg38ncRNAToClosestCodingGene_withNcAnnot_withOrientation.bedx

# For each pcRNA-Coding gene pair in a given orientation keep only the closest line
sort -k1,1 -k5,5 -k2,2 -k4,4n $BASEDIR/nc/hg38ncRNAToClosestCodingGene_withNcAnnot_withOrientation.bedx | sort -k1,1 -k5,5 -k2,2 -u > $BASEDIR/nc/hg38ncRNAToClosestCodingGene.txt


# Do the same of each mouse ncRNA
$PSTR closest -p $P $BASEDIR/coding/Gencode_mm10_coding.bed $BASEDIR/nc/mm10CombinedNCTranscripts.bed > $BASEDIR/nc/mm10CombinedNCTranscripts_closestCoding.bedx
$PSTR closest -p $P $BASEDIR/coding/Gencode_mm10_coding.bed $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed > $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID_closestCoding.bedx

# Merge closest gene annotation for CombinedNCTranscripts and allmRNAs
cat $BASEDIR/nc/mm10CombinedNCTranscripts_closestCoding.bedx $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID_closestCoding.bedx | sort -k1,1 -k2,2n -k3,3n > $BASEDIR/nc/mm10CombinedNCTranscripts+allmrna_spliced_NC_uniqID_closestCoding.bedx

# Join the closest coding transcripts to their gene ids
$PSTR join -idCol 1 --columns 2 $BASEDIR/data/GencodeM4_info_tr.txt $BASEDIR/nc/mm10CombinedNCTranscripts+allmrna_spliced_NC_uniqID_closestCoding.bedx > $BASEDIR/nc/mm10ncRNAToClosestCodingGene.txt

# Join the non-coding BED line
# Awk swaps coding and nc coords. Final out is:
# 1-5: Coding coords
# 6-10: nc coords
# 11-15: distances between coding and nc (pstr closest)
# 16: Gene id of coding
join -1 4 -2 13 -t $'\t' <(cut -f1-6 $BASEDIR/nc/mm10CombinedNCTranscripts+allmrna_spliced_NC_uniqID.bed | sort -k4,4) <(sort -k13,13 $BASEDIR/nc/mm10ncRNAToClosestCodingGene.txt) | awk 'BEGIN{OFS=FS="\t"}{print $7,$8,$9,$10,$12,$2,$3,$4,$1,$6,$19,$20,$21,$22,$23,$24}' > $BASEDIR/nc/mm10ncRNAToClosestCodingGene_withNcAnnot.bedx

# Do positional annotation of nc in mouse
posAnnot $BASEDIR/nc/mm10ncRNAToClosestCodingGene_withNcAnnot.bedx  >  $BASEDIR/nc/mm10ncRNAToClosestCodingGene_withNcAnnot_withOrientation.bedx

# Annotate the human syntenic of the coding gene closeset to each ncRNA
join -1 5 -2 3 -t $'\t' <(sort -k5,5 $BASEDIR/nc/mm10ncRNAToClosestCodingGene_withNcAnnot_withOrientation.bedx) <(sort -k3,3 $BASEDIR/data/hg38ToMm10_orthologs_mart.txt) | cut -f2,3,4,5,6 | sort -u> $BASEDIR/nc/mm10ncRNAToClosestCodingGene_humanSynt.txt


### Blast human promoter regions against mm10
# outfmt 6 indicates the output format. The columns are as follows:
# 	qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
$BLASTN -task blastn -db $BASEDIR/data/mm10.fa -out $BASEDIR/nc/hsaNcVsMm10_2.blast  -query $BASEDIR/nc/GCode_nc_splice_5pEnd_prom.fa -outfmt 6 -evalue 0.001 -num_threads $P -db_soft_mask 40 -lcase_masking

### Clean and filter the results of blastn of ncRNA promoter from human to mouse
# we keep only promoter with match >100nt and Eval<1e-10. We then print in BED format the mouse coords. Convert blast (1 based) results back into bed format (0 based)
# We keep multimapping promoters
awk 'BEGIN{OFS="\t"}$4>=100 && $11<1E-10 {if ($10>$9) {print $2,$9-1,$10,$1,0,"."} else {print $2,$10-1,$9,$1,0,"."}}' $BASEDIR/nc/hsaNcVsMm10_2.blast > $BASEDIR/nc/hsaNcVsMm10_2_gt100.bed


rm -rf $TMP
touch $BASEDIR/.identify_blast

