#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


P=16


#########################################################
#							#
#			HUMAN				#
#							#
#########################################################
mkdir -p $BASEDIR/data
### Select spliced noncoding RNA from human
mkdir -p $BASEDIR/nc
awk '$7==$8 && $10>1' $BASEDIR/data/Gencode-current_hsa.bed > $BASEDIR/nc/Gencode-current_hsa_ncSplice.bed
awk '$10>1' $BASEDIR/data/Gencode-current_hsa.bed > $BASEDIR/nc/Gencode-current_hsa_Spliced.bed

# Annotate all pseudogenes for later filtering
awk '$3~"transcript"{if($14 ~ "IG_C_pseudogene" || $14 ~  "IG_J_pseudogene" || $14 ~  "IG_V_pseudogene" || $14 ~  "polymorphic_pseudogene" || $14 ~  "processed_pseudogene" || $14 ~  "pseudogene" || $14 ~  "transcribed_processed_pseudogene" || $14 ~  "transcribed_unitary_pseudogene" || $14 ~  "transcribed_unprocessed_pseudogene" || $14 ~  "translated_processed_pseudogene" || $14 ~  "translated_unprocessed_pseudogene" || $14 ~  "TR_J_pseudogene" || $14 ~  "TR_V_pseudogene" || $14 ~  "unitary_pseudogene" || $14 ~  "unprocessed_pseudogene"){print $12}}' $BASEDIR/data/GencodeV21.gtf | perl -pe 's/[\";]//g' | perl -pe 's/\.[0-9]+$//' > $BASEDIR/nc/Gencode-current-pseudogenes.txt

$PSTR join -x $BASEDIR/nc/Gencode-current-pseudogenes.txt $BASEDIR/data/Gencode-current_hsa.bed | cut -f1-12 > $BASEDIR/nc/Gencode-current-pseudogenes.bed


# Create reference file containing only coding exons
mkdir -p $BASEDIR/coding
awk '$7!=$8' $BASEDIR/data/Gencode-current_hsa.bed > $BASEDIR/coding/Gencode-current_hsa_coding_with_untracked.bed

# If there are two identical Gencode transcripts for some reason
# cufflinks only tracks one. To make sure that we use the coding
# that is tracked, we de-duplicate Gencode based on the output of
# cuffmerge.
# To do so:
# 1: Extract IDs of tracked transcripts
# 2: join tracked trascript IDs with gencode
# 3: filter duplicated transcripts keeping the tracked one

cut -f14 $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bedx | perl -pe 's/(ENS.+)\.[0-9]+$/$1/' | awk '{print $1"\t1"}' > $BASEDIR/transcriptomes/cuffmerge/hsa/transcripts_in_merged.txt
join -1 4 -2 1 -a 1 -t $'\t' <(sort -k4,4  $BASEDIR/coding/Gencode-current_hsa_coding_with_untracked.bed) <(sort -k1,1 $BASEDIR/transcriptomes/cuffmerge/hsa/transcripts_in_merged.txt) | awk 'BEGIN{OFS=FS="\t"}{if($13!=1){$13=0};print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sort -k1,1 -k2,2n -k3,3n > $BASEDIR/coding/Gencode-current_hsa_coding_isTracked.bedx
cut -f1-12 $BASEDIR/coding/Gencode-current_hsa_coding_isTracked.bedx > $BASEDIR/coding/Gencode-current_hsa_coding.bed

$PSTR getCoding -p $P $BASEDIR/coding/Gencode-current_hsa_coding.bed > $BASEDIR/coding/Gencode-current_hsa_codingExons.bed


### Only those with no sense overlap to coding exon
$OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping $BASEDIR/coding/Gencode-current_hsa_codingExons.bed $BASEDIR/nc/Gencode-current_hsa_ncSplice.bed $BASEDIR/nc/Gencode-current_hsa_ncSplice_noCodeOlap.bed

###Integrate denovo transcriptomes. Check only spliced transcripts with no coding exon overlap and not similar to existing
# Process the merged transcriptomes of all human tissues (cuffmerge) and:
# 	1) Select those with more than 1 exon and
#	2) Discard transcript in supercontigs ($1!~/_/) and
# 	3) less than 20bp of overlap with a coding exon and
# 	4) less than 50% sense overlap with a spliced Gencode transcript (if the transcript overlaps an unspliced Gencode transcript we keep it, because it's spliced anyway

mkdir -p $BASEDIR/human/transcripts/
awk '$10>1 && $1!~/_/ && $1!="chrM"' $BASEDIR/transcriptomes/cuffmerge/hsa/merged.bed | $OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping -inFmt=bed $BASEDIR/coding/Gencode-current_hsa_codingExons.bed stdin stdout | $INTERSECT_BED -f 0.5 -s -split -v -a - -b $BASEDIR/nc/Gencode-current_hsa_Spliced.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $BASEDIR/human/transcripts/allNocode.bed

# Remove redundancy
$PSTR dedup -p $P --exEncomp $BASEDIR/human/transcripts/allNocode.bed > $BASEDIR/human/transcripts/allNocode_nr.bed

# Calculate coding potential with CPAT
bedtools getfasta -split -name -s -fi $BASEDIR/data/hg38.fa -bed $BASEDIR/human/transcripts/allNocode_nr.bed -fo $BASEDIR/human/transcripts/allNocode_nr.fa
cpat.py -g $BASEDIR/human/transcripts/allNocode_nr.fa -x $DATA/prebuilt_cpat_models/Human_Hexamer.tab -d $DATA/prebuilt_cpat_models/Human_logitModel.RData -o $BASEDIR/human/transcripts/allNocode_nr.cpat

# remove header line and sort by id
tail -n +2 $BASEDIR/human/transcripts/allNocode_nr.cpat | sort -k1,1 > $BASEDIR/human/transcripts/allNocode_nr.cpat.sort
sort -k4,4 $BASEDIR/human/transcripts/allNocode_nr.bed > $BASEDIR/human/transcripts/allNocode_nr.bed.sort

# Join bed file with cpat results
join -1 4 -2 1 -a 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3,2.4,2.5,2.6 $BASEDIR/human/transcripts/allNocode_nr.bed.sort $BASEDIR/human/transcripts/allNocode_nr.cpat.sort > $BASEDIR/human/transcripts/allNocode_nr.cpat.bed



### Create combined denovo and annotated dataset 
# 0.364 is the CPAT threshold for human as suggested by the CPAT paper
awk '$17<0.364' $BASEDIR/human/transcripts/allNocode_nr.cpat.bed | cut -f1-12 | cat - $BASEDIR/nc/Gencode-current_hsa_ncSplice_noCodeOlap.bed |  sort -k1,1 -k2,2n > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_withUTRolaps.bed

# Remove transcripts than have any (i.e. including UTR) overlap >50% with coding genes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_withUTRolaps.bed -b $BASEDIR/coding/Gencode-current_hsa_coding_with_untracked.bed > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_NOUTRolaps_withPseudogenes.bed
# Remove transcripts that overlap >50% pseudogenes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap_NOUTRolaps_withPseudogenes.bed -b $BASEDIR/nc/Gencode-current-pseudogenes.bed > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed

### Get list of transcripts with any real overlap to coding gene - not just coding exons - this is for later annotation
cut -f1-12 $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed | $OVERLAPSELECT inFmt=bed -strand -overlapBases=20  $BASEDIR/coding/Gencode-current_hsa_coding.bed stdin stdout | awk 'BEGIN{OFS="\t"}{print $4,1}' > $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodingGeneOlap.txt




#########################################################
#	     						#
#			MOUSE				#
#							#
#########################################################


# Create reference file containing only coding exons
mkdir -p $BASEDIR/coding
awk '$7!=$8' $BASEDIR/data/Gencode-current_mmu.bed > $BASEDIR/coding/Gencode_mm10_coding.bed
$PSTR getCoding $BASEDIR/coding/Gencode_mm10_coding.bed > $BASEDIR/coding/Gencode_mm10_codingExons.bed

# Annotate all pseudogenes for later filtering
awk '$3~"transcript"{if($14 ~ "IG_C_pseudogene" || $14 ~ "IG_D_pseudogene" || $14 ~ "IG_V_pseudogene" || $14 ~ "polymorphic_pseudogene" || $14 ~ "processed_pseudogene" || $14 ~ "pseudogene" || $14 ~ "transcribed_processed_pseudogene" || $14 ~ "transcribed_unprocessed_pseudogene" || $14 ~ "translated_processed_pseudogene" || $14 ~ "translated_unprocessed_pseudogene" || $14 ~ "TR_J_pseudogene" || $14 ~ "TR_V_pseudogene" || $14 ~ "unitary_pseudogene" || $14 ~ "unprocessed_pseudogene"){print $12}}' $BASEDIR/data/GencodeM4.gtf | perl -pe 's/[\";]//g' | perl -pe 's/\.[0-9]+$//' > $BASEDIR/nc/Gencode-current_mmu-pseudogenes.txt

$PSTR join -x $BASEDIR/nc/Gencode-current_mmu-pseudogenes.txt $BASEDIR/data/Gencode-current_mmu.bed | cut -f1-12 > $BASEDIR/nc/Gencode-current_mmu-pseudogenes.bed


# Process the merged transcriptomes of all mouse tissues (cuffmerge) and:
# 	1) Select those with more than 1 exon and
#	2) Discard transcript in supercontigs ($1!~/_/) and
# 	3) less than 20bp of overlap with a coding exon and
# 	4) less than 50% sense overlap with an annotated transcript (Gencode-current_mmu.bed)
mkdir -p $BASEDIR/mouse/transcripts
awk '$10>1 && $1!~/_/ && $1!="chrM"' $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bed | $OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping -inFmt=bed $BASEDIR/coding/Gencode_mm10_codingExons.bed stdin stdout | $INTERSECT_BED -f 0.5 -s -split -v -a - -b $BASEDIR/data/Gencode-current_mmu.bed| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $BASEDIR/mouse/transcripts/allNocode.bed

# Merge them in a single file and remove redundancy
$PSTR dedup -p $P --exEncomp $BASEDIR/mouse/transcripts/allNocode.bed > $BASEDIR/mouse/transcripts/allNocode_nr.bed


# This blocks produces a reference transcriptome in mouse. It will be used later to blast human ncRNAs against it.
# Extract transcripts from RNA-Seq data that are spliced
awk '$10>1 && $1!~/_/ && $1!="chrM"' $BASEDIR/transcriptomes/cuffmerge/mmu/merged.bed | sort -k1,1 -k2,2nr -k3,3nr -k6,6 -k7,7 -k8,8 -k11,11 -u - > $BASEDIR/mouse/transcripts/deNovoSpliced.bed
# Create file of reference transcriptomes and denovo transcriptomes
cat $BASEDIR/data/Gencode-current_mmu.bed $BASEDIR/mouse/transcripts/deNovoSpliced.bed | sort -k1,1 -k2,2nr -k3,3nr -k6,6 -k7,7 -k8,8 -k11,11 -u - >  $BASEDIR/mouse/transcripts/CombinedRefAnddeNovo.bed


# Calculate coding potential with CPAT
bedtools getfasta -split -name -s -fi $BASEDIR/data/mm10.fa -bed $BASEDIR/mouse/transcripts/allNocode_nr.bed -fo $BASEDIR/mouse/transcripts/allNocode_nr.fa
cpat.py -g $BASEDIR/mouse/transcripts/allNocode_nr.fa -x $DATA/prebuilt_cpat_models/Mouse_Hexamer.tab -d $DATA/prebuilt_cpat_models/Mouse_logitModel.RData -o $BASEDIR/mouse/transcripts/allNocode_nr.cpat

# remove header line and sort by id
tail -n +2 $BASEDIR/mouse/transcripts/allNocode_nr.cpat | sort -k1,1 > $BASEDIR/mouse/transcripts/allNocode_nr.cpat.sort
sort -k4,4 $BASEDIR/mouse/transcripts/allNocode_nr.bed > $BASEDIR/mouse/transcripts/allNocode_nr.bed.sort

# Join bed file with cpat results
join -1 4 -2 1 -a 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3,2.4,2.5,2.6 $BASEDIR/mouse/transcripts/allNocode_nr.bed.sort $BASEDIR/mouse/transcripts/allNocode_nr.cpat.sort > $BASEDIR/mouse/transcripts/allNocode_nr.cpat.bed


# From NovelOrfsTSAnnotated.bedx we are only taking non coding (i.e. no overlap with coding and CPAT < 0.44), while with the awk line below we are excluding Gencode coding genes.
# The awk '$10>1' selects only those with more than one exon
awk '$17<0.44' $BASEDIR/mouse/transcripts/allNocode_nr.cpat.bed | cut -f 1-12 | cat - $BASEDIR/data/Gencode-current_mmu.bed | awk '$7==$8' | awk '$10>1' | $OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping -inFmt=bed $BASEDIR/coding/Gencode_mm10_codingExons.bed stdin stdout | sort -k1,1 -k2,2n > $BASEDIR/nc/mm10CombinedNCTranscripts_withUTRolaps.bed

# Remove transcripts than have any (i.e. including UTR) overlap >50% with coding genes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/mm10CombinedNCTranscripts_withUTRolaps.bed -b $BASEDIR/coding/Gencode_mm10_coding.bed > $BASEDIR/nc/mm10CombinedNCTranscripts_NOUTRolaps_withPseudogenes.bed 
# Remove transcritps that have >50% overlap with Pseudogenes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/mm10CombinedNCTranscripts_NOUTRolaps_withPseudogenes.bed -b $BASEDIR/nc/Gencode-current_mmu-pseudogenes.bed > $BASEDIR/nc/mm10CombinedNCTranscripts.bed

# Filter all Genbank mRNAs to remove transcripts with 
# non canonical splice sites (-N) and merge exons that
# are separated by an intron less than 4bp (-Z).
# Output is GTF (-T)
if [[ ! -f $BASEDIR/data/mm10_allmrna_canonical.bed || $DATA/mm10_allmrna.gtf.gz -nt $BASEDIR/data/mm10_allmrna_canonical.bed ]]; then
	zcat $DATA/mm10_allmrna.gtf.gz > $BASEDIR/data/mm10_allmrna.gtf
	$GFFREAD $BASEDIR/data/mm10_allmrna.gtf -g $BASEDIR/data/mm10.fa -N -T -Z -o $BASEDIR/data/mm10_allmrna_canonical.gtf
	$PSTR gtfToBed -p $P -a geneID $BASEDIR/data/mm10_allmrna_canonical.gtf > $BASEDIR/data/mm10_allmrna_canonical.bed
fi


# Extract spliced all mRNAs with no overlap with coding Exons
awk '$10>1'  $BASEDIR/data/mm10_allmrna_canonical.bed | sort -k1,1 -k2,2n -k3,3n | $OVERLAPSELECT -strand -overlapBases=20 -nonOverlapping -inFmt=bed $BASEDIR/coding/Gencode_mm10_codingExons.bed stdin stdout > $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed

# Calculate coding potential with CPAT
bedtools getfasta -split -name -s -fi $BASEDIR/data/mm10.fa -bed $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed -fo $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.fa
cpat.py -g $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.fa -x $DATA/prebuilt_cpat_models/Mouse_Hexamer.tab -d $DATA/prebuilt_cpat_models/Mouse_logitModel.RData -o $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed.cpat

# Some transcript IDs in mm10_allmrna_spliced_noCodOverlap.bed are duplicate,
# therefore we can't sort the files and join them based on IDs
# The solution here is just to keep the files in the same order as they were 
# given to cpat and then paste the columns together.
tail -n +2 $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed.cpat  > $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed.cpat.nohead

paste $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.bed.cpat.nohead > $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.cpat.bedx

# As a control, after pasting make sure that the two ids are the same
if [[ ! $(awk 'toupper($4)!=toupper($13)' $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.cpat.bedx | wc -l ) -eq 0 ]]; then
	echo "There is an error parsing CPAT output. Exiting";
	exit 1;
fi;

cut -f1-12,18 $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.cpat.bedx > $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.cpat.bed

# Select nc all mRNAs
awk '$13<0.44' $BASEDIR/nc/mm10_allmrna_spliced_noCodOverlap.cpat.bed |  sort -k1,1 -k2,2n | cut -f1-12 > $BASEDIR/nc/mm10_allmrna_spliced_NC_withUTRolaps.bed

# Remove transcripts than have any (i.e. including UTR) overlap >50% with coding genes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/mm10_allmrna_spliced_NC_withUTRolaps.bed -b $BASEDIR/coding/Gencode_mm10_coding.bed > $BASEDIR/nc/mm10_allmrna_spliced_NC_NOUTRolaps_withPseudogenes.bed
# Remove transcritps that have >50% overlap with Pseudogenes
$INTERSECT_BED -v -s -split -f 0.5 -a $BASEDIR/nc/mm10_allmrna_spliced_NC_NOUTRolaps_withPseudogenes.bed -b $BASEDIR/nc/Gencode-current_mmu-pseudogenes.bed > $BASEDIR/nc/mm10_allmrna_spliced_NC.bed


# Add a serial number to the ID so that every ID is unique
# (some IDs are repeated because they map to multiple locations
# in the genome)
awk 'BEGIN{OFS=FS="\t"}{$4=$4"#"NR; print $0}' $BASEDIR/nc/mm10_allmrna_spliced_NC.bed > $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed
cat $BASEDIR/nc/mm10CombinedNCTranscripts.bed $BASEDIR/nc/mm10_allmrna_spliced_NC_uniqID.bed > $BASEDIR/nc/mm10CombinedNCTranscripts+allmrna_spliced_NC_uniqID.bed
md5sum $BASEDIR/nc/Gencode-currentPlusNovel_ncSplice_noCodeOlap.bed $BASEDIR/nc/mm10CombinedNCTranscripts.bed > $BASEDIR/.prepare


