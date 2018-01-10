#!/bin/bash
set -xe -o pipefail
export LC_ALL=C

source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/include.sh


OUT=$BASEDIR/downstream/GOcontrol

mkdir -p $OUT

join -1 4 -2 1 -t $'\t' <(cat $BASEDIR/coding/Gencode-current_hsa_coding_with_untracked.bed | sort -k4,4) <(cut -f1,2 $BASEDIR/data/GencodeV21_info_tr.txt | sort -k1,1) | \
	awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$13,$5,$6,$7,$8,$9,$10,$11,$12,$4-$3}' | \
        sort -k 4,4 -k13,13nr | sort -k4,4 -u | cut -f1-12 | sort -k1,1 -k2,2n -k3,3n > $OUT/Gencode-current_hsa_coding_longest.bed



join -1 1 -2 4 -t $'\t' <( cut -f13 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -u ) <(sort -k4,4 $OUT/Gencode-current_hsa_coding_longest.bed )| awk 'BEGIN{OFS=FS="\t"}{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n -k3,3n > $OUT/pcAssocCodingGenes.bed


# ENSG00000177706 is the first gene of chr7, therefore it has no -b feature
# ClosestBed sets the distance to -1, which seems a reasonable thing to do.
$CLOSEST -N -D ref -t first -id -a $OUT/pcAssocCodingGenes.bed -b $OUT/Gencode-current_hsa_coding_longest.bed | cut -f 4,25 > $OUT/pcAssocCodingGenes_upstream.bed

$CLOSEST -N -D ref -t first -iu -a $OUT/pcAssocCodingGenes.bed -b $OUT/Gencode-current_hsa_coding_longest.bed | cut -f 4,25 > $OUT/pcAssocCodingGenes_downstream.bed

join -1 1 -2 1 -t $'\t' <(sort -k1,1 $OUT/pcAssocCodingGenes_upstream.bed) <(sort -k1,1 $OUT/pcAssocCodingGenes_downstream.bed) | awk 'BEGIN{OFS=FS="\t"}function abs(x){return x<0?-x:x};{ print $1,$2,$3,abs($2)+$3}' > $OUT/pcAssocCodingGenes_islandSize.txt

$CLOSEST -N -D ref -t first -iu -a $OUT/Gencode-current_hsa_coding_longest.bed -b $OUT/Gencode-current_hsa_coding_longest.bed | cut -f 4,25 > $OUT/Gencode-current_hsa_coding_downstream.bed
$CLOSEST -N -D ref -t first -id -a $OUT/Gencode-current_hsa_coding_longest.bed -b $OUT/Gencode-current_hsa_coding_longest.bed | cut -f 4,25 > $OUT/Gencode-current_hsa_coding_upstream.bed

join -1 1 -2 1 -t $'\t' <(sort -k1,1 $OUT/Gencode-current_hsa_coding_upstream.bed) <(sort -k1,1 $OUT/Gencode-current_hsa_coding_downstream.bed) | awk 'BEGIN{OFS=FS="\t"}function abs(x){return x<0?-x:x};{ print $1,$2,$3,abs($2)+$3}' >$OUT/Gencode-current_hsa_coding_islandSize.txt

