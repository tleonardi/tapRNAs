#!/bin/bash
set -xe -o pipefail
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../include.sh

RESULTS=$BASEDIR/../results

(
echo -e "Total pcRNAs (redundant)\t" $(cut -f4 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)
echo -e "Unique promoters\t" $(cut -f26 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -u | wc -l)
echo -e "Unique protein coding genes\t" $(cut -f13 $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | sort -u | wc -l)
echo -e "pcRNAs overlapping UTRs\t" $(awk '$17>0' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)
echo -e "pcRNAs annotated in Gencode\t" $(awk '$4~/ENST/'  $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l) "(" $(awk '$4~/ENST/'  $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | cut -f26 | sort -u | wc -l) "promoters )"
echo -e "pcRNAs not annotated in Gencode (novel)\t" $(awk '$4~/TCONS/'  $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)  "(" $(awk '$4~/TCONS/'  $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | cut -f26 | sort -u | wc -l) "promoters )"
echo -e "pcRNAs not overlapping Gencode transcripts (<20%)\t" $(awk '$18>0' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)
echo -e "pcRNAs overlapping Gencode transcripts (>20%)\t" $(awk '$18==0' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)
echo -e "pcRNAs overlapping Gencode lincRNAs\t" $(awk '$19>0' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l)
echo -e "pcRNAs overlapping miRNAs\t" $(awk '$20!="_"' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | wc -l) "(" $(awk '$20!="_"' $BASEDIR/nc/posConNCwithCodingPartner_andContext_expr.bedx | cut -f26 | sort -u | wc -l)  "promoters )"
echo -e "Total number of human mapped reads\t" $((for i in $BASEDIR/transcriptomes/tophat/hsa/*; do grep Mapped $i/align_summary.txt; done) | awk '{TOT+=$3}END{print TOT}')
echo -e "Total number of mouse mapped reads\t" $((for i in $BASEDIR/transcriptomes/tophat/mmu/*; do grep Mapped $i/align_summary.txt; done) | awk '{TOT+=$3}END{print TOT}')
echo -e "pcRNAs that overlap a loop end point\t" $(cut -f1 $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | sed 's/,/\n/g' | wc -l)
echo -e "pcRNA promoters that overlap a loop end point\t" $(cat $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | wc -l)
echo -e "pcRNAs that overlap a TAD boundary (+/-10kb)\t" $(cut -f1 $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | sed 's/,/\n/g' | wc -l)
echo -e "pcRNA promoters that overlap a TAD boundary (+/-10kb)\t" $(cat $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | wc -l)
echo -e "pcRNAs that overlap a loop end point and DON'T overlap a TAD boundary\t" $(comm -23 <(cut -f1 $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | sed 's/,/\n/g' | sort) <(cut -f1 $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | sed 's/,/\n/g' | sort) | wc -l)
echo -e "pcRNAs that overlap a TAD boundary and DON'T overlap a loop end point\t" $(comm -13 <(cut -f1 $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | sed 's/,/\n/g' | sort) <(cut -f1 $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | sed 's/,/\n/g' | sort) | wc -l)
echo -e "pcRNAs that overlap a loop end point AND a TAD boundary\t" $(comm -12 <(cut -f1 $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | sed 's/,/\n/g' | sort) <(cut -f1 $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | sed 's/,/\n/g' | sort) | wc -l)
echo -e "pcRNAs that overlap a loop end point OR a TAD boundary\t" $(cat <(cut -f1 $BASEDIR/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt | sed 's/,/\n/g') <(cut -f1 $BASEDIR/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt | sed 's/,/\n/g') | sort -u | wc -l)
)>$RESULTS/gen_stat.txt

