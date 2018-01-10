SHELL := /bin/bash
DATE=`date +'%y-%m-%d'`

include scripts/include.sh


all: save report

download: $(BASEDIR)/.download
$(BASEDIR)/.download:
	@echo "Downloading data"
	scripts/download.sh
	scripts/download_fastq.sh

map: $(BASEDIR)/.map scripts/map.sh
$(BASEDIR)/.map: $(BASEDIR)/.download
	@echo "Mapping RNA-Seq data"
	scripts/map.sh

coverage: $(BASEDIR)/.coverage
$(BASEDIR)/.coverage: $(BASEDIR)/.map scripts/coverage.sh
	@echo "Calculating coverage (producing bw files)"
	scripts/coverage.sh

assemble: $(BASEDIR)/.assemble
$(BASEDIR)/.assemble: $(BASEDIR)/.map scripts/assemble.sh
	@echo "Assembling transcripts"
	scripts/assemble.sh

merge: $(BASEDIR)/.merge
$(BASEDIR)/.merge: $(BASEDIR)/.assemble scripts/merge.sh
	@echo "Merging assembled transcripts with cuffcompare"
	scripts/merge.sh

quantitate: $(BASEDIR)/.quant
$(BASEDIR)/.quant: $(BASEDIR)/.merge scripts/quant.sh
	@echo "Quantitating transcripts with cuffquant"
	scripts/quant.sh

cuffnorm: $(BASEDIR)/.cuffnorm 
$(BASEDIR)/.cuffnorm: $(BASEDIR)/.quant scripts/cuffnorm.sh
	@echo "Quantitating transcripts with cuffnorm"
	scripts/cuffnorm.sh

prepare: $(BASEDIR)/.prepare
$(BASEDIR)/.prepare: $(BASEDIR)/.merge scripts/prepare.sh
	@echo "Preparing annotations"
	scripts/prepare.sh

identify_blast: $(BASEDIR)/.identify_blast
$(BASEDIR)/.identify_blast: $(BASEDIR)/.prepare scripts/identify_blast.sh
	@echo "Identifying pcRNAs: BLAST step"
	scripts/identify_blast.sh

identify_find: $(BASEDIR)/.identify_find 
$(BASEDIR)/.identify_find: $(BASEDIR)/.identify_blast scripts/identify_find.sh
	@echo "Identifying pcRNAs: FIND step"
	scripts/identify_find.sh

identify: identify_blast identify_find

annotate: $(BASEDIR)/.annotate
$(BASEDIR)/.annotate: $(BASEDIR)/.identify_find scripts/annotate.sh
	@echo "Annotating pcRNAs"
	scripts/annotate.sh

hic: $(BASEDIR)/.hic 
$(BASEDIR)/.hic: $(BASEDIR)/.annotate scripts/other_analysis/hic.sh
	@echo "Analysing HiC data"
	scripts/other_analysis/hic.sh

tad: $(BASEDIR)/.tad
$(BASEDIR)/.tad: $(BASEDIR)/.hic scripts/other_analysis/hic_TADs.sh
	@echo "Analysing HiC-TAD data"
	scripts/other_analysis/hic_TADs.sh

histones: $(BASEDIR)/.histones
$(BASEDIR)/.histones: $(BASEDIR)/.annotate scripts/other_analysis/histones.sh R/report/report.html scripts/other_analysis/ $(BASEDIR)/.save
	 @echo "Making histone modification heatmaps"
	 scripts/other_analysis/histones.sh && touch $(BASEDIR)/.histones

save: $(BASEDIR)/.save
$(BASEDIR)/.save: $(BASEDIR)/.annotate $(BASEDIR)/.hic scripts/save_results.sh R/report/report.html scripts/other_analysis/gen_stat.sh
	@echo "Saving results"
	scripts/save_results.sh && scripts/other_analysis/gen_stat.sh

report: R/report/report.html R/Makefile
R/report/report.html: $(BASEDIR)/.annotate $(BASEDIR)/.cuffnorm $(BASEDIR)/.hic $(BASEDIR)/.tad
	@echo "Running R analysis"
	cd R && make

.PHONY: website
website: $(BASEDIR)/.save R/report/report.html R/report/nanostring.html $(BASEDIR)/.histones
	cp -R $(RESULTS)/* website/data
	mkdir -p website/data/GO && cp -R R/report/data/GO/* website/data/GO
	cp -R R/report/*-fig website
	cp R/report/report.html website
	cp R/report/nanostring.html website
	mkdir -p website/data/histones/posCat && cp $(BASEDIR)/downstream/histones/annotationByPosition/*.pdf website/data/histones/posCat
	mkdir -p website/data/histones/singleCat && cp $(BASEDIR)/downstream/histones/annotationAllPc/*.pdf website/data/histones/singleCat
	mkdir -p website/data/histones/ && cp $(BASEDIR)/downstream/histones/Coverage_H3K27me3/results.pdf website/data/histones/Coverage_H3K27me3.pdf
	mkdir -p website/kd-arrays && cp -R other_analysis/microarray-KD/results/* website/kd-arrays
	cp $(BASEDIR)/downstream/hic/pcCoverageByLoopsAllLines/pcCoverageByLoops_profile.p* website/hic-fig
	cp $(BASEDIR)/downstream/hic/pcCoverageByLoopsAllLines/pcCoverageByLoops_pcByPos_profile.p* website/hic-fig
	cp $(BASEDIR)/downstream/hic/pcCoverageByLoopsAllLines/pcCoverageByLoop_promoters_profile.p* website/hic-fig
	cp $(BASEDIR)/downstream/ctcf/pcCoverageByCtcfAllLines/CTCF_coverage_profile.p* website/hic-fig
	cp $(BASEDIR)/downstream/hic-tads/pcCoveragebyTADends/pcCoverageByTADs_profile.p* website/hic-fig

.PHONY: clean_website
clean_website:
	rm -rf website/data/*
	rm -rf website/*-fig
	rm -rf website/report.html

