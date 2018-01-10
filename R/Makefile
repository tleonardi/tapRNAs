all: report/report.html report/pc_annotation.txt report/nanostring.html

BASEDIR=/path/to/repo/analysis/


prepare_out=report/data/pc.Rdata report/data/hsa_nc_exp.Rdata report/data/hsa_nc_exp.txt report/data/hsa_cd_exp.Rdata report/data/mmu_nc_exp.Rdata report/data/mmu_nc_exp.txt report/data/mouse2human.Rdata .prepare

prepare_input=$(BASEDIR)/nc/posConNCwithCodingPartner_andContext_expr.bedx $(BASEDIR)/transcriptomes/cuffnorm/hsa/isoforms.fpkm_tracking $(BASEDIR)/transcriptomes/cuffnorm/mmu/isoforms.fpkm_tracking $(BASEDIR)/transcriptomes/cuffmerge/hsa/merged.bedx $(BASEDIR)/transcriptomes/cuffmerge/mmu/merged.bedx

HIC=/nfs/gns/homes/tl344/tom/projects/new-pcRNAs/analysis/downstream/hic/
HICTADS= /nfs/gns/homes/tl344/tom/projects/new-pcRNAs/analysis/downstream/hic-tads/
hic_input=$(HIC)/pcRNA_loop_endPoints_hmmGm12878_table.txt $(HIC)/pcRNA_loop_endPoints_hmmHmec_table.txt $(HIC)/pcRNA_loop_endPoints_hmmK562_table.txt $(HIC)/pcRNA_loop_endPoints_hmmH1hesc_table.txt $(HIC)/pcRNA_loop_endPoints_hmmHsmm_table.txt $(HIC)/pcRNA_loop_endPoints_hmmNhek_table.txt $(HIC)/pcRNA_loop_endPoints_hmmHepg2_table.txt $(HIC)/pcRNA_loop_endPoints_hmmHuvec_table.txt $(HIC)/pcRNA_loop_endPoints_hmmNhlf_table.txt $(HIC)/genome_loop_endPoints_hmmGm12878_table.txt $(HIC)/genome_loop_endPoints_hmmHuvec_table.txt $(HIC)/genome_loop_endPoints_hmmH1hesc_table.txt $(HIC)/genome_loop_endPoints_hmmK562_table.txt $(HIC)/genome_loop_endPoints_hmmHepg2_table.txt $(HIC)/genome_loop_endPoints_hmmNhek_table.txt $(HIC)/genome_loop_endPoints_hmmHmec_table.txt $(HIC)/genome_loop_endPoints_hmmNhlf_table.txt $(HICTADS)/TADOverlapAllLinesPromoter/TADs_contact_summary_10000kb.txt $(HICTADS)/TADOverlapAllLinesPromoter/TADs_contact_summary_all.txt


$(prepare_out): prepare_data.R functions.R $(prepare_input)
	@echo "Preparing data"
	mkdir -p report/data
	/homes/tl344/tom/bin/R312/bin/Rscript -e "setwd('report'); source('../prepare_data.R')" && touch .prepare

report/general_stat.md: general_stat.Rmd report/data/pc.Rdata
	@echo "Running general stat"
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../general_stat.Rmd')"

report/GO.md: GO.Rmd report/data/pc.Rdata
	@echo "Doing GO enrichment"
	mkdir -p report/data/GO && /homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../GO.Rmd')"

report/human_expAnalysis.md report/data/pc_human_expAnalysis.Rdata: human_expAnalysis.Rmd report/data/pc.Rdata report/data/hsa_nc_exp.Rdata report/data/hsa_cd_exp.Rdata
	@echo "Calculating human expression analysis"
	mkdir -p report/data/GO/byTissue && /homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../human_expAnalysis.Rmd')"

report/mouse_expAnalysis.md: mouse_expAnalysis.Rmd report/data/pc_human_expAnalysis.Rdata report/data/hsa_nc_exp.Rdata report/data/mmu_nc_exp.Rdata report/data/mouse2human.Rdata
	@echo "Calculating mouse expression analysis"
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../mouse_expAnalysis.Rmd')"

report/hic.md: hic.Rmd $(hic_input) $(BASEDIR)/downstream/hic/loopOverlapAllLines/loops_contact_summary.txt $(BASEDIR)/downstream/hic/loopOverlapAllLinesPromoter/loops_contact_summary.txt 
	@echo "Plotting results of HiC analysis"
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../hic.Rmd')"

report/hic_go.md: hic_go.Rmd report/hic.md
	@echo "HiC GO enrichment"
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report');knit('../hic_go.Rmd')"

report/cover.md: cover.Rmd report/general_stat.md
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report'); knit('../cover.Rmd')"

report/nanostring.md: nanostring.Rmd
	/homes/tl344/tom/bin/R312/bin/Rscript -e "library(knitr); setwd('report'); knit('../nanostring.Rmd')"

report/nanostring.html: report/nanostring.md template.htm
	cd report && pandoc nanostring.md --toc -t html5 -s --template ../template.htm -o nanostring.html

report/report.html report/report.md: report/cover.md report/general_stat.md report/GO.md report/human_expAnalysis.md report/mouse_expAnalysis.md report/hic.md template.htm
	@echo "Converting md files to html"
	cd report && pandoc cover.md general_stat.md GO.md human_expAnalysis.md mouse_expAnalysis.md hic.md hic_go.md --toc -t html5 -s --template ../template.htm -o report.html
	cd report && pandoc cover.md general_stat.md GO.md human_expAnalysis.md mouse_expAnalysis.md hic.md hic_go.md -t markdown_github --toc -s -o report.md

report/pc_annotation.txt: report/report.html save_results.R
	@echo Saving results
	/homes/tl344/tom/bin/R312/bin/Rscript -e "setwd('report'); source('../save_results.R')"

