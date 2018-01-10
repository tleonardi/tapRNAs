# Main pipeline for pcRNAs identification

### download.sh download_fastq.sh
Download annotations and RNA-Seq fastq files

### map.sh
Map the RNA-Seq data to the reference genomes

### coverage.sh
Produce bigWig coverage files for each RNA-Seq dataset

### assemble.sh
Run cufflinks to assemble each RNA-Seq dataset

### merge.sh
Run cuffmerge to produce reference transcriptomes for mouse and human

### quant.sh
Run cuffquant to quantitate each RNA-Seq dataset againts the assembled transcriptome

### cuffnorm.sh
Generate expression matrixes for human and mouse

### prepare.sh
Prepare all annotations for the identification of pcRNAs (make references for coding and non-coding, filter, calculate coding potential, etc...)

### identify_blast.sh
Blast human promoters to mouse to identify candidate pcRNAs

### identify_find.sh
Identify pcRNAs

### annotate.sh
Produce a detailed annotation of pcRNAs

### save_results.sh
Save useful files in the 'results' folder

### make_hub.sh
Make files for the genome browser hub

### other_analysis
Downstream analysis
