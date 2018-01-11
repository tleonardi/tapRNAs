[![DOI](https://zenodo.org/badge/116976101.svg)](https://zenodo.org/badge/latestdoi/116976101)

# Pipeline for the identification of pcRNAs and tapRNAs
_Paulo P. Amaral*, Tommaso Leonardi*, Namshik Han*, Emmanuelle Viré, Dennis K. Gascoigne, Raúl Arias-Carrasco, Magdalena Büscher, Luca Pandolfini, Anda Zhang, Stefano Pluchino, Vinicius Maracaja-Coutinho, Helder I. Nakaya, Martin Hemberg, Ramin Shiekhattar, Anton J. Enright, Tony Kouzarides_

## Abstract
### Background
The mammalian genome is transcribed into large numbers of long noncoding RNAs (lncRNAs), but the definition of functional lncRNA groups has proven difficult, partly due to their low sequence conservation and lack of identified shared properties. Here we consider promoter conservation and positional conservation as indicators of functional commonality.

### Results
We identify 665 conserved lncRNA promoters in mouse and human that are preserved in genomic position relative to orthologous coding genes. These positionally conserved lncRNA genes are primarily associated with developmental transcription factor loci with which they are coexpressed in a tissue-specific manner. Over half of positionally conserved RNAs in this set are linked to chromatin organization structures, overlapping binding sites for the CTCF chromatin organizer and located at chromatin loop anchor points and borders of topologically associating domains (TADs). We define these RNAs as topological anchor point RNAs (tapRNAs). Characterization of these noncoding RNAs and their associated coding genes shows that they are functionally connected: they regulate each other’s expression and influence the metastatic phenotype of cancer cells in vitro in a similar fashion. Furthermore, we find that tapRNAs contain conserved sequence domains that are enriched in motifs for zinc finger domain-containing RNA-binding proteins and transcription factors, whose binding sites are found mutated in cancers.

### Conclusions
This work leverages positional conservation to identify lncRNAs with potential importance in genome organization, development and disease. The evidence that many developmental transcription factors are physically and functionally connected to lncRNAs represents an exciting stepping-stone to further our understanding of genome regulation.

## Disclaimer
This Github repository contains the code of the bioinformatics analysis for our [paper on tapRNAs](https://www.biorxiv.org/content/early/2016/05/04/051052).
We published this repository to allow the scientific community to verify and reproduce our results, but the code is not to be intended as a standalone, general purpose tool. 
Although the majority of the scripts should run without modifications on any GNU/Linux distribution, some of them are tied to the system where they were developed (e.g. make use of IBM LSF queing system).

## Depedencies

* [BEDTools](https://github.com/arq5x/bedtools2)
* [Pinstripe](http://pinstripe.matticklab.com/)
* [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml)
* [Bowtie](bowtie-bio.sourceforge.net)
* [Cufflinks](cole-trapnell-lab.github.io/cufflinks)
* [Samtools](samtools.sourceforge.net/)
* [UCSC utilities](http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)
* [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Pfam Scan](http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/)
* [EMBOSS suite](http://emboss.sourceforge.net/index.html)
* [deepTools](https://github.com/deeptools/deepTools)
* [ImageMagick](http://www.imagemagick.org)
* [HMMER](http://hmmer.org/)
* [R](https://www.r-project.org/) 
* R packages: boot, data.table, dplyr, gdata, ggbiplot, ggplot2, gplots, gridExtra, MatchIt, org.Hs.eg.db, plyr, RColorBrewer, reshape, reshape2, scales, sfsmisc, topGO, VennDiagram 

## Citation
**Genomic positional conservation identifies topological anchor point (tap)RNAs linked to developmental loci**

_Paulo P. Amaral*, Tommaso Leonardi*, Namshik Han*, Emmanuelle Viré, Dennis K. Gascoigne, Raúl Arias-Carrasco, Magdalena Büscher, Luca Pandolfini, Anda Zhang, Stefano Pluchino, Vinicius Maracaja-Coutinho, Helder I. Nakaya, Martin Hemberg, Ramin Shiekhattar, Anton J. Enright, Tony Kouzarides_

bioRxiv 051052; doi: https://doi.org/10.1101/051052

