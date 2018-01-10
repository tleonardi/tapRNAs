# Pipeline for the identification of pcRNAs and tapRNAs

## Abstract
The mammalian genome is transcribed into large numbers of long noncoding RNAs (lncRNAs), but the definition of functional lncRNA groups has proven difficult, partly due to their low sequence conservation and lack of identified shared properties. Here we consider positional conservation across mammalian genomes as an indicator of functional commonality. We identify 665 conserved lncRNA promoters in mouse and human genomes that are preserved in genomic position relative to orthologous coding genes. The identified positionally conserved lncRNA genes are primarily associated with developmental transcription factor loci with which they are co-expressed in a tissue-specific manner. Strikingly, over half of all positionally conserved RNAs in this set are linked to distinct chromatin organization structures, overlapping the binding sites for the CTCF chromatin organizer and located at chromatin loop anchor points and borders of topologically associating domains (TADs). These topological anchor point (tap)RNAs possess conserved sequence domains that are enriched in potential recognition motifs for Zinc Finger proteins. Characterization of these non-coding RNAs and their associated coding genes shows that they are functionally connected: they regulate each other ′s expression and influence the metastatic phenotype of cancer cells in vitro in a similar fashion. Thus, interrogation of positionally conserved lncRNAs identifies a new subset of tapRNAs with shared functional properties. These results provide a large dataset of lncRNAs that conform to the ″extended gene″ model, in which conserved developmental genes are genomically and functionally linked to regulatory lncRNA loci across mammalian evolution.

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
* [Pfam Scan](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/)
* [EMBOSS suite](http://emboss.sourceforge.net/index.html)

## Citation
_Genomic positional conservation identifies topological anchor point (tap)RNAs linked to developmental loci_
*Paulo P Amaral, Tommaso Leonardi, Namshik Han, Emmanuelle Vire, Dennis K Gascoigne, Raul Arias-Carrasco, Magdalena Buscher, Anda Zhang, Stefano Pluchino, Vinicius Maracaja-Coutinho, Helder I Nakaya, Martin Hemberg, Ramin Shiekhattar, Anton J Enright, Tony Kouzarides*
bioRxiv 051052; doi: https://doi.org/10.1101/051052

