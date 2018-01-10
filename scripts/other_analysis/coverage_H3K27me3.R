library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library("topGO")
library("org.Hs.eg.db")
library(scales)

BASEDIR="/path/to/basedir"
load(paste0(BASEDIR,"/../R/report/data/pc.Rdata"))

Gm12878 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/Gm12878H3k27me3.map"), header=F) %>% mutate(Gm12878_K27=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
H1hesc <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/H1hescH3k27me3.map"), header=F) %>% mutate(H1hesc_K27=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
Hsmm <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/HsmmH3k27me3.map"), header=F) %>% mutate(Hsmm_K27=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
K562 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/K562H3k27me3.map"), header=F) %>% mutate(K562_K27=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)

Gm12878k4 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/Gm12878H3k4me3.map"), header=F) %>% mutate(Gm12878_K4=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
H1hesck4 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/H1hescH3k4me3.map"), header=F)  %>% mutate(H1hesc_K4=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
Hsmmk4 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/HsmmH3k4me3.map"), header=F) %>% mutate(Hsmm_K4=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)
K562k4 <- read.table(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/K562H3k4me3.map"), header=F) %>% mutate(K562_K4=V7) %>% dplyr::select(-V7,-V1,-V2,-V3,-V5,-V6)

pdf(paste0(BASEDIR,"/analysis/downstream/histones/Coverage_H3K27me3/results.pdf"), width=10)
df <- merge(Gm12878, H1hesc, by="V4") %>%  merge(., Hsmm, by="V4")  %>% merge(., K562, by="V4")  %>% merge(., Gm12878k4, by="V4")  %>% merge(., H1hesck4, by="V4")  %>% merge(., Hsmmk4, by="V4")  %>% merge(.,K562k4, by="V4")


df_m <- melt(df) %>% mutate(Line=gsub("(.+)_(.+)", "\\1", variable), Mod=gsub("(.+)_(.+)", "\\2", variable))
df_m_c <- dcast(df_m, V4+Line~Mod, value="value")
print(ggplot(df_m, aes(y=value, x=Mod, fill=Line)) + geom_boxplot() + ggtitle("H3k27me3 promoter coverage") + xlab("Cell type") + ylab("Mean promoter coverage"))

print(ggplot(df_m_c, aes(x=K4, y=K27)) + geom_point() + facet_wrap(~Line) + scale_x_log10() + scale_y_log10())


# K-Means clustering
esc <- df_m_c[df_m_c$Line=="H1hesc",]
esc$Cluster <- kmeans(log10(esc[,c("K4","K27")]+0.01),3)$cluster

# Hierarchical clustering
hc <- hclust(dist(log10(esc[,c("K4","K27")]+0.01)))
hc_cut <- cutree(hc, 5)
esc$Cluster2 <- hc_cut
# Rename the clusters to something more sensible
esc[hc_cut==1, "Cluster2"] <- 2
esc[hc_cut==2, "Cluster2"] <- 3
esc[hc_cut==3, "Cluster2"] <- 1

esc <- dplyr::select(pc, Name, Context) %>% merge(esc, ., by.x="V4", by.y="Name")

print(ggplot(esc, aes(x=K4, y=K27, colour=as.factor(Cluster2))) + geom_point() + scale_x_log10() + scale_y_log10())
print(ggplot(esc, aes(x=K4, y=K27, colour=as.factor(Cluster2))) + geom_point() + scale_x_log10() + scale_y_log10() + facet_wrap(~Context))

df_m_c_annot <- merge(df_m_c, esc[,c("V4","Cluster2")], by="V4")
print(ggplot(df_m_c_annot, aes(x=K4, y=K27, colour=as.factor(Cluster2))) + geom_point() + facet_wrap(~Line) + scale_x_log10() + scale_y_log10())

# Load and merge expression data
load(paste0("/../R/report/data/hsa_nc_exp_mean.Rdata"))
cluster2name <- dplyr::select(pc, Name, pcTCONSid) %>% merge(., esc, by.x="Name", by.y="V4") %>% dplyr::select(pcTCONSid, Cluster2)
hsa_nc_exp_mean_filt <- hsa_nc_exp_mean[apply(hsa_nc_exp_mean,1,max)>0,]
hsa_nc_exp_mean_m <- melt(as.matrix(hsa_nc_exp_mean_filt)) %>% merge(., cluster2name, by.x="Var1", by.y="pcTCONSid")

hsa_nc_exp_mean_m$Expressed <- "No"
hsa_nc_exp_mean_m[hsa_nc_exp_mean_m$value>0.1, "Expressed"] <- "Yes"
hsa_nc_exp_mean_m$Expressed <- as.factor(hsa_nc_exp_mean_m$Expressed)

print(ggplot(filter(hsa_nc_exp_mean_m, Cluster2!=5), aes(x=as.factor(as.character(Cluster2)), y=log10(value), colour=as.character(Cluster2))) + geom_point(alpha=0.2, position="jitter") + facet_wrap(~Var2) + geom_boxplot(alpha=0, colour="black", notch=T))
print(ggplot(hsa_nc_exp_mean_m, aes(x=as.character(Cluster2), fill=as.character(Expressed))) + geom_bar(position="dodge") + facet_wrap(~Var2))
print(ggplot(hsa_nc_exp_mean_m, aes(x=as.character(Cluster2), fill=as.character(Expressed))) + geom_bar(position="fill") + facet_wrap(~Var2))


GOres <- list()
# Get list of all genes with synteny in mouse
bgAnnot<- read.delim(paste(BASEDIR,"data/hg38ToMm10_orthologs_mart.txt", sep="/"), stringsAsFactors=F, header=F)
bgGenes <- unique(bgAnnot[grepl("ENSMUSG", bgAnnot$V3),"V1"])

for( i in 1:max(esc$Cluster2)){
	sel_pc <- filter(esc, Cluster2==i) %>% dplyr::select(V4)
	
	sel_names <- unique(filter(pc, pc$Name %in% sel_pc$V4) %>% dplyr::select(`Close GeneSym`))
	sel_names <- sel_names[,1]
	sel_ens <- unique(filter(pc, pc$Name %in% sel_pc$V4) %>% dplyr::select(`Closest Gene Gencode`))
	sel_ens <- sel_ens[,1]
	
	# Get the list of genes associated with pcRNAs
	selectedGenes <- sel_ens
	if(length(selectedGenes)>1){	
		# Make the geneList for topGO
		geneList <- factor(as.integer( bgGenes %in% selectedGenes))
		names(geneList) <- bgGenes
		
		# Create topGO object
		data <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, ID = "ensembl", mapping="org.Hs.eg.db", nodeSize=3)
		
		# Run test and correct p-vals
		resultFis.def <- runTest(data, statistic = "fisher")
		resultFis.clas<- runTest(data, algorithm="classic", statistic = "fisher")
		
		score(resultFis.def) <- p.adjust(score(resultFis.def), method="BH")
		score(resultFis.clas) <- p.adjust(score(resultFis.clas), method="BH")
		
		allRes <- GenTable(data, Fisherdef = resultFis.def, Fisherclas=resultFis.clas, orderBy="Fisherclas", topNodes=20)
		allRes$Fisherdef <- as.numeric(allRes$Fisherdef)
		allRes$Fisherclas<- as.numeric(allRes$Fisherclas)
		allRes$Term <- paste(allRes$Term,allRes$GO.ID)
		allRes$Term <- factor(allRes$Term, levels=allRes$Term[order(allRes$Significant/allRes$Annotated,decreasing=F)])
		
		allGO = genesInTerm(data)
		my_annotation= lapply(allGO,function(x) x[x %in% selectedGenes] )
		sig_tab <- allRes
		sig_tab$genes <- apply(allRes,1,function(x) paste(unique(unlist(pc[pc[,"Closest Gene Gencode"] %in% my_annotation[[x[1]]], "Close GeneSym"]), collapse=",")))
		GOres[[i]] <- sig_tab
		print(ggplot(sig_tab[sig_tab$Fisherclas<0.1,], aes(y=Significant/Annotated,x=Term,colour=Fisherclas, size=Significant)) + geom_point(aes(order=plyr::desc(Fisherclas))) + scale_colour_gradient(name="p-value", low="red",high="yellow")   + coord_flip() + scale_size(range = c(3, 10)) + ggtitle(paste("Cluster ",i, sep="")))
	}
}
dev.off()
