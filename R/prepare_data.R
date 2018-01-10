library(dplyr)
###########################################################
# 		Load pcRNA annotation
###########################################################
BASEDIR="/nfs/research2/enright/tom/projects/new-pcRNAs/analysis"
pc <- read.delim(paste(BASEDIR,"/nc/posConNCwithCodingPartner_andContext_expr.bedx", sep="/"), stringsAsFactors=F, header=F)
colnames(pc) <- c("Chrom","Chromst","Chromend","Name","Score","Strand","CDSSt","CDSEnd","RGB","BlockCount","BlockSize","BlockStarts","Closest Gene Gencode","Close GeneSym","Context","Dist (st-st)","OlapUTR","Novel Genc","linc Overlap","miRNA","Mouse conserved seq","Matched length similarity","Match Length","%transcript matching", "Mouse transcript", "PCpromoterID", "CDpromoterID")


# Set the Score field to the total exonic length
pc$Score<-unlist(lapply(lapply(strsplit(as.character(pc[,"BlockSize"]),","),as.numeric),sum))

# Add information on the Orientation (Sense vs Antisense)
pc[pc$Context=="DS-S" | pc$Context=="US-S" | pc$Context=="OLAP","Orientation"] <- "S"
pc[pc$Context=="Asense" | pc$Context=="BIDIR" | pc$Context=="DS-AS" | pc$Context=="US-AS","Orientation"] <- "AS"
pc$Context <-factor(pc$Context, levels=c("Asense","BIDIR","DS-AS","US-AS","DS-S","OLAP","US-S"))
###########################################################


###########################################################
# 		Load cuffnorm tracking files
###########################################################
hsa_tracking <- read.delim(paste(BASEDIR,"/transcriptomes/cuffnorm/hsa/isoforms.fpkm_tracking", sep="/"), stringsAsFactors=F, header=T)
mmu_tracking <- read.delim(paste(BASEDIR,"/transcriptomes/cuffnorm/mmu/isoforms.fpkm_tracking", sep="/"), stringsAsFactors=F, header=T)
# Extract transcript info
hsa_transcript_info <- hsa_tracking[,1:9]
mmu_transcript_info <- mmu_tracking[,1:9]

# Save the sequencing status
hsa_transcript_status <- hsa_tracking[,grepl("status",colnames(hsa_tracking))]
mmu_transcript_status <- mmu_tracking[,grepl("status",colnames(mmu_tracking))]
row.names(hsa_transcript_status) <- hsa_tracking$tracking_id
row.names(mmu_transcript_status) <- mmu_tracking$tracking_id
###########################################################


###########################################################
#	TRANSCRIPT ID MAPPING FOR HUMAN
###########################################################
# Load the BED file for hsa
hsa_bed <- read.delim(paste(BASEDIR,"transcriptomes/cuffmerge/hsa/merged.bedx",sep="/"),header=F, stringsAsFactors=F)
hsa_bed <- hsa_bed[,c("V4","V14")]
colnames(hsa_bed) <- c("cuff_id", "oId")

# Copy the old ID if it matches ENST
hsa_bed$oId <- gsub("(ENST[0-9]+)\\.[0-9]+$", "\\1", hsa_bed$oId)
hsa_bed$tracking_id <- hsa_bed$cuff_id
hsa_bed[grepl("ENST", hsa_bed$"oId"),"tracking_id"] <- hsa_bed[grepl("ENST", hsa_bed$"oId"),"oId"]

# Initialise a data frame to map pcRNA names to cufflinks IDs
pc2id <- data.frame(pc=pc$Name)

# Merge the tables
pc2id <- merge(pc2id, hsa_bed, by.x=1, by.y=3)
colnames(pc2id) <- c("pc", "pcTCONSid", "oId")
pc <- merge(pc, pc2id[,1:2], by.x="Name", by.y="pc", all.x=T)

# Repeat for coding
cd2id <- data.frame(cd=unique(pc$"Closest Gene Gencode"))

# Read table with Coding transcript id -> coding gene id
coding_tr2gene <- read.delim(paste(BASEDIR,"coding/GencodeV21_coding_trID2geneID.txt",sep="/"),header=F, stringsAsFactors=F)

# Make and ID tracking DF
cd2id <- merge(cd2id, coding_tr2gene, by.x=1, by.y=2)
colnames(cd2id) <- c("cd_gene", "cd_transcript")
cd2id <- merge(cd2id, hsa_bed, by.x="cd_transcript", by.y="tracking_id")
###########################################################


###########################################################
#	TRANSCRIPT ID MAPPING FOR MOUSE
###########################################################
mmu_bed <- read.delim(paste(BASEDIR,"transcriptomes/cuffmerge/mmu/merged.bedx",sep="/"),header=F, stringsAsFactors=F)
mmu_bed <- mmu_bed[,c("V4","V14")]
colnames(mmu_bed) <- c("cuff_id", "oId")

# Copy the old ID if it matches ENST
mmu_bed$oId <- gsub("(ENSMUST[0-9]+)\\.[0-9]+$", "\\1", mmu_bed$oId)
mmu_bed$tracking_id <- mmu_bed$cuff_id
mmu_bed[grepl("ENSMUST", mmu_bed$"oId"),"tracking_id"] <- mmu_bed[grepl("ENSMUST", mmu_bed$"oId"),"oId"]

# One human pc can be associated with multiple mouse transcripts
# Therefore, we split the list of mouse transcripts to obtain a comprehensive list
mouse_pc_names <- strsplit(pc$"Mouse transcript",",")
mouse2human <- data.frame(Mouse=unlist(mouse_pc_names), Human=rep(pc$pcTCONSid, vapply(mouse_pc_names, FUN=length, FUN.VALUE=integer(1))))

# Initialise a data frame to map mouse pcRNA names to cufflinks IDs
mmu_pc2id <- data.frame(pc=unlist(mouse_pc_names))

# Merge the tables
# In this merge mouse pcRNAs not annotated in Gencode or not in the RNA-Seq
# e.g. AK154964 ecc, will be lost.
mmu_pc2id <- unique(merge(mmu_pc2id, mmu_bed, by.x=1, by.y=3))
colnames(mmu_pc2id) <- c("Mouse", "mmu-pcTCONSid", "oId")
mouse2human <- merge(mouse2human, mmu_pc2id[,1:2], by="Mouse",  all.x=T)

# Create an expression matrix for ncRNAs (Transcript ~ Tissue)
hsa_nc_exp <- hsa_tracking[hsa_tracking$tracking_id %in% pc2id$pcTCONSid, grepl("tracking_id|FPKM", colnames(hsa_tracking))]
row.names(hsa_nc_exp) <- hsa_nc_exp$tracking_id
hsa_nc_exp <- hsa_nc_exp[,-1]

# and the same for coding genes
# 16 of the pc[,"Closest Gene Gencode"] are absent from hsa_tracking
# Therefore the will be 16 NA rows in hsa_cd_exp
hsa_cd_exp <- hsa_tracking[hsa_tracking$tracking_id %in% cd2id$cuff_id, grepl("tracking_id|FPKM", colnames(hsa_tracking))]
row.names(hsa_cd_exp) <- hsa_cd_exp$tracking_id
hsa_cd_exp <- hsa_cd_exp[,-1]
hsa_cd_exp <- merge(hsa_cd_exp, cd2id[,c("cuff_id", "cd_gene")], by.x=0, by.y="cuff_id")

hsa_cd_exp <- group_by(hsa_cd_exp, cd_gene) %>% dplyr::summarize(EScepA_FPKM=sum(EScepA_FPKM), EScypA_FPKM=sum(EScypA_FPKM), ESnupA_FPKM=sum(ESnupA_FPKM), GM12878_FPKM=sum(GM12878_FPKM), Hsmm_FPKM=sum(Hsmm_FPKM), K562_FPKM=sum(K562_FPKM), br_FPKM=sum(br_FPKM), cb_FPKM=sum(cb_FPKM), ht_FPKM=sum(ht_FPKM), kd_FPKM=sum(kd_FPKM), lv_FPKM=sum(lv_FPKM), ts_FPKM=sum(ts_FPKM))
hsa_cd_exp <- as.data.frame(hsa_cd_exp)
row.names(hsa_cd_exp) <- hsa_cd_exp$cd_gene
hsa_cd_exp <- hsa_cd_exp[,-1]



# Create an expression matrix for mouse pcRNAs
# Extract expression of coding genes
# 16 of the pc[,"Closest Gene Gencode"] are absent from hsa_tracking
# Therefore the will be 16 NA rows in hsa_cd_exp
mmu_nc_exp <- mmu_tracking[mmu_tracking$tracking_id %in% mmu_pc2id$"mmu-pcTCONSid", grepl("tracking_id|FPKM", colnames(mmu_tracking))]
row.names(mmu_nc_exp) <- mmu_nc_exp$tracking_id
mmu_nc_exp <- mmu_nc_exp[,-1]

# Set extremely low expression values to 0
hsa_nc_exp[hsa_nc_exp<0.001] <- 0
hsa_cd_exp[hsa_cd_exp<0.001] <- 0
mmu_nc_exp[mmu_nc_exp<0.001] <- 0


# Things to save
# pc
# hsa_nc_exp
# hsa_cd_exp
# mmu_nc_exp
# mouse2human

save(pc, file="data/pc.Rdata")
save(hsa_nc_exp, file="data/hsa_nc_exp.Rdata")
write.table(hsa_nc_exp, file="data/hsa_nc_exp.txt", quote=F, sep="\t")
save(hsa_cd_exp, file="data/hsa_cd_exp.Rdata")
save(mmu_nc_exp, file="data/mmu_nc_exp.Rdata")
write.table(mmu_nc_exp, file="data/mmu_nc_exp.txt", quote=F, sep="\t")
save(mouse2human, file="data/mouse2human.Rdata")








