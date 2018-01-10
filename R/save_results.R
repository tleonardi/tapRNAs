library(dplyr)
load("data/pc_human_expAnalysis.Rdata")

pc_clean <- mutate(pc, CodingID=`Closest Gene Gencode`, CodingSymbol=`Close GeneSym`, Dist=`Dist (st-st)`) %>% select(Name, CodingID, CodingSymbol, Context, Dist, OlapUTR, `Novel Genc`, `linc Overlap`, miRNA, PCpromoterID, pcTCONSid, `Mouse transcript`, tisSpecPC, maxTisPC, maxFpkmPC, tisSpecCD, maxTisCD , maxFpkmCD , pc_cd_rho, pc_cd_r) 


loops <- read.delim("/nfs/gns/homes/tl344/tom/projects/new-pcRNAs/analysis/downstream/hic/loopOverlapAllLinesPromoter/pc_with_loops_in_prom.txt", header=F, stringsAsFactors=F)
pl <- data.frame(Name=unique(unlist(strsplit(loops[,1], ","))), PromoterLoop="Yes", stringsAsFactors=F)

pc_clean <- merge(pc_clean, pl, by="Name", all.x=T)

pc_clean[is.na(pc_clean$PromoterLoop), "PromoterLoop"] <- "No"

tads <- read.delim("/nfs/gns/homes/tl344/tom/projects/new-pcRNAs/analysis/downstream/hic-tads/TADOverlapAllLinesPromoter/pc_with_TADs_in_prom_10000kb.txt", header=F, stringsAsFactors=F)
pt <- data.frame(Name=unique(unlist(strsplit(tads[,1], ","))), PromoterTADBoundary="Yes", stringsAsFactors=F)

pc_clean <- merge(pc_clean, pt, by="Name", all.x=T)
pc_clean[is.na(pc_clean$PromoterTADBoundary), "PromoterTADBoundary"] <- "No"

write.table(pc_clean, file="pc_annotation.txt", quote=F, sep="\t", row.names=F)
