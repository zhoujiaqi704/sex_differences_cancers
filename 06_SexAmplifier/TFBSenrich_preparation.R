######tfbs enrichment (LOLA)
######date: 20230328
######author: jqzhou

#library
library(LOLA)
library(data.table)
library(GenomicAlignments)
load("./source/Illumina450_dt.rda")

##TFBS
#load data 
ref<-read.csv("./reference/hg19/TFlist_258_Vorontsov.csv",row.names=1) 
#segment A&B
refA<-ref[which(ref$hg19_cistrome.A>0),]

#modify regionDB
filename <- paste(rownames(refA),".A.bed",sep="")
index= data.table(filename)
index$description="TFBS from Vorontsov et al."

write.table(index,file="./lola_enrichment/hg19/hg19_A/tfbs_encode_A/index.txt",sep="\t",col.names=T,row.names=F,quote=F)


coll<-data.frame(collector="Vorontsov",date="2018-10-23",source="cistrome_A",description="hg19_tfbs")
write.table(coll,file="./lola_enrichment/hg19/hg19_A/tfbs_encode_A/collection.txt",sep="\t",col.names=T,row.names=F,quote=F)

#select *.bed files to new fold
oriDir <- "./lola_enrichment/hg19/tfbs_encode/regions"
tarDir <- "./lola_enrichment/hg19/hg19_A/tfbs_encode_A/regions"

sapply(filename,function(x){file.copy(paste(oriDir,x,sep="/"),tarDir)})


#load regionDB
db=loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19/hg19_A")

##load regiondb
regionDB <- loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19/hg19_A",collection="tfbs_encode_A")

##################################################### +/- 200bp #################################################
##################################################### LIHC ######################################################	
#modify universeSets 
LIHC_dmpTable$cpg=rownames(LIHC_dmpTable)
LIHC_dmpTable= as.data.table(LIHC_dmpTable)
LIHC_bg= Illumina450_dt[match(LIHC_dmpTable$cpg,Illumina450_dt$Name),]
LIHC_bg[,cpgstart_pre:=ifelse(strand=="-",pos-200,pos-199),]
LIHC_bg[,cpgend_pre:=ifelse(strand=="-",pos+200,pos+201),]
LIHC_bg[,cpgID:=LIHC_bg$Name,]

gr_range = with(LIHC_bg,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre))) # +-200bp
gr_cpg = with(LIHC_bg,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))

overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]

LIHC_bg[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
LIHC_bg[,cpgend:=end(gr_range[overlap_red$subjectHit])]
LIHC_bg[,NsubjectHits:=overlap_red$NsubjectHits]

LIHC_bg_sub=LIHC_bg[,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
setnames(LIHC_bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))

LIHC_bg_gr=unique(with(LIHC_bg, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID)))

save(LIHC_bg,LIHC_bg_gr,file="./lola_enrichment/LOLA_LIHC_bg_gr_200bp.rda")



