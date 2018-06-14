#gffOrderExtend by Pierre Cauchy 2009. Part 1 of ChIP-Seq analysis. Produces sorted sgr file with third column corresponding to tag start (1) or tag end (0). To be used with corresponding gff2wig genome revision file (e.g. gff2wig_mm9.pl)
gffOrderExtend<-function(fileName, extSize){

gff=read.table(fileName, sep="\t")

gff_plus_elongated=gff[gff[,7]=="+",]

gff_plus_elongated[,5]=gff_plus_elongated[,4]+extSize

gff_minus_elongated=gff[gff[,7]=="-",]

gff_minus_elongated[,4]=gff_minus_elongated[,5]-extSize

gff_plus_minus_elongated=rbind(gff_plus_elongated,gff_minus_elongated)

gff_plus_minus_elongated$V8=gff_plus_minus_elongated[,8]=sub('chr','',gff_plus_minus_elongated[,1])

gff_plus_minus_elongated[gff_plus_minus_elongated[,8]=="M",]$V8=23

gff_plus_minus_elongated[gff_plus_minus_elongated[,8]=="X",]$V8=24

gff_plus_minus_elongated[gff_plus_minus_elongated[,8]=="Y",]$V8=25

gff_plus_minus_elongated$V8=as.numeric(gff_plus_minus_elongated$V8)

gff_plus_minus_elongated_ordered=gff_plus_minus_elongated[with(gff_plus_minus_elongated, order(V8, V4)), ]

gff_plus_minus_elongated_ordered$V8=NULL

gff_start=gff_plus_minus_elongated
gff_end=gff_plus_minus_elongated
rm(gff_plus_minus_elongated)
rm(gff)
gff_start$V2=NULL
gff_start$V3=NULL
gff_start$V5=NULL
gff_start$V6=NULL
gff_start$V7=NULL
gff_start$V8=1
gff_end$V2=NULL
gff_end$V3=NULL
gff_end$V4=NULL
gff_end$V6=NULL
gff_end$V7=NULL
gff_end$V8=0
gff_start$start=gff_start$V4
gff_end$start=gff_end$V5
gff_start$V4=NULL
gff_end$V5=NULL
gff_total=rbind(gff_start, gff_end)
rm(gff_start)
rm(gff_end)
gff_total$flag=gff_total$V8
gff_total$V8=NULL
gff_total_ordered=gff_total[with(gff_total, order(V1, start)), ]
rm(gff_total)
gff_total_ordered$start=format(gff_total_ordered$start, scientific=FALSE)
write.table(gff_total_ordered, paste(fileName, "_",extSize, "bp_elongated_ordered.sgr", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


}
