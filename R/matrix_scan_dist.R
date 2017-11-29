# reads RSATools matrix-scan output and outputs a list of motif positions for each fixed length sequence 

matrixScanDist=function(motifFileName)

{

# create output filename

filenamesplit=strsplit(motifFileName,".txt")
outfilename=paste(filenamesplit[length(filenamesplit)],"_distances.txt",sep="")

# read matrix-scan output

GC=read.table(motifFileName,comment.char=";",header=FALSE,sep="\t",quote="",fill=TRUE, skip=33)

# get sequence size/2, first line of new sequence is always -size
distance=-GC$V5[1][1]/2

cat(paste("\nsequence length = ", distance*2,"\n\n",sep=""))

# treat lines with -size apart, they are necessary though because not all sequences will have motifs, and we need all sequences for matrices

GC[(GC$V3=="START_END"),5]=999999;

# extract chrom number, start and end from UCSC sequence identifier

GC$V1=sub("[a-zA-Z0-9_'('')'|]*_range=","",GC$V1,perl = T)


GC$V9=sub("_5'['-=a-zA-Z0-9_|]*","",GC$V1,perl = T)

GC$V10=sub(":[0-9-_|]*","",GC$V9,perl = T)
GC$V11=sub("chr[a-zA-Z0-9]*:","",GC$V9,perl = T)

GC$V11=sub("-[0-9]*","",GC$V11,perl = T)
GC$V12=sub("chr[a-zA-Z0-9-]*:","",GC$V9,perl = T)
GC$V12=sub("[0-9]*-","",GC$V12,perl = T)
GC$V12=as.numeric(GC$V12)

# get midpoint of sequence

GC$V13=GC$V12-distance

# creates identifier with chrX:13456789 format

GC$V14=paste(GC$V10,":",GC$V13, sep="")

# sort by chrom and start, so that sequence IDs will be in alphabetical order

GC=GC[with(GC, order(V10, V13)), ]

# create primary key for each unique sequence

uniqueList=cbind(unique(GC$V14),1:length(unique(GC$V14)))
colnames(uniqueList)=c("coord","ID")

# assign all motifs with sequence key

z=merge(GC,uniqueList,by.x=14,by.y=1)

# distances are relative to 3' end, so need to be reduced by half of sequence length, except lines with 999999 distance

z[!(z$V3=="START_END"),16]=z[!(z$V3=="START_END"),6]+distance
z[(z$V3=="START_END"),16]=999999;

# sort again due to order possibly lost by merge function (sorting generally not necessary though, convenient to debug) and write output

z=z[with(z, order(V10, V13)), ]

write.table(data.frame(z[,15],z[,16],z[,1]),outfilename, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

cat(paste("Successfully wrote motif distances to ",outfilename, "\n\n", sep=""))

}


