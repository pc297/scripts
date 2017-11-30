############# foldChangeDiffPeak.R  by Pierre Cauchy, 2016
# filters peaks based on gaussian fit (the premise being that background peaks and true peaks will occur in different distributions)
# retrieves counts pairwise from wig files based on union of two peak files, and rank by fold change, creating equal size cdt files
# usage: main(file1,file2,wigFile1,wigFile2,peakDist=400,profileSize=2000,countsRange=400,fcCutoff=1,peakFilteringMethod="none",cutoffMethod="foldChange",blackList="",simpleRepeats="",genome="hg38",histSize=10,outdir="")
# requires mixtools, bedtools and Homer (both externally)

library(mixtools)

# gaussian.fit.remove.bg.peaks
# remove weak, background peaks called as peaks by MACS based on Gaussian fit model mix
gaussian.fit.remove.bg.peaks=function(file)
{
	
	summits=read.table(file,head=F,sep="\t")
	summitsmdl=normalmixEM(log(summits$V5,2),k=2)
	cutoff=summitsmdl$mu[1]+2*summitsmdl$sigma[1]
	sig.summits=summits[summits$V5>2^cutoff,]
	write.table(sig.summits,paste(strsplit(file,".bed")[[1]],"_significant.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
	# return absolute file name
	return(paste(strsplit(file,".bed")[[1]],"_significant.bed",sep=""))
	#return(paste(unlist(strsplit(strsplit(file,".bed")[[1]],"/"))[[length(unlist(strsplit(strsplit(file,".bed")[[1]],"/")))]],"_significant.bed",sep=""))
}

# median.remove.bg.peaks
# remove weak, background peaks called as peaks by MACS below median summit count
median.remove.bg.peaks=function(file)
{
	
	summits=read.table(file,head=F,sep="\t")
	sig.summits=summits[summits$V5>median(summits$V5),]
	write.table(sig.summits,paste(strsplit(file,".bed")[[1]],"_significant.bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
	# return absolute file name
	return(paste(strsplit(file,".bed")[[1]],"_significant.bed",sep=""))
	#return(paste(unlist(strsplit(strsplit(file,".bed")[[1]],"/"))[[length(unlist(strsplit(strsplit(file,".bed")[[1]],"/")))]],"_significant.bed",sep=""))
}


# merge.uneven.cdt, merges two cdt files of different row lenghts (sometimes occurring with annotatePeaks.pl and wig files
# splits resulting merge in two
# writes merged files with ".merged.cdt" extension

merge.uneven.cdt=function(file1,file2)

{

	#read cdt files
	x=read.table(file1,head=T,sep="\t")
	y=read.table(file2,head=T,sep="\t")

	#merge cdt files, filling up missing values between each other with NA
	x.y.union=merge(x,y,by.x="Gene",by.y="Gene",all.x=T,all.y=T)

	#give merged files rownames from "Gene"
	rownames(x.y.union)=x.y.union$Gene


	#delete Gene column
	x.y.union$Gene=NULL


	#get merged table column length and split table in two based on table length/2
	cdt.range=length(colnames(x.y.union))
	x.y.union.x=x.y.union[,1:(cdt.range/2)]
	x.y.union.y=x.y.union[,(1+cdt.range/2):cdt.range]

	# replace NA with 0
	x.y.union.x[is.na(x.y.union.x)]=0
	x.y.union.y[is.na(x.y.union.y)]=0

	# write files
	write.table(x.y.union.x, paste(strsplit(file1,".cdt"),".merged.cdt",sep="") ,col.names=NA,row.names=T,sep="\t",quote=F)
	write.table(x.y.union.y, paste(strsplit(file2,".cdt"),".merged.cdt",sep="") ,col.names=NA,row.names=T,sep="\t",quote=F)

}

main=function(file1,file2,wigFile1,wigFile2,peakDist=400,profileSize=2000,countsRange=400,fcCutoff=1,peakFilteringMethod="none",cutoffMethod="foldChange",blackList="",simpleRepeats="",genome="hg38",histSize=10,outdir="")
{
	# set HOMER environment variables
	Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/package/homer/","/package/homer/bin/","/package/bedtools2/bin/",sep=":"))
	
	if(nchar(outdir)>0)
	{
		system(paste("mkdir ", outdir,sep=""))
		setwd(outdir)
	}
	
	# default no filtering
	if(peakFilteringMethod=="none")
	{
		sigSummits1=file1
		sigSummits2=file2
	}

	# if gaussian filtering set, get significant MACS summits based on Gaussian mixture model, retaining high peaks only, represented by separate population
	if(peakFilteringMethod=="gaussian")
	{
		sigSummits1=gaussian.fit.remove.bg.peaks(file1)
		sigSummits2=gaussian.fit.remove.bg.peaks(file2)
	}
	
	# if median filtering set, selects only summits for which counts > median
	if(peakFilteringMethod=="median")
	{
		sigSummits1=median.remove.bg.peaks(file1)
		sigSummits2=median.remove.bg.peaks(file2)
	}
		
	# get relative summit names
	sigSummitsName1=(paste(unlist(strsplit(strsplit(sigSummits1,".bed")[[1]],"/"))[[length(unlist(strsplit(strsplit(sigSummits1,".bed")[[1]],"/")))]],"_significant.bed",sep=""))
	sigSummitsName2=(paste(unlist(strsplit(strsplit(sigSummits2,".bed")[[1]],"/"))[[length(unlist(strsplit(strsplit(sigSummits2,".bed")[[1]],"/")))]],"_significant.bed",sep=""))

	
	# compute union of significant peaks, with 2 summits considered a same summit if within peakDist. Optionally overlaps with blacklist.
	if(blackList=="")
	{
		system(paste("cat ",sigSummits1," ",sigSummits2," | bedtools sort -i - | bedtools merge -i - -d ",peakDist," > union_",sigSummitsName1,"_",sigSummitsName2,"_union_",peakDist,".bed",sep=""))
		unionName=paste("union_",sigSummitsName1,"_",sigSummitsName2,"_dist_",peakDist,".bed",sep="")
	}
	if((blackList!="")&(simpleRepeats==""))
	{
		system(paste("cat ",sigSummits1," ",sigSummits2," | bedtools sort -i - | bedtools merge -i - -d ",peakDist,"| bedtools intersect -a - -b ",blackList," -v > union_",sigSummitsName1,"_",sigSummitsName2,"_union_",peakDist,"_noBlackList.bed",sep=""))
		unionName=paste("union_",sigSummitsName1,"_",sigSummitsName2,"_dist_",peakDist,"_noBlackList.bed",sep="")
	
	}
	if((blackList!="")&(simpleRepeats!=""))
	{
		system(paste("cat ",sigSummits1," ",sigSummits2," | bedtools sort -i - | bedtools merge -i - -d ",peakDist," | bedtools intersect -a - -b ",blackList," -v | bedtools intersect -a - -b ",simpleRepeats," -v > union_",sigSummitsName1,"_",sigSummitsName2,"_union_",peakDist,"_noBlackList_noSimpleRepeats.bed",sep=""))
		unionName=paste("union_",sigSummitsName1,"_",sigSummitsName2,"_union_",peakDist,"_noBlackList_noSimpleRepeats.bed",sep="")
	}
	
	
	
	# read resulting union file and give name based on chr:start:end, score and strand to conform to bed6 norm
	unionSummits=read.table(unionName,sep="\t")
	unionSummits$V4=paste(unionSummits$V1,":",unionSummits$V2,"-",unionSummits$V3,sep="")
	unionSummits$V5=0
	unionSummits$V6="+"
	write.table(unionSummits,unionName,row.names=F,col.names=F,quote=F,sep="\t")
	wigName1=(unlist(strsplit(wigFile1,"/"))[[length(unlist(strsplit(wigFile1,"/")))]])
	wigName2=(unlist(strsplit(wigFile2,"/"))[[length(unlist(strsplit(wigFile2,"/")))]])

	# gunzip wig files if gzipped, create coverage matrices using HOMER annotatePeaks.pl as cdt files
	if(length(grep(".gz",wigFile1))>0)
	{
		cat("Decompressing gzipped wig file 1\n")
		system(paste("gunzip -k -f ",wigFile1,sep=""))
		wigName1=(strsplit(wigName1,".gz"))[[1]]
		wigFile1=(strsplit(wigFile1,".gz"))[[1]]
	}
	if(length(grep(".gz",wigFile2))>0)
	{
		cat("Decompressing gzipped wig file 2\n")
		system(paste("gunzip -k -f ",wigFile2,sep=""))
		wigName2=(strsplit(wigName2,".gz"))[[1]]
		wigFile2=(strsplit(wigFile2,".gz"))[[1]]

	}
	
	
	system(paste("annotatePeaks.pl ",unionName," ",genome," -hist ",histSize," -ghist -wig ",wigFile1," -size ",profileSize," > ",unionName,"_",wigName1,".cdt",sep=""))
	system(paste("annotatePeaks.pl ",unionName," ",genome," -hist ",histSize," -ghist -wig ",wigFile2," -size ",profileSize," > ",unionName,"_",wigName2,".cdt",sep=""))
	
	# perform left outer join on cdts as can have different lengths
	merge.uneven.cdt(paste(unionName,"_",wigName1,".cdt",sep=""),paste(unionName,"_",wigName2,".cdt",sep=""))
	
	# read cdts, compute total coverages in countsRange around centre
	cdt1=read.table(paste(unionName,"_",wigName1,".merged.cdt",sep=""),row.names=1,head=T)
	cdt2=read.table(paste(unionName,"_",wigName2,".merged.cdt",sep=""),row.names=1,head=T)
	
	cdt1$coverages=rowSums(cdt1[(round(length(colnames(cdt1))/2,0)+1-0.5*countsRange/histSize):(round(length(colnames(cdt1))/2,0)+1+0.5*countsRange/histSize)])
	cdt2$coverages=rowSums(cdt2[(round(length(colnames(cdt2))/2,0)+1-0.5*countsRange/histSize):(round(length(colnames(cdt2))/2,0)+1+0.5*countsRange/histSize)])
	
	# replace 0s by 1s for log2
	cdt1$coverages=replace(cdt1$coverages,cdt1$coverages==0,1)
	cdt2$coverages=replace(cdt2$coverages,cdt2$coverages==0,1)

	# compute log2 cdt2/cdt1 FC and sort cdts by FC
	cdt1$FC=log(cdt2$coverages/cdt1$coverages,2)
	cdt2$FC=log(cdt2$coverages/cdt1$coverages,2)
	cdt1=cdt1[order(cdt1$FC),]
	cdt2=cdt2[order(cdt2$FC),]
	
	# write sorted bed file for further operations, remove coverages and FC then write sorted cdts
	write.table(cbind((matrix(unlist(strsplit(rownames(cdt1),"[:-]")),ncol=3,byrow=T)),rownames(cdt1),cdt1$FC,"+"),paste(unionName,"_FC_sorted.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
	cdt1$coverages=NULL
	cdt1$FC=NULL
	cdt2$coverages=NULL
	cdt2$FC=NULL
	write.table(data.frame(rownames(cdt1),cdt1), paste(unionName,"_",wigName1,".merged.cdt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
	write.table(data.frame(rownames(cdt2),cdt2), paste(unionName,"_",wigName2,".merged.cdt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
	


}












