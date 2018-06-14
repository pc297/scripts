# bootstrapMotifs by Pierre Cauchy (c) 2014
# computes Z-score of motif co-occurrence enrichments in sample bed file versus randomised control file
# requires: R packages abind, gplots; external bedtools, Homer and pybedtools

# main function:
# sampleBedFileName, backgroundBedFileName, genomeChromSizesFileName as filenames with absolute or relative path
# workingDir, bedToolsPath, pyBedToolsPath as directory names with absolute or relative path
# motif list as list of motifs using absolute or relative filenames, e.g. c("/motifs/ap1.motif, will be copied as ./boostrapping/motifs/motif1.motif etc
# genome as string for Homer, e.g. "mm10", can be path to homer genome e.g. "/genomes/mm10"
# peakSize as size of peaks to map motifs: "given" for footprints, or 200 (default for Homer) for ChIP/DNaseI/ATAC-Seq peaks
# extSize as +/- extension size used for intersection of motifs: 25 for footprints defaults (then looks within +/- 25bp), recommend 100 for ChIP/DNaseI/ATAC-Seq peaks
# repetitions: 1000 iterations of random samplings of size recommended
# outputColorRampPalette as colorRampPalette(c("grey","white","brown"))(300)
bootstrapMotifs=function(sampleBedFileName, backgroundBedFileName, motifList="", genome, genomeChromSizesFileName, extSize=25, peakSize="given", repetitions=1000, workingDir=getwd(), homerPath=Sys.getenv("PATH"), bedToolsPath=Sys.getenv("PATH"), pyBedToolsPath=Sys.getenv("PATH"), outputColorRampPalette=colorRampPalette(c("grey","white","brown"))(300))
{
	library("abind")
	library("gplots")

	#system environment vars, files and directory preparation
	Sys.setenv(PATH=paste(Sys.getenv("PATH"),homerPath, paste(homerPath, "/bin",sep=""), bedToolsPath, pyBedToolsPath, sep=":"))
	setwd(workingDir)
	system("mkdir bootstrapping")
	system("mkdir bootstrapping/motifs")
	system("rm bootstrapping/motif*.bed",intern=T)
	motifNum=1
	for(motif in motifList)
	{
		system(paste("cp", motif, paste(workingDir, "/bootstrapping/motifs/motif",motifNum,".motif", sep=""),sep=" "))
		motifNum=motifNum+1
	}

	#read bed files for background and sample, sets sampleSize from sample size 
	backgroundBed=read.table(backgroundBedFileName,header=F,sep="\t")
	sampleSize=nrow(read.table(sampleBedFileName,header=F,sep="\t"))

	#bootstrapping doing repetitions iterations of random samplings of sampleSize size in background, motif mapping, create intersection matrix everytime
	for(i in 1:repetitions)
	{
		#create random intervals from control file
		rand=sample(nrow(backgroundBed))
		randBackgroundBed=backgroundBed[rand,]
		sampleRandBackgroundBed=head(randBackgroundBed,sampleSize)
		write.table(sampleRandBackgroundBed,"sample_rand.bed",col.names=F,row.names=F,sep="\t",quote=F)
		motifs=system("ls bootstrapping/motifs/",intern=T)
	
		#use homer annotatePeaks to map motifs in random regions; use length(motifs)-1) here as system commands ended with & thus running in parallel and not waiting completion
		for(j in 1:(length(motifs)-1))
		{
			motif_split=strsplit(motifs[j],"\\.")
			motif_name=motif_split[[1]][[1]]		
			system(paste("annotatePeaks.pl sample_rand.bed ",genome," -size ", peakSize, " -m bootstrapping/motifs/",motifs[j]," -mbed ", workingDir,"/bootstrapping/",motif_name,".bed -noann > bootstrapping/output_tmp.txt &",sep=""))
		}
		
		j=length(motifs)
		motif_split=strsplit(motifs[j],"\\.")
		motif_name=motif_split[[1]][[1]]	
			
		#use homer annotatePeaks to map motifs in random regions; last iteration here as previous system commands ended with & thus running in parallel and not waiting completion; this command has no & so waiting for completion
		system(paste("annotatePeaks.pl sample_rand.bed ",genome," -size ", peakSize, " -m bootstrapping/motifs/",motifs[j]," -mbed ", workingDir,"/bootstrapping/",motif_name,".bed -noann > bootstrapping/output_tmp.txt",sep=""))
	
		#extend motif beds to motif length -1
		for(j in 1:(length(motifs)-1))
		{
			motif_split=strsplit(motifs[j],"\\.")
			motif_name=motif_split[[1]][[1]]
			system(paste("bedtools slop -i ", workingDir, "/bootstrapping/",motif_name,".bed -g ", genomeChromSizesFileName, " -b ",extSize," > ", workingDir,"/bootstrapping/",motif_name,"_",extSize,"_extended.bed",sep=""))
		}

		#extend last motif to make sure system command of intersection is carried out after completion
		j=length(motifs)
		motif_split=strsplit(motifs[j],"\\.")
		motif_name=motif_split[[1]][[1]]
		system(paste("bedtools slop -i ", workingDir, "/bootstrapping/",motif_name,".bed -g ", genomeChromSizesFileName, " -b ",extSize," > ", workingDir,"/bootstrapping/",motif_name,"_",extSize,"_extended.bed",sep=""))

		#perform intersection
		system(paste("pybedtools intersection_matrix ", workingDir, "/bootstrapping/*_",extSize,"_extended.bed > ",workingDir, "/bootstrapping/intersection_round",i,".txt",sep=""))

		#remove motif bed files
		system("rm bootstrapping/motif*.bed",intern=T)
	}

	#perform actual intersection from sample file
	#use homer annotatePeaks to map motifs in random regions; use length(motifs)-1) here as system commands ended with & thus running in parallel and not waiting completion
	for(j in 1:(length(motifs)-1))
	{
		motif_split=strsplit(motifs[j],"\\.")
		motif_name=motif_split[[1]][[1]]
		system(paste("annotatePeaks.pl ", sampleBedFileName, " ", genome," -size ", peakSize, " -m bootstrapping/motifs/",motifs[j]," -mbed ", workingDir, "/", motif_name,".bed -noann > bootstrapping/output_tmp.txt &",sep=""))
	}
	
	#use homer annotatePeaks to map motifs in random regions; last iteration here as previous system commands ended with & thus running in parallel and not waiting completion; this command has no & so waiting for completion
	j=length(motifs)
	motif_split=strsplit(motifs[j],"\\.")
	motif_name=motif_split[[1]][[1]]	
	system(paste("annotatePeaks.pl ", sampleBedFileName, " ", genome," -size ", peakSize, " -m bootstrapping/motifs/",motifs[j]," -mbed ", workingDir, "/", motif_name,".bed -noann > bootstrapping/output_tmp.txt",sep=""))

	#extend motif beds to motif length -1
	for(j in 1:(length(motifs)-1))
	{
		motif_split=strsplit(motifs[j],"\\.")
		motif_name=motif_split[[1]][[1]]
		system(paste("bedtools slop -i ", workingDir, "/bootstrapping/",motif_name,".bed -g ", genomeChromSizesFileName, " -b ",extSize," > ", workingDir,"/bootstrapping/",motif_name,"_",extSize,"_extended.bed",sep=""))
	}

	#extend last motif to make sure system command of intersection is carried out after completion
	j=length(motifs)
	motif_split=strsplit(motifs[j],"\\.")
	motif_name=motif_split[[1]][[1]]
	system(paste("bedtools slop -i ", workingDir, "/bootstrapping/",motif_name,".bed -g ", genomeChromSizesFileName, " -b ",extSize," > ", workingDir,"/bootstrapping/",motif_name,"_",extSize,"_extended.bed",sep=""))
		
	#perform intersection
	system(paste("pybedtools intersection_matrix ", workingDir, "/bootstrapping/*_",extSize,"_extended.bed > ", workingDir, "/bootstrapping/intersection_actual.txt",sep=""))

	#remove motif bed files
	system("rm bootstrapping/motif*.bed",intern=T)

	#read intersections and bind in 3D array
	#read first intersection text file to initialise int_array
	int_array=read.table(paste(workingDir,"/bootstrapping/intersection_round1.txt",sep=""))
	cat("Read round 1 intersection\n")

	#recpirocate matrix as sometimes intersections off by small numbers due to multiple peak overlaps
	for(k in 1:length(motifs))
	{
		for(l in k:length(motifs))
		{
			int_array[k,l]=int_array[l,k]
		}
	}

	#read remaining intersection text files to fill int_array
	for(i in 2:repetitions)
	{
		tmp=read.table(paste(workingDir,"/bootstrapping/intersection_round",i,".txt",sep=""))
	
		#recpirocate matrix as sometimes intersections off by small numbers due to multiple peak overlaps
		for(k in 1:length(motifs))
		{
			for(l in k:length(motifs))
			{
				tmp[k,l]=tmp[l,k]
			}
		}	
		int_array=abind(int_array,tmp,along=3)
		cat(paste("Read round",i,"intersection\n",sep=" ")) 
	}

	#create arrays of mean and sd of random intersection
	mean_array=array(0,c(length(motifs),length(motifs)))
	sd_array=array(0,c(length(motifs),length(motifs)))
	for(i in 1:length(motifs))
	{
		for(j in 1:length(motifs))
		{
			mean_array[i,j]=mean(int_array[i,j,])
			sd_array[i,j]=sd(int_array[i,j,])
		}
	}
	motif_header=NULL
	intersection_actual=read.table(paste(workingDir,"/bootstrapping/intersection_actual.txt",sep=""),header=T,row.names=1,sep="\t")
	motif_names=colnames(intersection_actual)
	motif_names_short=NULL
	for(w in 1:length(motif_names))
	{
		motif_name_split=strsplit(motif_names[w],"_")
		motif_names_short=rbind(motif_names_short,motif_name_split[[1]][[1]])
	}
	for(m in 1:length(motif_names_short))
	{
		TF_name=system(paste("head -n 1 bootstrapping/motifs/",motif_names_short[m],".motif",sep=""),intern=T)
		TF_name_split=strsplit(TF_name,"\t")		
		motif_header=cbind(motif_header,TF_name_split[[1]][[2]])
	}	

	#calculate z scores
	intersection_actual=read.table(paste(workingDir,"/bootstrapping/intersection_actual.txt"),header=T,row.names=1,sep="\t")

	#recpirocate matrix as sometimes intersections off by small numbers due to multiple peak overlaps
	for(k in 1:length(motifs))
	{
		for(l in k:length(motifs))
		{
			intersection_actual[k,l]=intersection_actual[l,k]
		}
	}

	# z scores computed as z=((x+1)-(mean+1))/(sd+1) to avoid dividing by 0
	z_scores=((intersection_actual+1)-(mean_array+1))/(sd_array+1)
	colnames(z_scores)=motif_header
	rownames(z_scores)=motif_header
	write.table(z_scores,"bootstrapping/zscores.txt",quote=F,row.names=T,sep="\t",col.names=NA)
	pdf("bootstrapping/zscores.pdf")
	heatmap.2(z_scores,trace="none",density.info="none", col=outputColorRampPalette)
	dev.off()
}

