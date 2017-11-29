#!/usr/bin/perl
use strict;
open(INFILE,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
@ARGV[0]=~m/(.*)\./;
my $filename=$1."_fixed.wig";
my $expname=$1;
my @filenames=split(/\//,$expname);
my $size=@filenames;
$expname=@filenames[$size-1];
print "\n\n@ARGV[0]\n\n$filename\n\n";
open(OUTFILE,">$filename") or die "could not open $filename for writing";
my $chr;
my $coordstart=1;
my $wigstart;
my $step;
my $end=0;
while(my $line = <INFILE>)
{
	if(($line=~m/chrom/)&&(!($line=~m/#/)))
	{
		$line=~m/chrom=(chr[0-9_a-zA-Z]*) span=([0-9]*)/;
		$chr=$1;
		$wigstart=1;
		$step=$2;
		if($coordstart>1)
		{
			for(my $i=$coordstart;$i<$end;$i+=$step)
			{
				print OUTFILE "0\n";
			}
		}
		$coordstart=1;		


			print OUTFILE "track type=wiggle_0\nfixedStep chrom=$chr start=1 step=$2\n";
			
			if($chr eq "chr1")
			{
				$end=197195432;
			}
			if($chr eq "chr2")
			{
				$end=181748087;
			}
			if($chr eq "chr3")
			{
				$end=159599783 ;
			}
			if($chr eq "chr4")
			{
				$end=155630120 ;
			}
			if($chr eq "chr5")
			{
				$end=152537259 ;
			}
			if($chr eq "chr6")
			{
				$end=149517037;
			}
			if($chr eq "chr7")
			{
				$end=152524553;
			}
			if($chr eq "chr8")
			{
				$end=131738871;
			}
			if($chr eq "chr9")
			{
				$end=124076172;
			}
			if($chr eq "chr10")
			{
				$end=129993255;
			}
			if($chr eq "chr11")
			{
				$end=121843856;
			}
			if($chr eq "chr12")
			{
				$end=121257530;
			}
			if($chr eq "chr13")
			{
				$end=120284312;
			}
			if($chr eq "chr14")
			{
				$end=125194864;
			}
			if($chr eq "chr15")
			{
				$end=103494974;
			}
			if($chr eq "chr16")
			{
				$end=98319150;
			}
			if($chr eq "chr17")
			{
				$end=95272651;
			}
			if($chr eq "chr18")
			{
				$end=90772031;
			}
			if($chr eq "chr19")
			{
				$end=61342430;
			}
			if($chr eq "chrX")
			{
				$end=166650296;
			}
			if($chr eq "chrY")
			{
				$end=15902555;
			}
			if($chr eq "chrM")
			{
				$end=16299;
			}
			print "$chr\tend: $end\n";
			
			
		
		
	
		
						
	}
	
	
	
	
	else
	{
		
		if((!($line=~m/chrom/))&&(!($line=~m/track/))&&(!($line=~m/#/)))
		{
			$line=~m/([0-9]*)\t([0-9\.]*)/;
			# modified			
			#while(($coordstart<$1+1)&&($coordstart<$end))
			while(($coordstart<$1)&&($coordstart<$end))
			{
				print OUTFILE "0\n";
				$coordstart+=$step;
			}
			
				print OUTFILE "$2\n";
				$coordstart+=$step;
			
		}
			
	}
}
print "$coordstart\n";
if($coordstart<$end)
{
	for(my $i=$coordstart;$i<=$end;$i+=$step)
	{
		print OUTFILE "0\n"
	}
}



exit (0);
