#!/usr/bin/perl
# vstep2fixed_hg19.pl
# converts hg19 continuous variable step wig to fixed step wig
# fills up gaps in vstep wig with zeros: looks for coord of first line, fills up
# start to coord of first line with zeros, and fills end of annotated region 
# to end of chrom with zeros. Does not fill gaps between annotated regions
# with zeros as there should be none, if there are use wig_vstep2fixed.pl
# useful for UCSC downloaded tracks, usually in vstep
# Usage: vstep2fixed_hg19.pl <variable_step.wig>
use strict;
open(INFILE,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
@ARGV[0]=~m/(.*)\./;

if((@ARGV[0] eq "")||(@ARGV[0] eq "-h")||(@ARGV[0] eq "--help"))
{
	die "Usage: vstep2fixed_hg19.pl <variable_step.wig>\n";
}

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
	
if($line=~m/chrom/)
	{
		$line=~m/chrom=(chr[0-9_a-zA-Z]*) span=([0-9]*)/;
		$chr=$1;
		$wigstart=1;
		$step=$2;
		
			print OUTFILE "track type=wiggle_0\nfixedStep chrom=$chr start=1 step=$2\n";
			
			if($chr eq "chr1")
			{
				$end=249250621;
			}
			if($chr eq "chr2")
			{
				$end=243199373;
			}
			if($chr eq "chr3")
			{
				$end=1998022430 ;
			}
			if($chr eq "chr4")
			{
				$end=191154276 ;
			}
			if($chr eq "chr5")
			{
				$end=180915260 ;
			}
			if($chr eq "chr6")
			{
				$end=171115067;
			}
			if($chr eq "chr7")
			{
				$end=159138663;
			}
			if($chr eq "chr8")
			{
				$end=146364022;
			}
			if($chr eq "chr9")
			{
				$end=141213431;
			}
			if($chr eq "chr10")
			{
				$end=135534747;
			}
			if($chr eq "chr11")
			{
				$end=135006516;
			}
			if($chr eq "chr12")
			{
				$end=133851895;
			}
			if($chr eq "chr13")
			{
				$end=115169878;
			}
			if($chr eq "chr14")
			{
				$end=107349540;
			}
			if($chr eq "chr15")
			{
				$end=102531392;
			}
			if($chr eq "chr16")
			{
				$end=90354753;
			}
			if($chr eq "chr17")
			{
				$end=81195210;
			}
			if($chr eq "chr18")
			{
				$end=78077248;
			}
			if($chr eq "chr19")
			{
				$end=59128983;
			}
			if($chr eq "chr20")
			{
				$end=63025520;
			}
			if($chr eq "chr21")
			{
				$end=48129895;
			}
			if($chr eq "chr22")
			{
				$end=51304566;
			}
			if($chr eq "chrX")
			{
				$end=155270560;
			}
			if($chr eq "chrY")
			{
				$end=59373566;
			}
			if($chr eq "chrM")
			{
				$end=16571;
			}
			print "$chr\tend: $end\n";
			
			
		
		if($coordstart>1)
		{
			for(my $i=$coordstart;$i<=$end;$i+=$step)
			{
				print OUTFILE "0.000\n";
			}
		}
		$coordstart=1;
	
		
						
	}
	
	
	
	
	else
	{
		
		if((!($line=~m/chrom/)&&(!($line=~m/track/))))
		{
			$line=~m/([0-9]*)\t([0-9]*)/;
			while(($coordstart<$1)&&($coordstart<$end))
			{
				print OUTFILE "0.000\n";
				$coordstart+=$step;
			}
			print OUTFILE "$2.000\n";
			$coordstart+=$step;
		}
			
	}
}
print "$coordstart\n";
if($coordstart<$end)
{
	for(my $i=$coordstart;$i<=$end;$i+=$step)
	{
		print OUTFILE "0.000\n"
	}
}



exit (0);
