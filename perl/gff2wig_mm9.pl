#!/usr/bin/perl
# gff2wig_mm9.pl by Pierre Cauchy 2009
# Produces wig file from sorted sgr file having extended read start/end information as 3d column (1=start, 0=end)
# mm9 version
use strict;

@ARGV[0]=~m/(.*)\./;
open(GFF,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
my $step=@ARGV[1];
my $filename=$1."_bin$step.wig";


my $expname=$1;
my @filenames=split(/\\/,$expname);
my $size=@filenames;
$expname=@filenames[$size-1];
open(OUTFILE,">$filename") or die "could not open @ARGV[0] for writing";
print "\n\n@ARGV[0]\n\n$filename\n\n";






my $start=1;
my $chr="";








my $k=0;


my @chromarray=(["chr1",197195432],["chr10",129993255],["chr11",121843856],["chr12",121257530],["chr13",120284312],["chr13_random",400311],["chr14",125194864],["chr15",103494974],["chr16",98319150],["chr16_random",3994],["chr17",95272651],["chr17_random",628739],["chr18",90772031],["chr19",61342430],["chr1_random", 1231697],["chr2",181748087],["chr3",159599783],["chr3_random",41899],["chr4",155630120],["chr4_random",160594],["chr5",152537259],["chr5_random",357350],["chr6",149517037],["chr7",152524553],["chr7_random",362490],["chr8",131738871],["chr8_random",849593],["chr9",124076172],["chr9_random",449403],["chrM",16299],["chrUn_random",5900358],["chrX",166650296],["chrX_random",1785075],["chrY",15902555],["chrY_random",58682461]);



my $max1=@chromarray;

my $line=<GFF>;
$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
my $gffchr=$1;
my $gffstart=$2;
my $gffflag=$3;

print "start=$gffchr:$gffstart\n";
for(my $i=0;$i<$max1;$i++)
{
	print "$chromarray[$i][0]:1\nCurrent position: $gffchr:$gffstart\n";
	#if neg coords generated during extension
	if($gffstart<=0)
	{
		$gffstart=1;
	}
	
	
	
	print OUTFILE "track type=wiggle_0\nfixedStep chrom=$chromarray[$i][0] start=1 step=$step\n";
	
	
	
	my $score=0;
	for(my $j=0;$j+$step-1<$chromarray[$i][1];$j+=$step-1)
	{
		
		$j++;
		if(($chromarray[$i][0] eq "chr13")&&($gffchr eq "chr1"))
		{
			print "ERROR: $chromarray[$i][0]:$j-".($j+$step-1)."\nCurrent position: $gffchr:$gffstart\n";
			exit(0);
			
		}
			
				
				
				
				if(($chromarray[$i][0] eq $gffchr)&&($gffstart>=$j))
				{
				
					
					while(($gffstart<$j+$step)&&($chromarray[$i][0] eq $gffchr))
					{
						if($gffflag==1)
						{
					
							$score+=1;
						}
						if(($gffflag==0)&&($score>0))
						{
					
							$score-=1;
						}
						
						
												
						$line=<GFF>;
						$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
						$gffchr=$1;
						
						$gffstart=$2;
						$gffflag=$3;
						
						while($gffchr eq "chrL")
						{
							$line=<GFF>;
							$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
							$gffchr=$1;
						
							$gffstart=$2;
							$gffflag=$3;
						}
				
						if($gffstart<=0)
						{
							while($gffstart<=0)
							{
								print "Overflow: $gffchr:$gffstart is negative\n";
								print "Current interval: $chromarray[$i][0]:$j-".($j+$step-1)."\n";
								if($line=<GFF>)
								{	
									$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
									$gffchr=$1;
									$gffstart=$2;
									$gffflag=$3;
									while($gffchr eq "chrL")
									{
										$line=<GFF>;
										$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
										$gffchr=$1;
						
										$gffstart=$2;
										$gffflag=$3;
									}
									print "Next position read: $gffchr:$gffstart\n";
								}
							}
						}
						#check if some elongated tags are outside chromosomal coords
						if($gffstart>=1+$step*(int($chromarray[$i][1]/$step)))
						{
							while(($chromarray[$i][0] eq $gffchr))
							{
								print "Overflow: $gffchr:$gffstart is greater than chromosome maximum $chromarray[$i][0]:$chromarray[$i][1]\n";
								print "Current interval: $chromarray[$i][0]:$j-".($j+$step-1)."\n";
								if($line=<GFF>)
								{	
									$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
									$gffchr=$1;
									$gffstart=$2;
									$gffflag=$3;
									while($gffchr eq "chrL")
									{
										$line=<GFF>;
										$line=~m/(chr[0-9_a-zA-Z]*)\t[ ]*([-0-9]*)\t([0-9]*)/;
										$gffchr=$1;
						
										$gffstart=$2;
										$gffflag=$3;
									}


									print "Next position read: $gffchr:$gffstart\n";
								}
								else
								{
									exit(0);
								}
							}
							last;
						}
						
						
						
						
						
						
				
					}
				
				}
				
			print OUTFILE "$score.0\n";	
			}
			
			
			
		
			
			
			
		
		
		
		
	}

	
	
	
	
	
	
	
	
	
	































		

