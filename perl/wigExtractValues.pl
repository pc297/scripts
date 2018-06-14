#!/usr/bin/perl
# wigExtractValues.pl by Pierre Cauchy 2010
# takes sorted sgr coordinate file as chrX [0-9]* [+/-], retrieves values from wig file +/- window, every step bp

use strict;

if((@ARGV[0] eq "")||(@ARGV[0] eq "--help")||(@ARGV[0] eq "-h")||(@ARGV[1] eq "")||(@ARGV[2] eq "")||(@ARGV[3] eq ""))
{
	die "usage: wigExtractValues.pl wig_file coordinates_sgr_file window step";
}
open(INFILE,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
open(COORDS,"@ARGV[1]") or die "could not open @ARGV[1] for reading";

@ARGV[0]=~m/(.*)\./;
my $name1=$1;
@ARGV[1]=~m/\\([a-zA-Z0-9-_]*)\./;
@ARGV[1]=~m/\/([a-zA-Z0-9-_]*)\./;
my $name2=$1;
my $window=@ARGV[2];
@ARGV[3]=~m/([0-9]*)/;
my $step=$1;

my $filename=$name1."_".$name2."[-$window;$window]_bin$step.sgr";
my $expname=$1;
my @filenames=split(/\//,$expname);
my $size=@filenames;
$expname=@filenames[$size-1];
open(COORDS,"@ARGV[1]") or die "could not open @ARGV[0] for reading";



print "\n\n@ARGV[0]\n\n$filename\n\n";
open(OUTFILE,">$filename") or die "could not open $filename for writing";

my $oldcoordstart=0;
my $oldcoordchr="";
my @coordarray=();
my $wigstep;
my $chr;
my $start;
<INFILE>;

#recognise if file has no track header
my $header=<INFILE>;
if($header=~m/track/)
{
	<INFILE>=~m/chrom=(chr[0-9_a-zA-Z]*) start=([0-9]*) step=([0-9]*)/;
	$chr=$1;
	$start=$2;
	$wigstep=$3;
}
else
{
	$header=~m/chrom=(chr[0-9_a-zA-Z]*) start=([0-9]*) step=([0-9]*)/;
	$chr=$1;
	$start=$2;
	$wigstep=$3;
}

print "$chr\t$start\t$step\n";

#go through coordinate file
while(my $coord=<COORDS>)
{
	chomp($coord);
	$coord=~m/(chr[0-9a-zA-Z_-]*)\t([0-9]*)\t([+-])/;
	my $coordchr=$1;
	my $coordstart=$2;
	my $coordstrand=$3;
	if(($coordchr eq $oldcoordchr)&&($coordstart-$window<=$oldcoordstart+$window))
	{
	print "OVERFLOW";
	}
	else
	{
		for(my $i=$coordstart-$window;$i<$coordstart+$window+$step;$i+=$wigstep)
		{
		
			my $distance;
			if($coordstrand eq "+")
			{
				$distance=$i-$coordstart;
			}
			if($coordstrand eq "-")
			{
				$distance=$coordstart-$i;
			}
			
			my @coordline=($coordchr,$i,$coordstrand,$distance,$coordstart);
			push(@coordarray, [@coordline]);
		}
	}
	$oldcoordchr=$coordchr;
	$oldcoordstart=$coordstart;
	
}	
my $max=@coordarray;

my $score=0;
my $i=0;
print "max: $max\n";

print OUTFILE "COORD";
for(my $j=0;$j<=2*$window;$j+=$step)
{
	my $coord=$j-$window;
	print OUTFILE "\t$coord";
}

my $oldcoord="";
my $oldstrand="";
my $oldscore=0;
my @region=();

while(my $line=<INFILE>)
{
	
	if($line=~m/chrom/)
	{
		$line=~m/chrom=(chr[0-9_a-zA-Z]*) start=([0-9]*) step=([0-9]*)/;
		$chr=$1;
		$start=$2;
		$wigstep=$3;
		print "$chr\t$start\n";
	}
	
	if(!($line=~m/track/)&&!($line=~m/chrom/))
	{
		
		
		
		if(($chr eq $coordarray[$i][0])&&($start==$coordarray[$i][1]))
		{
			
			if(!("$coordarray[$i][0]:$coordarray[$i][4]" eq $oldcoord))
			{
				
				#if input step lower than wig step perform interpolation
				if($step<$wigstep)
				{
					
					
					if($oldstrand eq "+")
					{
						
						for(my $j=0;$j<2*$window/$wigstep;$j++)
						{
							my $slope=($region[$j+1]-$region[$j])/$wigstep;
							for(my $k=0;$k<$wigstep;$k+=$step)
							{
								
								my $interpolation=$region[$j] + $k*$slope;
								print OUTFILE "\t$interpolation";
				
							}
							
							
							
							
						}
						print OUTFILE "\t$region[2*$window/$wigstep]";
						
						
					}
					if($oldstrand eq "-")
					{
						
						
						for(my $j=2*$window/$wigstep;$j>0;$j--)
						{
							my $slope=($region[$j-1]-$region[$j])/$wigstep;
							for(my $k=0;$k<$wigstep;$k+=$step)
							{
								my $interpolation=$region[$j] + $k*$slope;
								print OUTFILE "\t$interpolation";
				
							}
						}
						print OUTFILE "\t$region[0]";
						
					}
					
					
					
					
				
				}
				
				
				
				
				else
				{
					
					
					
					if($oldstrand eq "+")
					{
						for(my $j=0;$j<=2*$window/$wigstep;$j++)
						{
							print OUTFILE "\t$region[$j]";
						}
					}
					
					
					if($oldstrand eq "-")
					{
						for(my $j=2*$window/$wigstep;$j>=0;$j--)
						{
							print OUTFILE "\t$region[$j]";
						}
					}
				
					
				
				}
				print OUTFILE "\n$coordarray[$i][0]:$coordarray[$i][4]";
				$oldstrand=$coordarray[$i][2];
				$oldcoord="$coordarray[$i][0]:$coordarray[$i][4]";
				$oldscore=0;
				@region=();
			}
			chomp($line);
			$line=~s/ //;
			$line=~s/\r//;
			push(@region, $line);
			
			print "$chr\t$start\t$line\t$coordarray[$i][2]\t$coordarray[$i][3]\n";
			$i++;
		}
		
		$start+=$wigstep;
		$oldscore=$line;
	}
	
}
