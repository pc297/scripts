#!/usr/bin/perl
# bedGraphToWig_mm9.pl by Pierre Cauchy
# Converts bedGraph to fixed step wig file, for mm9
use strict;
@ARGV[0]=~m/(.*)\./;
open(BDG,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
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
my @chromarray=(["chr1",197195432],["chr2",181748087],["chr3",159599783],["chr4",155630120],["chr5",152537259],["chr6",149517037],["chr7",152524553],["chr8",131738871],["chr9",124076172],["chr10",129993255],["chr11",121843856],["chr12",121257530],["chr13",120284312],["chr14",125194864],["chr15",103494974],["chr16",98319150],["chr17",95272651],["chr18",90772031],["chr19",61342430],["chrX",166650296],["chrY",15902555]);
my $max1=@chromarray;
my $line=<BDG>;
$line=<BDG>;
$line=~m/(chr[0-9_a-zA-Z]*)\t([0-9_a-zA-Z -]*)\t([0-9_a-zA-Z -]*)\t([-0-9\.]*)/;
my $bdgchr=$1;
my $bdgstart=$2;
my $bdgend=$3;
my $bdgscore=$4;
my $score=$bdgscore;
for(my $i=0;$i<$max1;$i++)
{
	print "$chromarray[$i][0]\n";
	print OUTFILE "track type=wiggle_0\nfixedStep chrom=$chromarray[$i][0] start=1 step=$step\n";
	for(my $j=0;$j<$chromarray[$i][1];$j+=$step-1)
	{
		$j++;
		
		if(!((($chromarray[$i][0] eq $bdgchr)&&($j+$step-1>=$bdgstart))))
		{
			print OUTFILE "0\n";
			
		}
		else
		{
			if(($chromarray[$i][0] eq $bdgchr)&&($j+$step-1>=$bdgstart))
			
			{	
				$score+=$bdgscore;
				if($bdgend>$j-1)
				{	
					while($bdgend>$j-1)
					{
						print OUTFILE "$score\n";
						$j+=$step;
					}
					$j-=$step;
					$score+=$bdgscore;
				}
				else
				{
					print OUTFILE "$score\n";
					$score=0;
				}
				
			}
			if($line=<BDG>)
			{
				$line=~m/(chr[0-9_a-zA-Z]*)\t([0-9_a-zA-Z -]*)\t([0-9_a-zA-Z -]*)\t([-0-9\.]*)/;
				$bdgchr=$1;
				$bdgstart=$2;
				$bdgend=$3;
				$bdgscore=$4;
				$score=$bdgscore;
						
			}
		}
	}
}
exit(0);

























		
