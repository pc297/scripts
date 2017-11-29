#!/usr/bin/perl
# matrixScanDistanceToCDT.pl
# Uses processed RSATools matrix scan input (from matrix_scan_dist-R) to create TreeView compatible matrix
# Usage: matrixScanDistanceCreateTreeViewMatrices.pl <distancefile> range_around_seq_midpoint binsize


# read args and create output file name
use strict;
@ARGV[0]=~m/(.*)\./;


my $filename=$1."_MATRIX.cdt";
my $range=@ARGV[1];
my $step=@ARGV[2];

if((@ARGV[0] eq "")||(@ARGV[1] eq "")||(@ARGV[2] eq "")||(@ARGV[0] eq "--help")||(@ARGV[1] eq "--help")||(@ARGV[0] eq "--help")||(@ARGV[0] eq "-h")||(@ARGV[1] eq "-h")||(@ARGV[0] eq "-h"))
{
	die "\nUsage: matrixScanDistanceCreateTreeViewMatrices.pl <distancefile> range_around_seq_midpoint binsize\n\n"
}

# open input distance file and output file

open(DISTANCES,"@ARGV[0]") or die "could not open @ARGV[0] for reading";
open(MATRIX, ">$filename") or die "could not open $filename for writing";

# create matrix of size 2*range_around_seq_midpoint, which equals sequence size, read distances file, and increments each position for each sequence.
# multiple sequence names 100% ok and expected, will read sequence primary key and distance (adjusted + half sequence length for matrix colums
# , e.g. if distance is -195bp, in sequences 400bp large, then index is -195+400/2=5

my @matrix;

# create array for unique sequence name
my @header;

# read distance file line by line
while(my $line=<DISTANCES>)
{
	# map sequence key, distance and sequence ID	
	$line=~m/([0-9]*)\t([0-9-]*)\t(chr[a-zA-Z0-9_:-]*)/;
	my $order=$1;
	my $distance=$2;
	my $linename=$3;

	
	# treat sequence identifiers differently, they are however necessary to get full heatmap as all sequences may not have motif
	if($distance<999999)

	
	{
		
		$matrix[$order][int(($distance+$range)/$step)]++;
		
	}
	
	# increment matrix at index indicated by sequence num and position
	else
	{	
		$matrix[$order][int(2*+$range+1)]++;
	}
	@header[$order]=$linename;
	
}


# output matrix size

my $matrix=@matrix;
print "\nmatrix size: $matrix\n";

# print first header lines of TreeView cdt format, each bin corresponding to column

print MATRIX "GID\torder\tNAME\tGWEIGHT\t";
for(my $j=0;$j<=2*$range;$j+=$step)
{
	my $coord=$j-$range;
	print MATRIX "$coord\t";
}
print MATRIX "\nAID\t\t\t\t";

for(my $j=0;$j<=2*$range;$j+=$step)
{
	my $coord=$j+1;
	print MATRIX "ARRAY$j"."X\t";
}

print MATRIX "\nEWEIGHT\t\t\t\t";
for(my $j=0;$j<=2*$range;$j+=$step)
{
	
	print MATRIX "1.000000\t";
}

# go through motif matrix and write

print MATRIX "\n";
for(my $i=1;$i<$matrix;$i++)
{
	print MATRIX "GENE"."$i"."X\t@header[$i]\t$i\t1.000000\t";
	for(my $j=0;$j<=(2*$range)/$step;$j++)
	{
		
		if($matrix[$i][$j] eq "")
		{
			print MATRIX "0\t";
		}
		else
		{
			print MATRIX "$matrix[$i][$j]\t";
		}
	}
	print MATRIX "\n";
}

close(MATRIX);

die "\n\nSuccessfully wrote TreeView CDT file $filename\n\n"

