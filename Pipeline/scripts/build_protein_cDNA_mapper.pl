#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: perl build_protein_cDNA_mapper.pl <filename>

EOF

my $inFile = $ARGV[0];
my $org = substr($inFile,0,-4);
my $outFile = "$org\_mapper.csv";

open(IN, '<', $inFile) or die "Could not open '$inFile' $! \n";
open (OUT, ">",$outFile);
printf OUT "Protein_ID".","."cDNA_ID".","."\n";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>(ENS\w+\d+)\|(ENS\w+\d+)/){
	  printf OUT $1.",".$2.","."\n";	
	}
}
close(IN);
close (OUT);
system("mv $outFile data/mappers/");


