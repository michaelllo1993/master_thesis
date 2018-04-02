#!/usr/bin/perl
use warnings;
use strict;
my $usage = <<EOF;

USAGE: perl seq2csv.pl organism

EOF

my $inFile = $ARGV[0];
system("sed -i '/Sequence\ unavailable/d' $inFile");
my $org = substr($inFile,0,-4);
my $out_name="$org.csv";

open(IN, '<', $inFile) or die "Could not open '$inFile' $! \n";
open (OUT, ">$out_name");
print OUT "ensembl_peptide_id".","."ensembl_transcript_id".","."peptide";
while (my $line = <IN>){
	chomp $line;
	if (($line =~ /^>(ENS\w+)\|(ENS\w+)/)) {
	  print OUT "\n".$1.",".$2.",";
	}
	elsif ($line =~ /^>([A-Z]+\S*)\|/){
	  print OUT "\n".$1.",";
	}
	elsif (($line !~ /^>/) && ($line !~ /^Sequence/)) {
	  print OUT $line
	}
}
