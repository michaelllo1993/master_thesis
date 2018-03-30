#!/usr/bin/perl
use warnings;
use strict;
my $usage = <<EOF;

USAGE: perl seq2csv.pl <filename>

EOF

my $inFile = $ARGV[0];
my $org = substr($inFile,0,-4);
my $out_name="$org.csv";

open(IN, '<', $inFile) or die "Could not open '$inFile' $! \n";
open (OUT, ">$out_name") or die "Couldn't open: $!";
print OUT "ensembl_transcript_id".","."ensembl_peptide_id".","."START".","."STOP".","."cdna";
while (my $line = <IN>){
	chomp $line;
	if (($line =~ /^>(ENS\w+)\|(ENS\w+)\|(\S+)\|(\S+)/)) {
	  print OUT "\n".$1.",".$2.",".$3.",".$4.",";
	}
	elsif (($line =~ /^>(|)/) || ($line =~ /^>([A-Z]+\S*)\|/)) {
	  print OUT "\n".$1.","."".","."".","."".",";
	}
	elsif (($line =~ /^>(ENS\w+)/) || ($line =~ /^>(|)/) || ($line =~ /^>([A-Z]+\S*)\|/)) {
	  print OUT "\n".$1.","."".","."".","."".",";
	}
	elsif (($line !~ /^>/) && ($line !~ /^Sequence/)) {
	  print OUT $line
	}
}
close(IN);
close(OUT);
