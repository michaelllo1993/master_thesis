#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: ensembl_peptide_seq_cut.pl organism

EOF

my $file = "$ARGV[0].txt" || die $usage;
print "\nFile $file loaded successfully\n";
my $seq = "PEPTIDE SEQUENCES";
my $ids = "";
my $out_name="shortened_seq_$ARGV[0]";
my $counter = 0;
my $num = 0; 
open(my $raw_data, '<', $file) or die "Could not open '$file' $! \n";
open (OUT, ">$num$out_name");	
while (my $line = <$raw_data>){
	chomp $line;
	if ($line =~ /^>ENS(\w*)P/){
		$seq = substr $seq, 0, 69;
		print OUT "$seq\n";
		$ids = $line;
		print OUT "$ids\n";
		$seq = "";
		$counter++;
	}
	elsif (($line =~ /^[A-Z]+/) && ($line !~ /^Sequence/)) {
		$seq = $seq . $line;
	}
	if ($counter > 1000){
		$num++;
		$counter = 0;
		close (OUT);
		open (OUT, ">$num$out_name");
	}
}
$seq = substr $seq, 0, 69;
print OUT "$seq\n";
close (OUT);

print "\nNumber of files: $num \n\n";