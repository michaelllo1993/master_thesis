#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: 01_prep_and_run_signalp.pl filename

EOF

#This script prepares the protein data (input file) for the SignalP4.1 analysis and performs one

my $file_name = "$ARGV[0]" || die $usage;
my @splitted = split '\.', $file_name;
my @organism_name = split '\/', $splitted[0];
my $organism_name = $organism_name[scalar(@organism_name)-1];
my $seq = "PEPTIDE SEQUENCES";
my $ids = "";
my $out_name="tmp/shortened_seq_$organism_name";
my $counter = 0;
my $num=0;
open(my $raw_data, '<', $file_name) or die "Could not open	 '$file_name' $! \n";
open (OUT, ">$out_name$num") or die "Can't write to file '$out_name$num' [$!]\n";;
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
	if ($counter > 9999){
		$seq = substr $seq, 0, 69;
		print OUT "$seq\n";
		$num++;
		$counter = 0;
		close (OUT);
		open (OUT, ">$out_name$num");
	}
}
$seq = substr $seq, 0, 69;
print OUT "$seq\n";
close (OUT);

system("cat tmp/shortened_seq_$organism_name* > results/SP/ensembl_parsed_$organism_name.txt");

my $file = $out_name;

for (my $no = 0; $no <= $num; $no++) {
	system("./software/signalp-4.1/signalp -t euk -f summary -u 0.34 -U 0.34 $file$no");
	system("rm $file$no");
}
