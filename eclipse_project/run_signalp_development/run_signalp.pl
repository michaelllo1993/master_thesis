#!/usr/bin/perl
use warnings;
use strict;
#TEST IT!
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

print "\nNumber of files to process: $num \n\n";


my $file = "shortened_seq_$ARGV[0]" || die $usage;


for (my $no = 0; $no <= $num; $no++) {
        system("/home/mstolarczyk/Programy/signalp/signalp-4.1/signalp -t euk -f summary -u 0.34 -U 0.34 $no$file > $no$ARGV[0].out");
        print "$ARGV[0] $num done \n";
}

system("cat *$ARGV[0].out > $ARGV[0]_all.out");
system("rm *$ARGV[0].out");
system("cat *shortened_seq_$ARGV[0] > ensembl_parsed_$ARGV[0]");
system("rm *shortened_seq_$ARGV[0]");
print "\nResults saved in $ARGV[0]_all.out\n\n";


