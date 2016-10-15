#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: run_signalp.pl organism no_of_files

EOF

my $file = "shortened_seq_$ARGV[0]" || die $usage;
my $no = $ARGV[1] || die $usage;

for (my $num = 0; $num <= $no; $num++) {
	system("/bi/home/mstolarczyk/Programs/signalp-4.1/signalp -t euk -f summary -u 0.34 -U 0.34 $num$file > $num$ARGV[0].out");
	print "$ARGV[0] $num done \n";
}

system("cat *$ARGV[0].out > $ARGV[0]_all.out");
system("rm *$ARGV[0].out");
system("cat *shortened_seq_$ARGV[0] > ensembl_parsed_$ARGV[0]");
system("rm *shortened_seq_$ARGV[0]");
print "\nResults saved in $ARGV[0]_all.out\n\n";