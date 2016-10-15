#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: 03_SignalP_positives.pl organism 

EOF

my $sp_out = "$ARGV[0]_all.out" || die $usage;
my $id;
my $from;
my $to;
my $decision;
my $TrainSet;

open(IN,"$sp_out") || die("SignalP4 output file doesn't exist\n");
open(OUT1,">SignalP4_positives_$ARGV[0].out");

while(<IN>)
{
	if ($_ =~ /# SignalP-4.1\s+(\w+)\s+.*/)
	{
		$TrainSet = $1;
	}
	elsif ($_ =~ /\s+D\s+(\d+)\-(\d+)\s+\d+\.\d+\s+\d+\.\d+\s+(\w+).*/)
	{
		$from = $1;
		$to = $2;
		$decision = $3;
	}
	elsif ($_ =~ /Name=(\w+).*/)
	{
		$id = $1;
		if ($decision =~ /^YES/)
		{
			print OUT1 $id."\tsp4\t".$from."\t".$to."\t".$TrainSet."\n";
		}
	}
} 
print "\n$ARGV[0] SignalP positives saved to SignalP4_positives_$ARGV[0].out\n\n"; 

close(OUT1);
close(IN);
