#!/usr/bin/perl
use warnings;
use strict;

my $sp_out = "$ARGV[0]";
my $id_prot;
my $id_cdna;
my $from;
my $to;
my $decision;
my $TrainSet;

open(IN,"$sp_out") || die("SignalP4 output file doesn't exist\n");

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
	elsif ($_ =~ /Name=(\w+)\|(\w+).*/)
	{
		$id_prot = $1;
		$id_cdna = $2;
		if ($decision =~ /^YES/)
		{
			print $id_prot."\t".$id_cdna."\tsp4\t".$from."\t".$to."\t".$TrainSet."\n";
		}
	}
} 

close(IN);

