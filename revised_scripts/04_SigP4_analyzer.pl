#!/usr/bin/perl
use warnings;
use strict;


# Usage

my $usage = <<EOF;

USAGE: SigP_analyzer.pl organism 

EOF

my $sp4 = "SignalP4_positives_$ARGV[0].out" || die $usage;
my $ensembl_parsed = "shortened_seq_$ARGV[0]" || die $usage;


my %signalp=();

open (IN, "$sp4") || die ("File missing.");

while(<IN>)
{
	chomp;									# delete the sign of the end of the line (like <ENTER>)
	my @entry = split("\t", $_);						# divide each line (wherever you find a tab) into separate words
	$signalp{$entry[0]}{cl_site} = $entry[3];				# assign second word from the B split (protid = key) to a value (cl_site)  ==> mozna dodac wiecej zmiennych ktore nas interesuja
	$signalp{$entry[0]}{decision} = 1;
}
close(IN);


my %fasta=();
my $id="";
my $id_t="";
my $seq="";

open (IN, "$ensembl_parsed") || die ("File missing.");

while(<IN>)
{
	chomp;
	if ($_ =~ /^\>(\w+)\|(\w+)\|/)
	{	 
		if($1 ne $id)
		{
			if($id ne "")
			{
				$fasta{$id}{seq}=$seq;
				$fasta{$id}{id_t}=$id_t;				
			}  
			$id = $1;
			$id_t = $2;
			$seq="";
		}               
	}
	else
	{	
		$seq.=$_; 
	}
}

$fasta{$id}{seq}=$seq;										# the last value will be added (because the last one is not printed in the loop anymore)
$fasta{$id}{id_t}=$id_t;

close(IN);



open (OUT1, ">>sigp_ALL_$ARGV[0].out") || die ("File missing.");
open (OUT2, ">>sigpL_ALL_$ARGV[0].out") || die ("File missing.");

my $seq_cl="";
my @tab="";
my $pos;
my $i;

for my $id (keys %fasta)
{
	if($signalp{$id}{decision})									# if it is TRUE - now I use id and not temp_protid because temp_protid was only local above
	{												# for SigP4 it's always TRUE - I don't print seq that don't have sig peptide
	$seq_cl= substr($fasta{$id}{seq},0, $signalp{$id}{cl_site});					# from $seq take only the part from 0 to cl_site (= signal peptide seq)
	print OUT1 $id."\t".$fasta{$id}{id_t}."\t".$fasta{$id}{seq}."\t".$signalp{$id}{cl_site}."\n";
	if(@tab = $seq_cl =~ /(L{5,$signalp{$id}{cl_site}})/g)						# globaly search for the pattern (max length = cl_site pos) in each sp seq
	{
		print OUT2 $id."\t".$fasta{$id}{id_t}."\t".$fasta{$id}{seq}."\t".$signalp{$id}{cl_site}."\t";	# if we want only sigP length print $seq_cl instead of $fasta{$id}{seq}
		$pos=0;
		for($i=0; $i < scalar(@tab); $i++)
		{
			print OUT2 length($tab[$i])."\t".index($seq_cl,$tab[$i],$pos)."\t";	# prints length of the L-repeat (length($tab[$i])) and its position (index($seq_cl,$tab[$i],$pos)) in the sequence
			$pos += (index($seq_cl,$tab[$i],$pos)+length($tab[$i]));
		}
		print OUT2 "\n";
	}
	}
}

close (OUT1);
close (OUT2);


