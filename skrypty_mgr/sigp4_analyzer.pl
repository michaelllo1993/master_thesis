#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

my $usage = <<EOF;

USAGE: sigp4_analyzer.pl organism 

EOF

my $dir = getcwd;
my $sp_out = "$dir/$ARGV[0]_all.out" || die $usage;
my $id;
my $from;
my $to;
my $decision;
my $TrainSet;

open(IN,"$sp_out") || die("SignalP4 output file doesn't exist\n");
open(OUT1,">$dir/SignalP4_positives_$ARGV[0].out");

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
 

close(OUT1);
close(IN);

print "\nSignalP positives saved to $dir/SignalP4_positives_$ARGV[0].out. \nCalcualtions are continued ...\n\n";

my $sp4 = "$dir/SignalP4_positives_$ARGV[0].out" || die $usage;
my $ensembl_parsed = "$dir/ensembl_parsed_$ARGV[0]" || die $usage;


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
my $id1="";
my $id_t="";
my $seq="";

open (IN, "$ensembl_parsed") || die ("File missing.");

while(<IN>)
{
	chomp;
	if ($_ =~ /^\>(\w+)\|(\w+)/)
	{	 
		if($1 ne $id1)
		{
			if($id1 ne "")
			{
				$fasta{$id1}{seq}=$seq;
				$fasta{$id1}{id_t}=$id_t;				
			}  
			$id1 = $1;
			$id_t = $2;
			$seq="";
		}               
	}
	else
	{	
		$seq.=$_; 
	}
}

$fasta{$id1}{seq}=$seq;										# the last value will be added (because the last one is not printed in the loop anymore)
$fasta{$id1}{id_t}=$id_t;

close(IN);



open (OUT1, ">>$dir/sigp_ALL_$ARGV[0].out") || die ("File missing.");
open (OUT2, ">>$dir/sigpL_ALL_$ARGV[0].out") || die ("File missing.");

my $seq_cl="";
my @tab="";
my $pos;
my $i;

for my $id1 (keys %fasta)
{
	if($signalp{$id1}{decision})									# if it is TRUE - now I use id and not temp_protid because temp_protid was only local above
	{												# for SigP4 it's always TRUE - I don't print seq that don't have sig peptide
	$seq_cl= substr($fasta{$id1}{seq},0, $signalp{$id1}{cl_site});					# from $seq take only the part from 0 to cl_site (= signal peptide seq)
	print OUT1 $id1."\t".$fasta{$id1}{id_t}."\t".$fasta{$id1}{seq}."\t".$signalp{$id1}{cl_site}."\n";
	if(@tab = $seq_cl =~ /(L{5,$signalp{$id1}{cl_site}})/g)						# globaly search for the pattern (max length = cl_site pos) in each sp seq
	{
		print OUT2 $id1."\t".$fasta{$id1}{id_t}."\t".$fasta{$id1}{seq}."\t".$signalp{$id1}{cl_site}."\t";	# if we want only sigP length print $seq_cl instead of $fasta{$id}{seq}
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

print "\nResults saved to $dir/sigp_ALL_$ARGV[0].out and $dir/sigpL_ALL_$ARGV[0].out\n\n"; 