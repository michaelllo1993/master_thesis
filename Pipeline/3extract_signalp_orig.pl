#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;
USAGE: (name).pl (organism)_signalp_positives.out ensembl_parsed_(organism) config file
EOF

my $sp4_positives = "$ARGV[0]" || die $usage;
my $ensembl_parsed = "$ARGV[1]" || die $usage;

my $File = $ARGV[2];
open (my $CONFIG, $File);

my %User_Preferences;
while (my $line = <$CONFIG>) {
		chomp $line;
		if ($line =~ /^(\w+)\s*=\s*(\w+)/){
	    my $var = $1;
			my $value = $2;
	    $User_Preferences{$var} = $value;
		}
}
close ($CONFIG);

my $tag = $User_Preferences{SAAR};
my $tag_threshold = $User_Preferences{SAAR_threshold};

my %signalp=();

open (IN, "$sp4_positives") || die ("File missing.");

while(<IN>)
{
	chomp;									# delete the sign of the end of the line (like <ENTER>)
	my @entry = split("\t", $_);						# divide each line (wherever you find a tab) into separate words
	$signalp{$entry[0]}{cl_site} = $entry[3];				# assign second word from the B split (protid = key) to a value (cl_site)
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
	if ($_ =~ /^\>(\w+)\|(\w+)/)
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

my $seq_cl="";
my @tab="";
my $pos;
my $i;

if($tag =~ "all"){
	for my $id (keys %fasta)
	{
		if($signalp{$id}{decision})							# if it is TRUE - now I use id and not temp_protid because temp_protid was only local above
		{										# for SigP4 it's always TRUE - I don't print seq that don't have sig peptide
			$seq_cl= substr($fasta{$id}{seq},0, $signalp{$id}{cl_site});		# from $seq take only the part from 0 to cl_site (= signal peptide seq)
			print $id."\t".$fasta{$id}{id_t}."\t".$fasta{$id}{seq}."\t".$signalp{$id}{cl_site}.",";
		}
	}
}
elsif($tag =~ /^\w{1}/){
	for my $id (keys %fasta)
	{
		if($signalp{$id}{decision})				# if it is TRUE - now I use id and not temp_protid because temp_protid was only local above
		{
			$seq_cl= substr($fasta{$id}{seq},0, $signalp{$id}{cl_site});
			if(@tab = $seq_cl =~ /(${tag}{${tag_threshold},$signalp{$id}{cl_site}})/g)						# globally search for the pattern (max length = cl_site pos) in each sp seq
			{
				print $id.",".$fasta{$id}{id_t}.",".$fasta{$id}{seq}.",".$signalp{$id}{cl_site}.",";	# if we want only sigP length print $seq_cl instead of $fasta{$id}{seq}
				$pos=0;
				for($i=0; $i < scalar(@tab); $i++)
				{
					print length($tab[$i]).",".index($seq_cl,$tab[$i],$pos).",";	# prints length of the L-repeat (length($tab[$i])) and its position (index($seq_cl,$tab[$i],$pos)) in the sequence
					$pos += (index($seq_cl,$tab[$i],$pos)+length($tab[$i]));
				}
				print "\n";
			}
		}
	}
}
