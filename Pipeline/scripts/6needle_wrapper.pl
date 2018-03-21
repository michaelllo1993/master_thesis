#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: 6needle_wrapper.pl <extracted_sigp_organism.out> <organism_ids_mapper.csv> <codes_mapper> <organisms in the right order> <region of intrest (SP/WHOLE)>

EOF
my $dir = getcwd;

my $sigpL_ALL_FILE = shift @ARGV;
system("sed -i 's/\"//g' $sigpL_ALL_FILE");
my $ids_mapper = shift @ARGV;
my $codes_mapper = shift @ARGV;
my $region = shift @ARGV;
my @organisms = @ARGV;

my %ids;
my $i = 0;

# Importing Ensemble ids mapper to hash of arrays
open (IN0, "+<", $ids_mapper);
my $j = 0;
my %HoA;
while (my $line1 = <IN0>) {
	chomp $line1;
	my @fields1 = split "," , $line1;
	for (my $i = 0; $i < scalar(@organisms); $i++){
		push @{$HoA{$organisms[$i]}}, $fields1[$i+1]; #Plus one because first column in the mapper is the ID
	}
}
close(IN0);


# Get letter codes of all analyzed organisms
open (IN, "+<", $codes_mapper);
my @string;
while (my $line = <IN>) {
	chomp $line;
	my @fields = split "," , $line;
		$string[$fields[0]] = $fields[1];
	}
}
close(IN);

# Reading FASTA files into hash
my %fasta;
my $id="";
my $seq="";

for (my $j = 0; $j < scalar(@organisms); $j++) {
	#Establishing the sequence files names: in case of SP alignment - ensembl_parsed_*, in case of WHOLE sequence alignment - *.txt input files
	my $prefix = "";
	my $suffix = "";
	if($region =~ /^SP$/){
		$prefix = "results/SP/ensembl_parsed_";
		$suffix = ".txt";
	} elsif ($region =~ /^WHOLE$/){
		$prefix = "data/proteinSequences/";
		$suffix = ".txt";
	}
	my $parsed_FILE = "$dir/$prefix$organisms[$j]$suffix" || die $usage;
	#fixing files - deleting lines with faulty characters i.e. "_" or entries without sequences
	system("sed -i '/\_/d' $parsed_FILE");
	system("sed -i '/>|/d' $parsed_FILE");
	open (IN, "$parsed_FILE") || die ("Parsed file is missing.");
	while(<IN>) {
		chomp;
		if ($_ =~ /^\>(\w+)\|/) {
			if($1 ne $id) {
				if($id ne "")
				{
					$fasta{$organisms[$j]}{$id}=$seq;
				}
				$id = $1;
				$seq="";
			}
		}
		else {
			if ($_ !~ /^PEPTIDE/) {
				$seq.=$_;
			}
		}
	}
	$fasta{$organisms[$j]}{$id}=$seq; # the last value will be added (because the last one is not printed in the loop anymore)
	close(IN);
}
# Same here, depending on the region of interest (SP or WHOLE) the variable $sigpL_ALL_FILE holds extracted sigp results or WHOLE AA sequences
# Reading sigpL_ALL_FILE or organism_of_interest.out into hash
open (INs, "$sigpL_ALL_FILE") || die ("File missing.");

my %sigpL_ALL_L;

while (my $line = <INs>) {
	chomp $line;
	my @fields = split(",",$line);
	if($fields[0]){
		if($fields[0] =~ /^ENS/){
			$sigpL_ALL_L{$fields[0]}{seq} = $fields[2];
		}
	}
}

close(INs);

my $seq_sigp;
my $id_o;
my $index;

foreach my $ens_id (sort(keys %sigpL_ALL_L)) {
	if ($ens_id =~ /^$string[0]\d/) {
		open (OUT, ">","$ens_id\_L_ALL_temp.fasta") || die("File missing.");
		print OUT ">".$ens_id."\n".$sigpL_ALL_L{$ens_id}{seq};
		close (OUT);
		my $l_seq_sigp = length($sigpL_ALL_L{$ens_id}{seq});
		open (OUT2, ">","$ens_id\_ORTHO_L_ALL_temp.fasta") || die("File missing.");
		$index = first_index {$_ eq $ens_id} @{$HoA{$organisms[0]}};
		for (my $on = 1; $on < scalar(@organisms); $on++) { #starts from 1 because first oraganism ($organism[0]) is the organism of interest (from first column of .tsv file)
			$id_o = $HoA{$organisms[$on]}[$index];
			if(($id_o !~ /^NULL/) && exists($fasta{$organisms[$on]}{$id_o})) {
				$seq_sigp = $fasta{$organisms[$on]}{$id_o};
				print OUT2 ">".$id_o."\n".$seq_sigp."\n";
				$seq_sigp = "";
			}
		}
		close(OUT2);
		my $exit_stat = system ("needle -outfile temp_L_ALL_$organisms[0]_Ortho.needle -gapopen 10.0 -gapextend 0.5 -aformat markx3 $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
		system ("cat temp_L_ALL_$organisms[0]_Ortho.needle");
		# if ($exit_stat != 0) {
		# 	print $ens_id."_L_ALL\t".$exit_stat."\n";
		# }
		system ("rm $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
	}
}
system ("rm ENS*");
system ("rm temp*");
