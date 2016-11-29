#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
use Data::Dumper;
use List::MoreUtils 'first_index';

# <organism_and_orthologues_protein_ids.tsv> must have the same column order as <organism_of_interest>+<orthologues>
my $usage = <<EOF;

USAGE: Needle_ALL_wrapper.pl <organism_and_orthologues_protein_ids.tsv> <sigpL_ALL.out> <organism_of_interest> <orthologues>

EOF
my $dir = getcwd;
print "\n\nStarting analysis in: " . $dir . "\n\n";

my %ids;
my $i = 0;
my $Ortho_data_set = shift(@ARGV);$Ortho_data_set = "$dir\/$Ortho_data_set";
system("sed -i 's/\"//g' $Ortho_data_set");
my $sigpL_ALL_FILE = shift(@ARGV);$sigpL_ALL_FILE = "$dir\/$sigpL_ALL_FILE";
system("sed -i 's/\"//g' $sigpL_ALL_FILE");
my @organisms;

while (@ARGV){
	my $element = shift(@ARGV);
	push @organisms, $element;
}
print "Organisms in the analysis:\n";
foreach my $n (@organisms) {
	print $n . "\n";
}

open (IN, "$Ortho_data_set") || die("Orthologues ids file missing.");
my $ncols = scalar(@organisms);
while (my $line = <IN>) {
	chomp $line;
	my @fields = split("\t",$line);
	if (exists($fields[0])){
		if ($fields[0] =~ /^ENS/) {
			$i++;
			for (my $column = 0; $column < $ncols; $column++){
				if (!($fields[$column])) {
					$ids{$i}{$organisms[$column]} = "NULL";
				}
				else {
					$ids{$i}{$organisms[$column]} = $fields[$column];
				}
			}
		}
	}
}


print "\n\n";
print($Ortho_data_set);
print "\n\n";

# Get letter codes of all analysed organisms
close(IN);
my $remember = $i;
my @string;
for(my $i = 0; $i < scalar(@organisms); $i++){
	print $i . "\n";
	my $rnd = int(rand($remember));
	if($ids{$rnd}{$organisms[$i]} eq "NULL"){
		until($ids{$rnd}{$organisms[$i]} ne "NULL"){
			print $ids{$rnd}{$organisms[$i]} . "\n";
			$rnd = int(rand($remember));
		}
		$string[$i] = $ids{$rnd}{$organisms[$i]};
	}
	else {
		$string[$i] = $ids{$rnd}{$organisms[$i]};
	}
	print "$string[$i]\n";
	$string[$i] =~ s/[0-9]//g;
	print "$string[$i]\n";
}


# Creating Ensemble ids mapper
open (OUT0, ">","Ensembl_ids_mapper.csv");
my @keys = keys %ids;

my $iter = 0;
for (my $l = 1; $l <= scalar(@keys); $l++) {
	print OUT0 "$iter";
	my $norg = scalar(@organisms);
	for (my $o = 0; $o < $norg; $o++) {
		print OUT0 "," . $ids{$l}{$organisms[$o]};
	}
	print OUT0 "\n"; $iter++;
}
close (OUT0);
print "\nEnsembl_ids_mapper.csv created in $dir/Ensembl_ids_mapper.csv\n";

# Improting Ensemble ids mapper to hash of arrays
open (IN0, "+<", "$dir/Ensembl_ids_mapper.csv");
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


# Reading FASTA files into hash
my %fasta;
my $id="";
my $seq="";

for (my $j = 0; $j < scalar(@organisms); $j++) {
	my $parsed_FILE = "ensembl_parsed_$organisms[$j]" || die $usage;
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

# Reading sigpL_ALL_FILE into hash
open (INs, "$sigpL_ALL_FILE") || die ("File missing.");

my %sigpL_ALL_L;



while (my $line = <INs>) {
	chomp $line;
	my @fields = split(",",$line);
	if($fields[0] =~ /^ENS/){
		$sigpL_ALL_L{$fields[0]}{seq} = $fields[2];
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
		$index = first_index {$_ eq $ens_id} @{$HoA{$organisms[0]}}; print $index . "\n";
		for (my $on = 1; $on < scalar(@organisms); $on++) { #starts from 1 because first oraganism ($organism[0]) is the organism of interest (from first column of .tsv file)
			$id_o = $HoA{$organisms[$on]}[$index];
			if(($id_o !~ /^NULL/) && exists($fasta{$organisms[$on]}{$id_o})) {
				$seq_sigp = $fasta{$organisms[$on]}{$id_o}; print "$fasta{$organisms[$on]}{$id_o}\n";
				print OUT2 ">".$id_o."\n".$seq_sigp."\n";
				$seq_sigp = "";
			}
		}
		close(OUT2);
		my $exit_stat = system ("needle -outfile temp_L_ALL_$organisms[0]_Ortho.needle -gapopen 10.0 -gapextend 0.5 -aformat markx3 $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
		system ("cat temp_L_ALL_$organisms[0]_Ortho.needle >> ALL_$organisms[0]_Ortho_all.needle");
		if ($exit_stat != 0) {
			print $ens_id."_L_ALL\t".$exit_stat."\n";
		}
		system ("rm $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
	}
}
system ("rm ENS*");
print "\nAlignment results saved to: $dir/ALL_$organisms[0]_Ortho.needle\n";

