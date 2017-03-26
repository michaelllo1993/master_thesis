#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: perl prep2spl_length_analysis.pl Ensembl_ids_mapper.csv sigp(L)_ALL.out organism_of_interest other_organisms (Ensembl exclusive code)

EOF

# Improting Ensemble ids mapper to arrays
my $ens_ids_mapper = shift(@ARGV) || die "ids mapper file is missing";
my $file = shift(@ARGV) || die "sigp_ALL.out is missing";
my $i = 0;
my @organisms;

while (@ARGV){
	my $element = shift(@ARGV);
	push @organisms, $element;
}

my %mapper;
open (my $mapper_file, '<', $ens_ids_mapper) or die "Could not open '$ens_ids_mapper' $! \n";
while (my $line = <$mapper_file>) {
	chomp $line; 
	my @fields = split(',',$line);
	for (my $org = 0; $org < (scalar(@organisms)); $org++){
		$mapper{$organisms[$org]}[$i] = $fields[$org+1] # in $fields[$org+1] the index is inceremented because first column in .csv is ID.
	}
	$i++; 
}
close($mapper_file);

my %tmp_data;
open (my $sigp, '<', $file) or die "Could not open '$file' $! \n";
open (OUT, ">","spl_length_analysis_prepped_file_$organisms[0]_Ortho.csv");
while (my $line1 = <$sigp>) {
	chomp $line1;
	if ($line1 =~ /^($organisms[0]P\d+),($organisms[0]T\d+),(\w+),(\d+),(\d+),(\d+)/) {
		my $index = first_index{/$1/} @{$mapper{$organisms[0]}};
		for(my $org = 0; $org < (scalar(@organisms)); $org++){
			open (my $sigp2, '<', $file) or die "Could not open '$file' $! \n";
			while (my $line2 = <$sigp2>) {
				chomp $line2;
				if ($line2 =~ /^($mapper{$organisms[$org]}[$index]),(\S+),(\w+),(\d+),(\d+)/){
					$tmp_data{$1}{sp_length} = $3;
					$tmp_data{$1}{lsaar_length} = $4;
				}
			}
			close($sigp2);	
			unless($tmp_data{$mapper{$organisms[$org]}[$index]}{sp_length}){
				$tmp_data{$mapper{$organisms[$org]}[$index]}{sp_length} = "NA";
				$tmp_data{$mapper{$organisms[$org]}[$index]}{lsaar_length} = "NA";
			}
		}
		for(my $org = 0; $org < (scalar(@organisms)); $org++){
			unless($org == 0){
				print OUT "," . $mapper{$organisms[$org]}[$index].",".$tmp_data{$mapper{$organisms[$org]}[$index]}{sp_length}.",".$tmp_data{$mapper{$organisms[$org]}[$index]}{lsaar_length};
			}
			else{
				print OUT $mapper{$organisms[$org]}[$index].",".$tmp_data{$mapper{$organisms[$org]}[$index]}{sp_length}.",".$tmp_data{$mapper{$organisms[$org]}[$index]}{lsaar_length};
			}
		}
		print OUT "\n";
	}
}
close(OUT);
close($sigp);