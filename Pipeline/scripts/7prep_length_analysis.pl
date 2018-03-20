#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: perl prep2spl_length_analysis.pl .cfg extracted_sigp_all.out organism_ids_mapper.csv

EOF

my $file = shift @ARGV || die "sigp_ALL.out is missing";
my $ids_mapper = shift @ARGV || die "organism_ids_mapper is missing";;
my @orgs = @ARGV;
open (IN0, "+<", $ids_mapper);
my $j = 0;
my $i = 0;
my %HoA;
while (my $line1 = <IN0>) {
	chomp $line1;
	my @fields1 = split "," , $line1;
	for (my $i = 0; $i < scalar(@orgs); $i++){
		push @{$HoA{$orgs[$i]}}, $fields1[$i+1]; #Plus one because first column in the mapper is the ID
	}
}
close(IN0);

# Get letter codes of all analyzed organisms
my $remember = $i;
my @string;
my $multiplier = 1;
for(my $i = 0; $i < scalar(@orgs); $i++){
	my $rnd = int(rand($remember));
	if(${$HoA{$orgs[$i]}}[$rnd] eq "NULL"){
		until(${$HoA{$orgs[$i]}}[$rnd] ne "NULL"){
                        $multiplier++;
			$rnd = int((rand($remember))*$multiplier);
		}
		$string[$i] = ${$HoA{$orgs[$i]}}[$rnd];
	}
	else {
		$string[$i] = ${$HoA{$orgs[$i]}}[$rnd];
	}
	$string[$i] =~ s/[0-9]//g;
}
my @organisms = @string;

my %tmp_data;
open (my $sigp, '<', $file) or die "Could not open '$file' $! \n";
while (my $line1 = <$sigp>) {
	chomp $line1;
	if ($line1 =~ /^($string[0]\d+),(\w+\d+),\w+,\d+,\d+,\d+/) {
		my $index = first_index{/$1/} @{$HoA{$orgs[0]}};
		for(my $org = 0; $org < (scalar(@orgs)); $org++){
			open (my $sigp2, '<', $file) or die "Could not open '$file' $! \n";
			while (my $line2 = <$sigp2>) {
				chomp $line2;
				if ($line2 =~ /^($HoA{$orgs[$org]}[$index]),(\S+),(\w+),(\d+),(\d+)/){
					$tmp_data{$1}{sp_length} = $3;
					$tmp_data{$1}{lsaar_length} = $4;
				}
			}
			close($sigp2);
			unless($tmp_data{$HoA{$orgs[$org]}[$index]}{sp_length}){
				$tmp_data{$HoA{$orgs[$org]}[$index]}{sp_length} = "NA";
				$tmp_data{$HoA{$orgs[$org]}[$index]}{lsaar_length} = "NA";
			}
		}
		for(my $org = 0; $org < (scalar(@organisms)); $org++){
			unless($org == 0){
				print "," . $HoA{$orgs[$org]}[$index].",".$tmp_data{$HoA{$orgs[$org]}[$index]}{sp_length}.",".$tmp_data{$HoA{$orgs[$org]}[$index]}{lsaar_length};
			}
			else{
				print $HoA{$orgs[$org]}[$index].",".$tmp_data{$HoA{$orgs[$org]}[$index]}{sp_length}.",".$tmp_data{$HoA{$orgs[$org]}[$index]}{lsaar_length};
			}
		}
		print  "\n";
	}
}
close($sigp);
