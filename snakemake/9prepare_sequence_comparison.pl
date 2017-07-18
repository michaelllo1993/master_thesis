#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

my $usage = <<EOF;

USAGE: perl 9prepare_sequence_comparison.pl organism_prepapre_sequence_comparison.cfg organism_needle_out organism_ids_mapper.csv

EOF

my %ids;
my $i = 0;
my $dir = getcwd;
my @organisms;
my $file = shift(@ARGV) || die "ALL_Organism_Ortho.needle is missing";
my $Ortho_data_set = shift(@ARGV);$Ortho_data_set = "$dir\/$Ortho_data_set";
system("sed -i 's/\"//g' $Ortho_data_set");

#Reading config file
my $File = shift @ARGV;
open (my $CONFIG, $File);

my %User_Preferences;
while (my $line = <$CONFIG>) {
		chomp $line;
		if ($line =~ /^(\w+)\s*=\s*(\w+\.*\w*)/){
	    my $var = $1;
			my $value = $2;
	    $User_Preferences{$var} = $value;
		}
}
close ($CONFIG);

my @organisms;
foreach my $key (sort(keys %User_Preferences)) {
  if ($key =~ /^org/){
    push @organisms, $User_Preferences{$key};
  }
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


# Get letter codes of all analyzed organisms
close(IN);
my $remember = $i;
my @string;
for(my $i = 0; $i < scalar(@organisms); $i++){
	my $rnd = int(rand($remember));
	if($ids{$rnd}{$organisms[$i]} eq "NULL"){
		until($ids{$rnd}{$organisms[$i]} ne "NULL"){
			$rnd = int(rand($remember));
		}
		$string[$i] = $ids{$rnd}{$organisms[$i]};
	}
	else {
		$string[$i] = $ids{$rnd}{$organisms[$i]};
	}
	$string[$i] =~ s/[0-9]//g;
}













my $ids_mapper = shift @ARGV;
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
my $remember = $i;
my @string;
for(my $i = 0; $i < scalar(@organisms); $i++){
	my $rnd = int(rand($remember));
	if(${$HoA{$organisms[$i]}}[$rnd] eq "NULL"){
		until(${$HoA{$organisms[$i]}}[$rnd] ne "NULL"){
			$rnd = int((rand($remember))*10);
		}
		$string[$i] = ${$HoA{$organisms[$i]}}[$rnd];
	}
	else {
		$string[$i] = ${$HoA{$organisms[$i]}}[$rnd];
	}
	$string[$i] =~ s/[0-9]//g;
}


open (my $needle, '<', $file) or die "Could not open '$file' $! \n";

while (<$needle>) {
	chomp;
	if (/^>($string[0]\d+)/) {
		print "\n".$1.",";
	}
	for(my $i=1; $i < scalar(@organisms); $i++){
		if (/^>($string[$i]\d+)/) {
			print ",".$1.",";
		}
	}
	if ((/^(-+)/)||(/^(\w+)/)) {
		print "".$_;
	} 	
}
close($file);
