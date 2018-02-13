#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

my $usage = <<EOF;

USAGE: perl prep2comp.pl ALL_organism_Ortho_(all/sp).needle <ortho_dataset.tsv> <organism_of_interest> <other_organisms>

EOF

my %ids;
my $i = 0;
my $dir = getcwd;
my @organisms;
my $file = shift(@ARGV) || die "ALL_Hsap_Ortho.needle is missing";
my $Ortho_data_set = shift(@ARGV);$Ortho_data_set = "$dir\/$Ortho_data_set";
system("sed -i 's/\"//g' $Ortho_data_set");

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

open (my $needle, '<', $file) or die "Could not open '$file' $! \n";

open (OUT, ">","prep2comp_$organisms[0]_all.csv");

while (<$needle>) {
	chomp;
	if (/^>($string[0]\d+)/) {
		print OUT "\n".$1.",";
	}
	for(my $i=1; $i < scalar(@organisms); $i++){
		if (/^>($string[$i]\d+)/) {
			print OUT ",".$1.",";
		}
	}
	if ((/^(-+)/)||(/^(\w+)/)) {
		print OUT "".$_;
	} 	
}
close(OUT);
close($file);

print "\n\nResult saved in $dir/prep2comp_$organisms[0]_all.csv\n\n";
