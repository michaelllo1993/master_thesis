#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

# the order of organisms in tsv file and in cfg file must be the same
my $usage = <<EOF;

USAGE: 5create_ens_ids_mapper.pl <ortho_dataset> <organisms in the same order as in the dataset - first organism of interest and then others in alphabetical order>

EOF
my $dir = getcwd;

my %ids;
my $i = 0;
my $Ortho_data_set = shift(@ARGV);$Ortho_data_set = "$dir\/$Ortho_data_set";
system("sed -i 's/\"//g' $Ortho_data_set");

my @organisms = @ARGV;

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


# Creating Ensemble ids mapper
my @keys = keys %ids;
my $iter = 1;
print "0,";
print join(',', @organisms );
print "\n";
for (my $l = 1; $l <= scalar(@keys); $l++) {
	print "$iter";
	my $norg = scalar(@organisms);
	for (my $o = 0; $o < $norg; $o++) {
		print "," . $ids{$l}{$organisms[$o]};
	}
	print "\n"; $iter++;
}

