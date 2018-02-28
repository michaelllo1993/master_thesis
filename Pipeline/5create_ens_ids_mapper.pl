#!/usr/bin/perl
use warnings;
use strict;
use Cwd;

# the order of organisms in tsv file and in cfg file must be the same
my $usage = <<EOF;

USAGE: 5create_ens_ids_mapper.pl .cfg 

EOF
my $dir = getcwd;

#Reading config file

my $File = shift @ARGV;
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

my @organisms;
foreach my $key (sort(keys %User_Preferences)) {
  if ($key =~ /^org/){
    push @organisms, $User_Preferences{$key};
  }
}

my %ids;
my $i = 0;
my $Ortho_data_set = shift(@ARGV);$Ortho_data_set = "$dir\/$Ortho_data_set";
system("sed -i 's/\"//g' $Ortho_data_set");

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
my $iter = 0;
for (my $l = 1; $l <= scalar(@keys); $l++) {
	print "$iter";
	my $norg = scalar(@organisms);
	for (my $o = 0; $o < $norg; $o++) {
		print "," . $ids{$l}{$organisms[$o]};
	}
	print "\n"; $iter++;
}

