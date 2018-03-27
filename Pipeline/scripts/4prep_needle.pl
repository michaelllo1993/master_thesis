#!/usr/bin/perl
use warnings;
use strict;

my $input_file = $ARGV[0];
open(my $file, '<', $input_file) or die "Could not open '$input_file' $! \n";

while (my $line = <$file>){
	chomp $line;
	if($line =~ /^>(ENS\S+)\|(ENS\S+)/){
		print "\n" . $1 . "," . $2 . ",";
	}
	else{
		print $line;
	}
}
close($file);
