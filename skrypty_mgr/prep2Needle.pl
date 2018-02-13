#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: prep2Needle_development.pl input_file_name

EOF

my $input_file = "$ARGV[0].txt";
open(my $file, '<', $input_file) or die "Could not open '$input_file' $! \n";
open (OUT, ">$ARGV[0].out");

while (my $line = <$file>){
	chomp $line;
	if($line =~ /^>(ENS\S+)\|(ENS\S+)/){
		print OUT "\n" . $1 . "," . $2 . ",";
	}
	else{
		print OUT $line;
	}
}
close(OUT);
close($file);