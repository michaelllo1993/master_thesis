#!/usr/bin/perl
use warnings;
use strict;

my $usage = <<EOF;

USAGE: perl 09_prep2comp.pl  ALL_Hsap_Ortho_(all/sp).needle

EOF

my $file = $ARGV[0] || die "ALL_Hsap_Ortho.needle is missing";
open (my $needle, '<', $file) or die "Could not open '$file' $! \n";

if ($ARGV[0] =~ /ALL_Hsap_Ortho_all.needle/) {
	open (OUT, ">","prep2comp_all.csv");
}
else {
	open (OUT, ">","prep2comp_sp.csv");
}
while (<$needle>) {
	chomp;
	if (/^>(ENSP0\d+)/) {
		print OUT "\n".$1.",";
	}
	if ((/^>(ENSMUSP\d+)/)||(/^>(ENSPTRP\d+)/)||(/^>(ENSGALP\d+)/)||(/^>(ENSBTAP\d+)/)||(/^>(ENSXETP\d+)/)) {
		print OUT ",".$1.",";
	}
	if ((/^(-+)/)||(/^(\w+)/)) {
		print OUT "".$_;
	} 	
}
close(OUT);
close($file);

if ($ARGV[0] =~ /ALL_Hsap_Ortho_all.needle/) {
	print "\n\nResult saved to: prep2comp_all.csv\n\n";
}
else {
	print "\n\nResult saved to: prep2comp_sp.csv\n\n";
}
