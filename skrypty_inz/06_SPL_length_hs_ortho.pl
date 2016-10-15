#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: perl 07_SPL_length_hs_ortho.pl  sigp(L)_ALL.out

EOF

# Improting Ensemble ids mapper to arrays
open (IN, "+<", "/bi/ala/work/mstolarczyk/working_directory/Ensembl_ids_mapper.csv");
my $i = 0;
my @homo_sapiens;
my @pan_troglodytes;
my @mus_musculus; 
my @bos_taurus; 
my @gallus_gallus; 
my @xenopus_tropicalis;

while (my $line1 = <IN>) {
	chomp $line1;
	my @fields = split "," , $line1;
	$homo_sapiens[$i] = $fields[1]; 
	$pan_troglodytes[$i] = $fields[2]; 
	$mus_musculus[$i] = $fields[3]; 
	$bos_taurus[$i] = $fields[4]; 
	$gallus_gallus[$i] = $fields[5]; 
	$xenopus_tropicalis[$i] = $fields[6];
	$i++; 
}
close(IN);

my $file = $ARGV[0] || die "sigp_ALL.out is missing";
open (my $sigp, '<', $file) or die "Could not open '$file' $! \n";
if ($ARGV[0] =~ /sigpL_ALL.out/) {
	open (OUT, ">","SPL_length_Hsap_Ortho.csv");
	my $ii = 0;
	while (<$sigp>) {
		my $pt = "";
		my $mm = ""; 
		my $bt = ""; 
		my $gg = ""; 
		my $xt = "";
		my $ens_id = "";
		my $index = "";
		chomp;
		if (/^(ENSP0\d+)\s(ENST\d+)\s(\w+)\s(\d+)\s(\d+)/) {
			$ii++;
			my $index = first_index{/$1/} @homo_sapiens; print"$ii\n";
			my $pt = $pan_troglodytes[$index];
			my $mm = $mus_musculus[$index];
			my $bt = $bos_taurus[$index];
			my $gg = $gallus_gallus[$index];
			my $xt = $xenopus_tropicalis[$index];
			my $ptl = "NA"; my $ptll = "NA";
			my $mml = "NA"; my $mmll = "NA";
			my $btl = "NA"; my $btll = "NA";
			my $ggl = "NA"; my $ggll = "NA";
			my $xtl = "NA"; my $xtll = "NA";
			print OUT "\n".$1.",".$4.",".$5;
			open (my $sigp, '<', $file) or die "Could not open '$file' $! \n";
			while (my $line2 = <$sigp>) {
				chomp $line2;
				if ($line2 =~ /^($pt)\s(\S+)\s(\w+)\s(\d+)\s(\d+)/){
					$ptl = $4; $ptll = $5;
				}
				if ($line2 =~ /^($mm)\s(\S+)\s(\w+)\s(\d+)\s(\d+)/) {
					$mml = $4; $mmll = $5;
				}
				if ($line2 =~ /^($bt)\s(\S+)\s(\w+)\s(\d+)\s(\d+)/) {
					$btl = $4; $btll = $5;
				}
				if ($line2 =~ /(^$gg)\s(\S+)\s(\w+)\s(\d+)\s(\d+)/) {
					$ggl = $4; $ggll = $5;
				}
				if ($line2 =~ /(^$xt)\s(\S+)\s(\w+)\s(\d+)\s(\d+)/) {
					$xtl = $4; $xtll = $5;
				}
				
			} print OUT ",".$pt.",".$ptl.",".$ptll.",".$mm.",".$mml.",".$mmll.",".$bt.",".$btl.",".$btll.",".$gg.",".$ggl.",".$ggll.",".$xt.",".$xtl.",".$xtll;
		}
	}
close(OUT); 
}

else {
	open (OUT, ">","SP_length_Hsap_Ortho.csv");
	my $ii = 0;
	while (<$sigp>) {
		my $pt = "";
		my $mm = ""; 
		my $bt = ""; 
		my $gg = ""; 
		my $xt = "";
		my $ens_id = "";
		my $index = "";
		chomp;
		if (/^(ENSP0\d+)\s(ENST\d+)\s(\w+)\s(\d+)/) {
			$ii++;
			print OUT "\n".$1.",".$4;
			my $index = first_index{/$1/} @homo_sapiens; print"$ii\n";
			my $pt = $pan_troglodytes[$index];
			my $mm = $mus_musculus[$index];
			my $bt = $bos_taurus[$index];
			my $gg = $gallus_gallus[$index];
			my $xt = $xenopus_tropicalis[$index];
			my $ptl = "NA";
			my $mml = "NA";
			my $btl = "NA";
			my $ggl = "NA";
			my $xtl = "NA";
			open (my $sigp, '<', $file) or die "Could not open '$file' $! \n";
			while (my $line2 = <$sigp>) {
				chomp $line2;
				if ($line2 =~ /^($pt)\s(\S+)\s(\w+)\s(\d+)/){
					$ptl = $4;
				}
				if ($line2 =~ /^($mm)\s(\S+)\s(\w+)\s(\d+)/) {
					$mml = $4;
				}
				if ($line2 =~ /^($bt)\s(\S+)\s(\w+)\s(\d+)/) {
					$btl = $4;
				}
				if ($line2 =~ /(^$gg)\s(\S+)\s(\w+)\s(\d+)/) {
					$ggl = $4;
				}
				if ($line2 =~ /(^$xt)\s(\S+)\s(\w+)\s(\d+)/) {
					$xtl = $4;
				}
				
			} print OUT ",".$pt.",".$ptl.",".$mm.",".$mml.",".$bt.",".$btl.",".$gg.",".$ggl.",".$xt.",".$xtl;
	}
}
	close(OUT);
}


if ($ARGV[0] =~ /sigpL_ALL.out/) {
	print "\n\nResult saved to: SPL_length_Hsap_Ortho.csv\n\n";
}
else {
	print "\n\nResult saved to: SP_length_Hsap_Ortho.csv\n\n";
}



