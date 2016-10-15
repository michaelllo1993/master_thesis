#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils 'first_index';


my $usage = <<EOF;

USAGE: 05_signal_peptide_only.pl <human_and_orthologues_protein_ids.tsv> <organisms> <sigpL_ALL.out>

EOF

my %ids;
my $i = 0;
my @organisms = ("homo_sapiens","pan_troglodytes","mus_musculus","bos_taurus","gallus_gallus","xenopus_tropicalis");
my $Ortho_data_set = $ARGV[0] || die $usage;
open (IN, "$Ortho_data_set") || die("Orthologues ids file missing.");

while (my $line = <IN>) {
	chomp $line;
	my @fields = split "\t" , $line;
	if (($fields[0]) && ($fields[0] =~ /^ENSP/)) {
	$i++;
	$ids{$i}{$organisms[0]} = $fields[0];
	if (!($fields[1])) {
		$ids{$i}{$organisms[1]} = "NULL";}
	else { $ids{$i}{$organisms[1]} = $fields[1];}
	if (!($fields[2])) {
		$ids{$i}{$organisms[2]} = "NULL";}
	else { $ids{$i}{$organisms[2]} = $fields[2];}
	if (!($fields[3])) {
		$ids{$i}{$organisms[3]} = "NULL";}
	else { $ids{$i}{$organisms[3]} = $fields[3];}
	if (!($fields[4])) {
		$ids{$i}{$organisms[4]} = "NULL";}
	else { $ids{$i}{$organisms[4]} = $fields[4];}
	if (!($fields[5])) {
		$ids{$i}{$organisms[5]} = "NULL";}
	else { $ids{$i}{$organisms[5]} = $fields[5];}
	}
}

# Creating Ensemble ids mapper
open (OUT0, ">","Ensembl_ids_mapper.csv");
my $keys = keys %ids;
my $iter = 0;
for (my $l = 1; $l <= $keys; $l++) {
  print OUT0 "$iter";
  for (my $o = 0; $o <= 5; $o++) {
  print OUT0 "," .$ids{$l}{$organisms[$o]};
  } print OUT0 "\n"; $iter++;
}
close (OUT0);
print "\nEnsembl_ids_mapper.csv created in /bi/ala/work/mstolarczyk/working_directory\n";

# Improting Ensemble ids mapper to arrays
open (IN0, "+<", "/bi/ala/work/mstolarczyk/working_directory/Ensembl_ids_mapper.csv");
my $ii = 0;
my @homo_sapiens;
my @pan_troglodytes;
my @mus_musculus; 
my @bos_taurus; 
my @gallus_gallus; 
my @xenopus_tropicalis;

while (my $line1 = <IN0>) {
	chomp $line1;
	my @fields1 = split "," , $line1;
	$homo_sapiens[$ii] = $fields1[1]; 
	$pan_troglodytes[$ii] = $fields1[2]; 
	$mus_musculus[$ii] = $fields1[3]; 
	$bos_taurus[$ii] = $fields1[4]; 
	$gallus_gallus[$ii] = $fields1[5]; 
	$xenopus_tropicalis[$ii] = $fields1[6];
	$ii++; 
}
close(IN0);
	
my %fasta;
my $id="";
my $seq="";
#my $g=0;

for (my $j = 1; $j <= 6; $j++) {
	my $parsed_FILE = "sigp_ALL_$ARGV[$j].out" || die $usage;
	open (IN, "$parsed_FILE") || die ("Parsed file is missing.");

	while(<IN>) {	
		chomp;
		if ($_ =~ /^(\S+)(\s+)(\S+)(\s+)(\w+)(\s+)(\d+)/) {	 
		#$g++;
			$fasta{$ARGV[$j]}{$1}=$5;
		#print "\n$g";
		}  
	} 
	
close(IN);
print "\n\nseq: $fasta{pan_troglodytes}{ENSPTRP00000018010}\n\n";
}


my $sigpL_ALL_FILE = $ARGV[7] || die $usage;
open (IN, "$sigpL_ALL_FILE") || die ("File missing.");

my %sigp_ALL_L;

while(<IN>) {

	if ($_ =~ /^(.+)\s(.+)\s(\w+)\s(\d+)\s(\d+)\s(\d+)/) {
		$sigp_ALL_L{$1}{seq} = $3;				# where $1 is the {Ensemble}
	}
}


close (IN);
my $ens_id;
my $seq_sigp;
my @id_o;
for $ens_id (keys %sigp_ALL_L) {		
		
		if ($ens_id =~ /^ENSP00000\d{6}/){				
			
			open (OUT, ">","$ens_id\_L_ALL_temp.fasta") || die("File missing.");
			print OUT ">".$ens_id."\n".$sigp_ALL_L{$ens_id}{seq};
			close (OUT);
			my $l_seq_sigp = length($sigp_ALL_L{$ens_id}{seq});
			open (OUT2, ">","$ens_id\_ORTHO_L_ALL_temp.fasta") || die("File missing.");
			my $index = first_index {/$ens_id/} @homo_sapiens; 
			$id_o[1] = $pan_troglodytes[$index];
			$id_o[2] = $mus_musculus[$index];
			$id_o[3] = $bos_taurus[$index];
			$id_o[4] = $gallus_gallus[$index];
			$id_o[5] = $xenopus_tropicalis[$index];
			
			for (my $on = 1; $on <= 5; $on++) { 
			    if (exists($fasta{$organisms[$on]}{$id_o[$on]})) {
				  $seq_sigp = $fasta{$organisms[$on]}{$id_o[$on]};
				  print OUT2 ">".$id_o[$on]."\n".$seq_sigp."\n";
				  $seq_sigp = "";
				}	
			}
			close(OUT2);

		my $exit_stat = system ("needle -outfile temp_L_ALL_Hsap_Ortho.needle -gapopen 10.0 -gapextend 0.5 -aformat markx3 $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
		system ("cat temp_L_ALL_Hsap_Ortho.needle >> ALL_Hsap_Ortho_sp.needle");
		if ($exit_stat != 0)
		{
			print $ens_id."_L_ALL\t".$exit_stat."\n";
		}
		system ("rm $ens_id\_L_ALL_temp.fasta $ens_id\_ORTHO_L_ALL_temp.fasta");
		}
		
}
system ("rm ENSP00*");
print "\nAlignment results saved to: ALL_Hsap_Ortho.needle\n";