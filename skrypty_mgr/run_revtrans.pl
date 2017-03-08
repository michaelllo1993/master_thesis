#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
use Data::Dumper;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: run_revtrans.pl <organism of interest Latin name> <organism exclusive ENS code>

EOF
my $code = $ARGV[1];
my $revtrans = "/home/mstolarczyk/Programy/RevTrans/RevTrans-1.4/revtrans.py";
my $dir = getcwd;
my $output_file = "$dir/$ARGV[0]_revtrans.out";
my $output_file_csv1 = "$dir/$ARGV[0]_revtrans1.csv";
my $output_file_csv2 = "$dir/$ARGV[0]_revtrans2.csv";
my $output_file_csv = "$dir/$ARGV[0]_revtrans.csv";

my $prep2comp_file = "$dir/prep2comp_$ARGV[0]_all.csv";
print "\n\nStarting analysis in: " . $dir . "\n\n";

system("cat prep2revtrans* > prep2revtrans_tmp.csv");

my $prep2revtrans_file = "$dir/prep2revtrans_tmp.csv";
my @prep2revtrans_ids;
my @prep2revtrans_seqs;
open (IN, $prep2revtrans_file) || die("$prep2revtrans_file is missing");
my $i = 0;
while (my $line = <IN>){
	chomp $line;
	my @fields = split ("," , $line);
	$prep2revtrans_ids[$i] = $fields[0];
	$prep2revtrans_seqs[$i] = $fields[1];
	$i++;
}
close(IN);

open (IN1, $prep2comp_file) || die("$prep2comp_file is missing");
my $cnt = 1;
while (my $line = <IN1>){
	chomp $line;
	if ($line =~ /^EN/){
		my @fields = split(",", $line);
		my $prot_tmp_file = "prot_tmp_$fields[0]";
		open (OUT, ">", $prot_tmp_file) || die("$prot_tmp_file is missing");
		my $cDNA_tmp_file = "cDNA_tmp_$fields[0]";
		open (OUT1,">", $cDNA_tmp_file) || die("$cDNA_tmp_file is missing");
		my $index_interest = first_index {$_ eq $fields[0]} @prep2revtrans_ids;
		my $index_other = first_index {$_ eq $fields[2]} @prep2revtrans_ids;
		print OUT ">$fields[0]\n$fields[1]\n>$fields[2]\n$fields[3]";
		print OUT1 ">$fields[0]\n$prep2revtrans_seqs[$index_interest]\n>$fields[2]\n$prep2revtrans_seqs[$index_other]";
		close(OUT);
		close(OUT1);
		my $revtrans_output_tmp = "tmp_$fields[0].out";		
		system("$revtrans $cDNA_tmp_file $prot_tmp_file > $revtrans_output_tmp");
		system("cat $revtrans_output_tmp >> $output_file");
		system("echo @ >> $output_file");
		$cnt++; print"$cnt \n";
		system("rm *tmp*");
	}
}
close(IN1);
print("DONE.\n");

#Parsing the output; writing to the .csv file
#Firstly, only human? ids and sequences are written
my $status;
open (IN2, $output_file) || die("$output_file is missing");
open (OUT2, ">", $output_file_csv1) || die("$output_file_csv1 is missing");
print OUT2 "$ARGV[0]_id" . "," . "$ARGV[0]_sequence" . "\n";
while (my $line = <IN2>){
	chomp $line;
	if ($line =~ /^>($code\d+)/){
		print OUT2 $1 . ",";
		$status = 1;
	}
	elsif ((($line =~ /^[A-Z]/) || ($line =~ /^-/)) && ($status == 1)){
		print OUT2 $line;
		$status = 1;
	}
	elsif ($line =~ /^@/){
		print OUT2 "\n";
		$status = 0;
	}
	else {
		$status = 0;
	}
}
close(IN2);
close(OUT2);

#Secondly, other organisms' ids and sequences are written
open (IN2, $output_file) || die("$output_file is missing");
open (OUT2, ">", $output_file_csv2) || die("$output_file_csv2 is missing");
print OUT2 "other_id" . "," . "other_sequence" . "\n";
while (my $line = <IN2>){
	chomp $line;
	if (($line !~ /^>($code\d+)/) && ($line =~ /^>(\w+\d+)/)) {
		print OUT2 $1 . ",";
		$status = 1;
	}
	elsif ((($line =~ /^[A-Z]/) || ($line =~ /^-/)) && ($status == 1)){
		print OUT2 $line;
		$status = 1;
	}
	elsif ($line =~ /^@/){
		print OUT2 "\n";
		$status = 0;
	}
	else {
		$status = 0;
	}
}
close(IN2);
close(OUT2);

#Finally, the two files are merged together
system("paste -d ',' $output_file_csv1 $output_file_csv2 > $output_file_csv");
system("sed -i '/^,/ d' $output_file_csv");
system("rm $output_file_csv1 $output_file_csv2 $output_file");
