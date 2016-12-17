#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
use Data::Dumper;
use List::MoreUtils 'first_index';

my $usage = <<EOF;

USAGE: run_revtrans.pl <organism of interest Latin name>

EOF

my $revtrans = "/home/mstolarczyk/Programy/RevTrans/RevTrans-1.4/revtrans.py";
my $dir = getcwd;
my $output_file = "$dir/revtrans.out";
my $prep2comp_file = "$dir\/$ARGV[0]";
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
my $size = @prep2revtrans_ids;
print $size;
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
		system("echo @ >> revtrans.out");
		$cnt++; print"$cnt of $size\n";
	}
system("rm *tmp*");
}
close(IN1);
system("rm *tmp*");
