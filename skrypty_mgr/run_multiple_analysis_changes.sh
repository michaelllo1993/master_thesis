#! /bin/bash

# USAGE: bash run_multiple_analysis_changes.sh infile <eg. prep2comp_homo_sapiens_all.csv> organisms_exclusive_codes <eg. ENSPTR ENSMMUP ENSGGOP ENSMUSP ENSRNOP ENSBTAP ENSGALP ENSXETP>
codes=("$@")
echo ${codes[@]}
file=${codes[0]}
echo $file
unset codes[0]
echo ${codes[@]}
for code in ${codes[@]}
do
	echo "Running $code"
	Rscript reverse_changes_analysis.R $file $code
done  
