#! /bin/bash

# This script takes prep2comp file as the input and outputs the tmp_org_code file with sequences to be compared in changes_analysis.R.
# Basically it filtrates only the specified organism from the fourth column of the file
file=$1
org_code=$2
OLDIFS=$IFS
IFS=","
while read org1 seq1 org2 seq2
do
	if [[ "$org2" == *"$org_code"* ]]
	then
		echo $org1,$seq1,$org2,$seq2 >> tmp_$org_code.csv
	fi
done < "$file"
IFS=$OLDIFS

