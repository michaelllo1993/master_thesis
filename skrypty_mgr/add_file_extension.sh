#!/bin/bash
ext=$1;
a=($(ls)); 
for filename in ${a[@]}
do 
	echo $filename
	new_filename=$filename".$ext"
	mv $filename $new_filename
done

