#!/bin/bash
array=(`ls *_ids_mapper.cfg`)
echo "${array[@]}"

len=${#array[*]}

i=0
while [ $i -lt $len ]; do
	echo "$i: ${array[$i]}"
	cp 
done

