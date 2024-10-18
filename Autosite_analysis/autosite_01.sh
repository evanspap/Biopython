#!/bin/bash
for i in ../../Autosite_2/*/*_fp_*.pdb
do
	filename=$(basename -- "$i")
	dir=$(dirname "$i")
	name="${filename%.*}"
	protein="${name%%_*}"
	pocketnt="${name: -3}"
	#echo $i
	#echo $dir
	#echo $filename
	#echo $name
	echo -n "$protein , $pocketnt , "
	echo
	#python -W ignore ../../pdb-seq.py $i > ./SEQ/${name%.*}.seq
	
done
