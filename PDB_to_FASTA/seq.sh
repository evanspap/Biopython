#!/bin/bash
for i in *.pdb
do
	name=$(basename -- "$i")
	python -W ignore ../../pdb-seq.py $i > ./SEQ/${name%.*}.seq
	
done
