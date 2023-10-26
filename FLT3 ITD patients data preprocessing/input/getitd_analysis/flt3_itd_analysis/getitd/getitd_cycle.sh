#!/bin/bash

ls ../input/ | grep _r1.fastq | while read line
do

	patient=$(echo $line | cut -d '_' -f1)
	python3 getitd.py -reference ../my_anno_files/custom_ref1.txt -anno ../my_anno_files/my_anno1.txt -minscore_alignments 0.4 -require_indel_free_primers False -min_insert_seq_length 6 $patient ../input/$patient'_r1.fastq' ../input/$patient'_r2.fastq'


done
