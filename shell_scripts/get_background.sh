#!/bin/bash
# usage: get_background bedgraph_dir output_dir
# note: point refdir at directory containing "NOT_CEN_TEL_labeled.bed"

refdir=/home/manning/scripts/sujin_testing/scripts_dani
d=$1
o=$2
bgs=$(ls ${d}/*.bedgraph)

echo file average > ${o}/bg_averages.txt
for bg in ${bgs}; do
	raw=$(basename $bg)
	b=${raw%.bedgraph}
	# mh comment: how NOT_CEN_TEL_labeled.bed made unknown
	bedtools map -o sum -c 4 -a ${refdir}/NOT_CEN_TEL_labeled.bed -b $bg > ${o}${b}_not_in.bedgraph
	# cat the resulting bedgraph | print end - start in new column | print read depth / interval
	cat ${o}${b}_not_in.bedgraph | awk '{print $0"\t"$3-$2}' > ${o}${b}_not_in_averages.bedgraph
	awk '{sum1+=$5; sum2+=$6;} END {print FILENAME,sum1/sum2}' ${o}${b}_not_in_averages.bedgraph >> ${o}/bg_averages.txt
	#awk '{sum+=$7} END {print FILENAME,sum/NR}' ${o}${b}_not_in_averages.bedgraph

	#rm ${o}${b}_not_in_averages.bedgraph
	#rm ${o}${b}_not_in.bedgraph
done
