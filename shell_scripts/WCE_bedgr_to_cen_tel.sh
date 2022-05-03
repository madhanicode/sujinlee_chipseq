d=/home/manning/scripts/sujin_testing/bedgraphs_WCE
#point refdir at directory for tel_windows_sorted.bed and cen_windows_sorted.bed
refdir=/home/manning/scripts/sujin_testing
files=$(ls $d/*.bedgraph)
o=/home/manning/scripts/sujin_testing/temp_outputs

for f in ${files}; do
    #trim label from end of string
	b=${f%.bedgraph}
    #trim leading /../ from string
    b=${b##*/}

	bedtools map -o sum -c 4 -a ${refdir}/tel_windows_sorted.bed -b ${f} > ${o}/${b}_tel.bedgraph
	bedtools map -o sum -c 4 -a ${refdir}/cen_windows_sorted.bed -b ${f} > ${o}/${b}_cen.bedgraph
done
