dir1=/home/manning/scripts/sujin_testing/temp_outputs
dir2=/home/manning/scripts/sujin_testing/bedgraph_cen_tel

files=$(ls $dir1/*.bedgraph)

for f in ${files}; do
    #trim leading /../ from string
    b=${f##*/}
    diff ${dir1}/${b} ${dir2}/${b}
    echo finished
done