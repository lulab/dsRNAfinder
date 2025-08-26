#!/bin/bash

#rediTab=$1 # REDItools tbl results as input 
#chrSize=$2 # chrsize
#gnFa=$3 # gn/tx fa
#binBed=$4 # gn/tx bin bed
#EERNum=$5 # cutoff of of EER num in $binSize bin
#binSize=$6 # length/size of bin
#extSize=$7 # extended length each side of merged EER bins, length/size of extended flank, usually half of binSize
#minLen=$8 # min length/size cutoff of merged,extended EER cluster region
#maxLen=$9 # max length/size cutoff of merged,extended EER cluster region
#tmpDir=${10} # tmp dir

rediTab=$1

chrSize=$2
gnFa=$3
binSize=50 
binBed=$4
tmpDir=$5

EERNum=$6
mergeSize=$7
extSize=$8
minLen=$9
maxLen=$10


echo "rediTab: ${rediTab}"
echo "chrSize: ${chrSize}"
echo "gnFa: ${gnFa}"
echo "binBed: ${binBed}"
echo "EERNum: ${EERNum}"
echo "binSize: ${binSize}"
echo "mergeSize: ${mergeSize}"
echo "extSize: ${extSize}"
echo "minLen: ${minLen}"
echo "maxLen: ${maxLen}"
echo "tmpDir: ${tmpDir}"

#while (( $# > 0 ))    # or [ $# -gt 0 ]
#do
#	    echo "$1"
#	        shift
#	done
	
source /BioII/lulab_b/baopengfei/mambaforge/bin/activate REDItools
export PATH=/BioII/lulab_b/zhanqing/cfRNAseq/tools/bedops/bin:$PATH
echo "start getEERcluster at `date`"

## filter EER-cluster bins
### count&filter
/BioII/lulab_b/zhanqing/cfRNAseq/tools/REDItools/bin/TableToGFF.py \
	-i ${rediTab} \
	-b 300000 -T ${tmpDir} \
	-o ${rediTab}.gff
gff2bed < ${rediTab}.gff \
	> ${rediTab}.bed
cut -f1-6 ${rediTab}.bed > ${rediTab}.bed6


cat ${rediTab}.bed6 \
	| sort -k1,1 -k2,2n \
	| bedtools coverage -s -sorted -counts \
		-a ${binBed} -b - \
	| awk 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,$4,$7,$6}}' \
	| awk -v num=${EERNum} '($5>=num) {print $0}' \
	> ${rediTab}_sort.bed6.count
#-g ${chrSize} \

## (optional: append editing sites position to Tab/bed)
#bedtools intersect ... 
#bedtools merge ... -c 2 -o collapse


## merge&extend EER-cluster bins
#note that extended region has potential editing sites, we neglect for current version
bedtools merge -s -d ${mergeSize} -c 4,5,6 -o distinct,sum,distinct -i ${rediTab}_sort.bed6.count \
	> ${rediTab}_sort.bed6.count.merge

# bedtools slop, increase the size of each feature by defined bases (start-extsize, end+extsize)
bedtools slop -s -l ${extSize} -r ${extSize} -g ${chrSize} \
	-i ${rediTab}_sort.bed6.count.merge \
	> ${rediTab}_sort.bed6.count.merge.ext

## filter length & add name
awk -v minL=${minLen} -v maxL=${maxLen} '(($3-$2)>=minL) && (($3-$2)<=maxL) {print $0}' ${rediTab}_sort.bed6.count.merge.ext \
	| awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3"_"$6 "\t" $5 "\t" $6}' \
	> ${rediTab}_sort.bed6.count.merge.ext.filterLen
echo "records num: `wc -l ${rediTab}_sort.bed6.count.merge.ext.filterLen`"

rm ${rediTab}.bed ${rediTab}.bed6 ${rediTab}_sort.bed6.count
echo "end getEERcluster at `date`"

