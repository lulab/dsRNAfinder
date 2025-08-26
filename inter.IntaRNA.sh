

outDir=$1
cores=$2
pair_start=$3
pair_end=$4


condaDir="/BioII/lulab_b/baopengfei/mambaforge/bin"
source $condaDir/activate IntaRNA

#pair_num=`tail -n +2 ${outDir}/intersection.txt | wc -l | cut -d " " -f 1` # not count header
#echo "total intercellular pairs: ${pair_num}"


mkdir -p ${outDir}/result/MYPAIRMINE/

for i in `seq ${pair_start} ${pair_end}`
do
echo "Begin intercellular pairs $i in ${pair_start}:${pair_end} at `date`"
IntaRNA --acc=N --mode=H --model=X --outMode=C --outOverlap=N --outCsvCols '*' --threads ${cores} \
	--out=pMinE:${outDir}/result/MYPAIRMINE/MYPAIRMINE-t${i}q${i}.csv \
	-q ${outDir}/fa/query_${i}.fa \
	-t ${outDir}/fa/target_${i}.fa \
	> ${outDir}/result/IntaRNA_${i}.txt
done
echo "END IntaRNA"

# output
# --outMode：设置输出模式，如 C（CSV 格式）或 D（详细模式）。
# --outCsvCols：指定 CSV 输出文件的列

# 能量模型
# -m 或 --model：指定能量模型 
# --model=H, Note, due to the low run-time requirement of the heuristic prediction mode (--mode=H), heuristic IntaRNA interaction predictions are widely used to screen for interaction in a genome-wide scale. 
# --model=X, This default model of IntaRNA predicts the single-site interaction I with minimal free energy. 
# --model=M, If you are more interested in specific details of an interaction site or of two relatively short RNA molecules, you should investigate the exact prediction mode (--mode=M) providing the global minimum free energy interaction.
# https://github.com/BackofenLab/IntaRNA?tab=readme-ov-file#IntaRNA
# --seedBP：设置种子区域的最小碱基对数（默认 7）。
# --seedMaxE：设置种子区域的最大能量阈值。

# --acc=N, turn off accessibility consideration

#IntaRNA --acc=N --mode=H --model=X --outMode=D --outOverlap=N --threads ${cores} \
#	-q ${outDir}/fa/query_${i}.fa \
#	-t ${outDir}/fa/target_${i}.fa \
#	> ${outDir}/result/MYPAIRMINE/structure-t${i}q${i}.txt


#head -n 1 ${outDir}/result/IntaRNA_1.txt > ${outDir}/result/IntaRNA.txt; tail -n +2 -q ${outDir}/result/IntaRNA_*.txt >> ${outDir}/result/IntaRNA.txt # concate multi output with the same header
#rm ${outDir}/result/IntaRNA_*.txt
#support move than exterior-exterior interaction
#support multiple matrix output for visualization 




