# dsRNAfinder

![](./dsRNAfinder.webp)

## Preparation 

### environment

### reference
```shell
gnFa="hg38.fa"
chrSize="chrom.size"
binSize=50
binBed="hg38.bins.50.bed"
```


## RNA editing calling
### REDItools
```shell
# REDItools
python2.7 reditools.py \
      -f ${input} \
      -r hg38.fa \
      -S -s 1 -ss 5 -mrl 20 -q 10 -bq 20 -C -T 2 -os 5 \
      -o {output}

# samtools mpileup
perl adapted_Query_Editing_Level.GRCh37.20161110.pl ${input} ${output}
```

### Editing sites filter
```shell
# REDItools
Rscript 1.REDItools_filter.R \
			--inTab ${input} \
			--outTab ${output} \
			--covQ30 3 \
			--MeanQ 0.0 \
			--AllSubs TRUE 
# samtools mpileup
Rscript 1.pileupAG_filter.R \
					--inTab ${input} \
					--outTab ${output} \
					--snp TRUE \
					--cov 3 --edited 1
```

## EER cluster
### get EER cluster
```shell
EERnum=1
mergeSize=1000
extSize=25
minLen=50
maxLen=3000
bash 2.getEERcluster.sh ${input} \
			${chrSize} ${gnFa} ${binSize} ${binBed} ${tmpDir}
			${EERnum} ${mergeSize} ${extSize} ${minLen} ${maxLen}
```

### get consensus EER cluster
```shell
cov_num=2
cores=8
file="consensusEERcluster.bed"
bash 3.getEERconsensus.sh \
	EERcluster/ \
	EERconsensus/ \
	${file} ${cov_num} ${chrSize}
```

## Intra-dsRNA identification
### filter dsRIP coverage
```shell
# intra dsRNA
bed="consensusEERcluster.bed"
dsRIPid="dsRIP_id"
input="dsRIP_hg38_dedup.bedgraph"
bedtools coverage -a ${bed} -b ${input} -s > ${output.intra} 

# filter intra- dsRNA
Rscript 4.EER_coverage_filter.R \
					--ratio 0.5 \
					--freq 0.5 \
					--coc 3 \
					--output intra_RNAfold/consensusEER_covfil.bed
```
### RNA fold
```shell
bedtools getfasta -nameOnly -s \
			 -fi ${gnFa} \
       -bed intra_RNAfold/consensusEER_covfil.bed \
       > intra_RNAfold/consensusEER_covfil.bed.fa

shuffle_times=50
seed=1234
cores=16
python 5.rnafold_dinushuffle_parallel.py \
        intra_RNAfold/consensusEER_covfil.bed.fa \
        ${shuffle_times} ${seed} \
        intra_RNAfold/consensusEER_covfil.bed.fa.csv \
        ${cores} \
				> intra_RNAfold/consensusEER_covfil.bed.fa.log 2>&1
rm $output/intra_RNAfold/${file}.fa_perm
```
### RNA fold filter
```shell
6.RNAfold_filterIntra.R
```

## Inter-dsRNA identification
### get sense-antisense pairs
```shell
file="consensusEERcluster.bed"
Rscript 4.getSenseAntisensePair.R \
	--inBed ${file} \
	--outDir inter_IntaRNA
```

### filter reads coverage
```shell
# inter dsRNA
bed_pos="intersection_overlap.+.bed"
bed_neg="intersection_overlap.-.bed"
bedtools coverage -a ${bed_pos} -b ${input} -s > ${output.inter_pos} 
bedtools coverage -a ${bed_neg} -b ${input} -s > ${output.inter_neg}
      
# filter intra- and inter- dsRNA
Rscript 5.EER_cov_filter.R
```

### IntaRNA
```shell
cores=16
start=n1
end=n2
bash 8.IntaRNA.sh inter_IntaRNA $cores $n1 $n2
head -n 1 $output/result/IntaRNA_1.txt > $output/IntaRNA.txt
for i in `ls $output/result | grep -v MYPAIRMINE | grep -v log | cut -d "." -f 1 | cut -d "_" -f 2 | sort -n`
do
        echo $i
        tail -n +2 -q $output/result/IntaRNA_${i}.txt >> $output/IntaRNA.txt 
done
```

### IntaRNA filter
```shell
```






