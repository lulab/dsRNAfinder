



indir=$1
outdir=$2
file=$3
cov_num=$4


#cov_num=1 # 1/800=0.0013 # 3/800=0.0038
chrSize=$5

strand="+"
cat ${indir}/*/${file} \
            | bedtools sort \
            | bedtools genomecov -i - -strand ${strand} -g $chrSize -bg \
            | awk -v s=${strand} 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,"X",$4,s}}' \
            | awk -v c=${cov_num} '$5 >= c' \
            | sort -k1,1 -k2,2n | bedtools merge -s -c 2,3,5,6 -o collapse,collapse,collapse,collapse \
            | awk 'BEGIN{{OFS="\t";FS="\t"}}
            {{split($4,a,/,/); split($5,b,/,/); split($6,c,/,/); split($7,d,/,/);
            cov=0.0;for(i=1;i<=length(a);i++){{cov+=c[i]*(b[i]-a[i]);}}
            cov /= $3-$2;
            print $1,$2,$3,$1":"$2"-"$3"_"d[1],cov,d[1]
            }}' > ${outdir}/${file}.consensus.${strand}
strand="-"
cat ${indir}/*/${file} \
            | bedtools sort \
            | bedtools genomecov -i - -strand ${strand} -g $chrSize -bg \
            | awk -v s=${strand} 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,"X",$4,s}}' \
            | awk -v c=${cov_num} '$5 >= c' \
            | sort -k1,1 -k2,2n | bedtools merge -s -c 2,3,5,6 -o collapse,collapse,collapse,collapse \
            | awk 'BEGIN{{OFS="\t";FS="\t"}}
            {{split($4,a,/,/); split($5,b,/,/); split($6,c,/,/); split($7,d,/,/);
            cov=0.0;for(i=1;i<=length(a);i++){{cov+=c[i]*(b[i]-a[i]);}}
            cov /= $3-$2;
            print $1,$2,$3,$1":"$2"-"$3"_"d[1],cov,d[1]
            }}' > ${outdir}/${file}.consensus.${strand}
#0.0012,+: 788
#0.0012,-: 825
cat ${outdir}/${file}.consensus.[+-] | sort -k1,1 -k2,2n > ${outdir}/${file}
rm ${outdir}/${file}.consensus.[+-]





