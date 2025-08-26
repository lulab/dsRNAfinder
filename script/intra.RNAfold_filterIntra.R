#!/usr/bin/env Rscript

# filter intramolecular/pairing records from RNAfold results by removing records with low intramolecular-pairing base fraction
# last 20231101
# by b.p.f@qq.com
options(stringsAsFactors = F)
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='filter intra-stable records for RNAfold-permutation result')
parser$add_argument('--inFdr', type='character', required=TRUE,
                    help='input RNAfold-permutation FDR csv table from rnafold_dinushuffle_parallel.py')
parser$add_argument('--inBed', type='character', required=TRUE,
                    help='input bed file used for filtering')
parser$add_argument('--intraPairFrac', type='double', default=0.3,
                    help='intramolecular-pairing base fraction cutoff. default: >=0.3, range [0-1.0]')
parser$add_argument('--intraPairLen', type='integer', default=15,
                    help='intramolecular-pairing base length cutoff. default: >=15, range [0-]')
parser$add_argument('--intraMfeFdr', type='double', default=0.1,
                    help='intramolecular-stable MFE FDR. default: <=0.1, range [0-1.0]')
parser$add_argument('--adjMFE', type='double', default=-0.35,
                    help='intramolecular-stable length-adjusted MFE. default: <=-0.35, range [-00-+00]. ref: 2023,Bioinfo,dsRID')
parser$add_argument('--outBed', type='character', required=TRUE,
                    help='output bed that removed records with low intramolecular-pairing base fraction, length, or high FDR.')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}



fdr <- read.csv(inFdr,sep = ",",header = F)
fdr$V1 <- gsub(">","",fdr$V1,fixed = T)
fdr$V1 <- gsub("(+)","",fdr$V1,fixed = T)
fdr$V1 <- gsub("(-)","",fdr$V1,fixed = T)
#table(fdr$V2<=intraFdr) # ~50%

bed <- read.table(inBed,sep = "\t",header = F)
if(ncol(bed)==6){
  message("input bed is bed6")
}else if(ncol(bed)==12){
  message("input bed is bed12")
}

#below seem to handle vec input
query.1 <- lengths(regmatches(fdr$V4, gregexpr("(", fdr$V4, fixed = T))) # query: (
query.2 <- lengths(regmatches(fdr$V4, gregexpr(")", fdr$V4, fixed = T))) # query: ), should be 0?
query.total <- nchar(fdr$V4)
all(query.1==query.2)
query.pair.frac <- (query.1+query.2)/query.total

fdr$intraPairLen <- query.1+query.2
fdr$intraPairFrac <- query.pair.frac
hist(fdr$intraPairLen)
hist(fdr$intraPairFrac) 
str(fdr)

fdr$adj.MFE <- fdr$V3/nchar(fdr$V4) # ref: 2023,Bioinfo,dsRID
hist(fdr$adj.MFE) # adj.MFE
hist(fdr$V2)  # intraMfeFdr

table(fdr$intraPairFrac >= intraPairFrac)
table(fdr$intraPairLen >= intraPairLen)
table(fdr$V2 <= intraMfeFdr)
table(fdr$adj.MFE <= adjMFE)


fdr[1:4,]
fdr2[fdr2$V1=="chr11:65439375-65442875_+",]

# ————————————————————————————
# filter inter pair frac
filter <- (fdr$intraPairFrac >= intraPairFrac) & 
  (fdr$intraPairLen >= intraPairLen) & 
  (fdr$adj.MFE <= adjMFE) 
  #(fdr$V2 <= intraMfeFdr)   ### 不做intraMFEfdr<0.1 过滤掉80%的dsRNA
# table(filter)

fdr2 <- fdr[filter,]
dim(fdr2)

# write reduced table + Pairing Frac.
write.table(bed[bed$V4 %in% fdr2$V1,], paste0(outBed,""),quote=F,row.names=F,col.names=F,sep="\t")
