#!/usr/bin/env Rscript

# filter sense-antisense pairing intervals and extract fasta for blast/RNAcofold/IntaRNA
# last 20231101
# by b.p.f@qq.com

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("seqinr"))


parser <- ArgumentParser(description='filter sense-antisense pairing intervals and extract fasta for blast/RNAcofold/IntaRNA')
parser$add_argument('--inBed', type='character', required=TRUE,
                    help='input genomic BED6/BED12 file path')
#parser$add_argument('--inFasta', type='character', #required=TRUE,
#                    help='input fasta file path, should contain all records in bed file')
#parser$add_argument('--gnFasta', type='character',
#                    help='genome fasta file path, required if no inFasta provided.')
parser$add_argument('--outDir', type='character', default='.',
                    help='output dir, default: .')
#parser$add_argument('--bedtoolsPath', type='character', 
#                    help='bedtools path')
# parser$add_argument('--cores', type="integer", default=4,
#                    help='multi-cores processing, default: half of auto detect')
args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}


# 20250218
options(bedtools.path = bedtoolsPath) 
suppressPackageStartupMessages(library("bedtoolsr"))

tmp <- read.table(inBed,header = F,sep = "\t")
if(ncol(tmp)==6){
  isBed6 <- TRUE 
  isBed12 <- FALSE  
  message("input is bed6")
}else if(ncol(tmp)==12){
  isBed6 <- FALSE  
  isBed12 <- TRUE  
  message("input is bed12")
}
tmp <- NULL


# get intersect sense-antisense records
outBed <- bedtoolsr::bt.intersect(S = T, wo = T, # sorted = T,
                        a = inBed,
                        b = inBed
                          ) # S: different strand
# Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported
# 一个bed文件内所有RNA 能正反互补（different strand）的情况
# 输出A的bed6和B的bed6

outBed_overlap <- bedtoolsr::bt.intersect(S = T, 
                        a = inBed,
                        b = inBed
                          ) # S: different strand
# 输出重叠区域的坐标, A的ID

table(outBed$V6==outBed_overlap$V6)
table(outBed$V4==outBed_overlap$V4)


outBed_clean = outBed %>% 
  filter(V6=="+") %>% 
  mutate(V14 = paste0(V4,"|",V10))

outBed_overlap = outBed_overlap %>%
  filter(V6=="+") %>%
  mutate(V4 = outBed_clean$V14)

outBed_overlap_neg = outBed_overlap %>%
  mutate(V6="-")

dim(outBed)
dim(outBed_clean)
dim(outBed_overlap)
outBed_clean[1:4,]
outBed_overlap[1:4,]


dim(outBed)
if(nrow(outBed)>0){
  dir.create(paste0(outDir,"/fa/"),showWarnings = T,recursive = T)
  
  ## write intersection/pairing table
  #bed12: 25 cols; bed6: 13 cols
  write.table(outBed, paste0(outDir,"/intersection.txt"), 
              quote=F,row.names=F,col.names=T,sep="\t")
  
  write.table(outBed_clean, paste0(outDir,"/intersection_clean.txt"), 
              quote=F,row.names=F,col.names=T,sep="\t")
  
  write.table(outBed_overlap, paste0(outDir,"/intersection_overlap.+.bed"), 
              quote=F,row.names=F,col.names=F,sep="\t")
  
  write.table(outBed_overlap_neg, paste0(outDir,"/intersection_overlap.-.bed"), 
              quote=F,row.names=F,col.names=F,sep="\t")
}

