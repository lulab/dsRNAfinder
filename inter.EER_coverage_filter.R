options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser(description='filter Tab from REDItools')
parser$add_argument('--ratio',  type='double', default=0.5,
                    help='coverage ratio')
parser$add_argument('--freq', type='double', default=0.5,
                    help='sample frequency')
parser$add_argument('--cov', type='integer', default=3,
                    help='coverage reads')
parser$add_argument('--input', type='character', default="./",
                    help='dsRIP coverage')
parser$add_argument('--output', type='character', default="./",
                    help='output file')

args <- parser$parse_args()

ratio = args$ratio
freq = args$freq
cov = args$cov
sampleid <- args$sampleid
input = args$input
output = args$output
gnFa = args$gnFa

bed_pos="inter_IntaRNA/intersection_overlap.+.bed"
bed_neg="inter_IntaRNA/intersection_overlap.-.bed"
bed_clean="inter_IntaRNA/intersection_clean.txt"

output2="inter_IntaRNA/intersection_overlap_covfil.bed"
output3="inter_IntaRNA/intersection_overlap_covmat.txt"
output4="inter_IntaRNA/intersection_clean_covfil.txt"

{
  bed <- read.table(bedpos,sep="\t",header=F,check.names=F)
  bed[1:4,]
  colnames(bed) <- c("chr","start","end","name","cov","strand")
  bed$strand = "."
  dim(bed)
  
  ids = read.table(sampleid,sep="\t",header=F,check.names=F)[,1]
  tmp1 <- read.table(paste0(input,"/",ids[1],".+"),sep="\t",header=F)
  mat <- data.frame(inter_dsRNA=tmp1$name)
  for (id in ids){
    print(id)
    tmp1 <- read.table(paste0(input,"/",id,".+"),sep="\t",header=F)
    tmp2 <- read.table(paste0(input,"/",id,".-"),sep="\t",header=F)
    colnames(tmp1) <- c("chr","start","end","name","cov","strand","coverage","cov_bp","ref_bp","cov_ratio")
    colnames(tmp2) <- c("chr","start","end","name","cov","strand","coverage","cov_bp","ref_bp","cov_ratio")
    
    table(tmp1$name==tmp2$name)
    
    ## interRNA 正负链的intersection区域都>0.8
    tmp1 <- tmp1 %>% mutate(filter_cov=ifelse(cov_ratio >= ratio, coverage, 0))
    tmp2 <- tmp2 %>% mutate(filter_cov=ifelse(cov_ratio >= ratio, coverage, 0))
    #tmp <- data.frame(a=ifelse(tmp1$filter_cov >= cov & tmp2$filter_cov >= cov, tmp1$filter_cov, 0),
    #                  b=ifelse(tmp1$filter_cov >= cov & tmp2$filter_cov >= cov, tmp2$filter_cov, 0))
    tmp <- data.frame(a=tmp1$filter_cov, b=tmp2$filter_cov)
    tmp5 <- data.frame(a=apply(tmp,1,min))
    colnames(tmp5) = id
    
    mat <- cbind(mat, tmp5)
    
  }
  
  mat <- mat %>% column_to_rownames("inter_dsRNA")
  dim(mat)
  mat[1:2,1:2]
  
  # coverage > 2, min
  keep = (rowSums(mat >= cov) > (ncol(mat) * freq)) # 30%的样本/6个样本中coverage大于>3
  table(keep)
  mat.fil <- mat[keep,]
  bed.fil <- bed[keep,]
  
  bed.fil[1:4,]
  # 820
  inter = read.table(bed_clean,sep="\t",header=T,check.names=F)
  inter.fil = inter[keep,]
  table(inter.fil$V14==bed.fil$name)
  
  
  write.table(bed.fil, output2, sep="\t",quote=F,row.names=F,col.names=F)
  #sort -k1,1V -k2,2n -k3,3n
  write.table(mat.fil, output3, sep="\t",quote=F,row.names=T,col.names=T)
  write.table(inter.fil, output4, sep="\t",quote=F,row.names=F,col.names=F)
  
  
}



#——————————————————————————
# get pair fasta

# use coverage filter bed ———— 820

outBed=output4
outBed = read.table(outBed, sep="\t",header=F,check.names=F)
outBed[1:4,]
{
  gnFasta=gnFa
  tmp.bed <- outBed[,7:12]
  colnames(tmp.bed) <- colnames(outBed)[1:6]
  outBed2 <- as.data.frame( rbind(outBed[,1:6], tmp.bed) )
  rm(tmp.bed)
  outBed2[1:4,]
  dim(outBed2)
  # 820*2 = 1640
  ######outBed2 <- outBed2[duplicated(outBed2),] # 筛选出能发生不同RNA正反链互补的RNA
  
  bedtoolsr::bt.getfasta(nameOnly = T, s = T, split = T,
                         fi = gnFasta,
                         output = paste0(output,"/tmp.fa"),
                         bed = outBed2)
  inFa <- seqinr::read.fasta(paste0(output,"/tmp.fa"))
  #file.remove(paste0(output,"/tmp.fa"))
  
  names(inFa) <- gsub("(+)","",names(inFa),fixed = T)
  names(inFa) <- gsub("(-)","",names(inFa),fixed = T)
  
  queryIdx <- "V4"
  queryFa <- inFa[outBed[[queryIdx]]] # inFa[outBed$V4]
  
  targetIdx <- "V10"
  targetFa <- inFa[outBed[[targetIdx]]]
  
  message("query fasta (+):",length(queryFa))
  message("target fasta (-):",length(targetFa))
  
  i=28
  #16,17,18,21,20,24,25,28
  #for (i in 1:length(queryFa)){
  
  seqinr::write.fasta(sequences = queryFa[outBed[[queryIdx]][i]], names = outBed[[queryIdx]][i],
                      file.out = paste0(output,"/inter_IntaRNA/fa/query_",i,".fa"),
                      open = "w", nbchar = 3000000000)
  seqinr::write.fasta(sequences = targetFa[outBed[[targetIdx]][i]], names = outBed[[targetIdx]][i],
                      file.out = paste0(output,"/inter_IntaRNA/fa/target_",i,".fa"),
                      open = "w", nbchar = 3000000000)
  #}
  
}

