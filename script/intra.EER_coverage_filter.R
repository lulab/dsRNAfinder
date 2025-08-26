options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser(description='filter Tab from REDItools')
parser$add_argument('--ratio',  type='double', default=0.5,
                    help='coverage ratio')
parser$add_argument('--freq', type='double', default=0.5,
                    help='sample frequency')
parser$add_argument('--cov', type='integer', default=3,
                    help='coverage reads')
parser$add_argument('--bed', type='character', default="consensusEERcluster.bed",
                    help='input files')
parser$add_argument('--input', type='character', default="./",
                    help='dsRIP coverage')
parser$add_argument('--output', type='character', default="consensusEER_covfil.bed",
                    help='output file')

args <- parser$parse_args()
ratio <- args$ratio
freq <- args$freq
cov <- args$cov
bed <- args$bed
input <- args$input
sampleid <- args$sampleid
output <- args$output

#ratio = 0.5
#freq = 0.5
#cov = 3
#output1 = paste0(dir,"intra_RNAfold/consensusEER_covfil.bed")



{
  bed <- read.table(bed,sep="\t",header=F,check.names=F)
  bed[1:4,]
  colnames(bed) <- c("chr","start","end","dsRNAid","cov","strand")
  bed <- bed %>% mutate(length=end-start)
  table(bed$length)
  dsrna1=bed
  dim(dsrna1)
  
  # 18 J2 samples
  ids = read.table(sampleid,sep="\t",header=F,check.names=F)[,1]
  tmp1 <- read.table(paste0(input,"/",ids[1],".EERcluster.coverage"),sep="\t",header=F)
  mat <- data.frame(EERcons=tmp1$name)
  for (id in ids){
    print(id)
    tmp1 <- read.table(paste0(input,"/",ids[1],".EERcluster.coverage"),sep="\t",header=F)
    colnames(tmp1) <- c("chr","start","end","name","cov","strand","coverage","cov_bp","ref_bp","cov_ratio")
    
    ## dsRIP J2 若cov < 0.5 则不算count数 
    tmp1[tmp1$cov_ratio < ratio,"coverage"] = 0
    tmp1. <- data.frame(a=tmp1$coverage)
    colnames(tmp1.) = id
    mat <- cbind(mat, tmp1.)
  }
  mat <- mat %>% column_to_rownames("EERcons")
  dim(mat)
  mat[1:2,1:2]
  keep = (rowSums(mat >= cov) > (ncol(mat) * freq))   # 50%的样本/9个样本中coverage大于>3
  table(keep)
  mat.fil <- mat[keep,]
  dim(mat.fil)
  

  
  bed2 = data.frame(chr=str_split(rownames(mat.fil),":",simplify=T)[,1],
                    start=str_split(rownames(mat.fil),":|-",simplify=T)[,2],
                    end=str_split(rownames(mat.fil),":|-|_",simplify=T)[,3],
                    strand=str_split(rownames(mat.fil),"_",simplify=T)[,2])
  bed2 = bed2 %>% transmute(chr=chr, start=start, end=end, name=rownames(mat.fil), cov=".", strand=strand)
  bed2 = bed2 %>% filter(chr %in% paste0("chr",c(1:22,"X","M")))
  dim(bed2)    
  bed2[12:43,]
  
  write.table(bed2,output,sep="\t",quote=F,row.names=F,col.names=F)          

}



