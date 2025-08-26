#!/usr/bin/env Rscript

# filter intermolecular/pairing records from IntaRNA results by removing records with low intermolecular-pairing base fraction
# last 20231101
# by b.p.f@qq.com

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='filter sense-antisense pairing records for IntaRNA result')
parser$add_argument('--inFile', type='character', required=TRUE,
                    help='input intermolecular interaction table from IntaRNA')
parser$add_argument('--interPairFrac', type='double', default=0.3,
                    help='intermolecular-pairing base fraction cutoff. default: >=0.3, range [0-1.0]')
parser$add_argument('--interPairLen', type='integer', default=15,
                    help='intermolecular-pairing base length cutoff. default: >=15, range [0-]')
parser$add_argument('--outFile', type='character', required=TRUE,
                    help='output intermolecular interaction table (with reduced columns) that removed records with low intercellular-pairing base fraction. detailed columns explained at https://github.com/BackofenLab/IntaRNA')
parser$add_argument('--outBed1', type='character', required=TRUE,
                    help='output bed that removed records with low intermolecular-pairing base fraction, length.')
parser$add_argument('--outBed2', type='character', required=TRUE,
                    help='output bed that removed records with low intermolecular-pairing base fraction, length.')
parser$add_argument('--outBed3', type='character', required=TRUE,
                    help='output bed that removed duplicated records.')
parser$add_argument('--outBed4', type='character', required=TRUE,
                    help='output bed that removed duplicated records.')                    

args <- parser$parse_args()

for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}
 


inFile="IntaRNA.txt"
interPairFrac <- 0.1
interPairLen <- 15
outFile="IntaRNA_filter.txt"
outBed1="IntaRNA_query.bed1"
outBed2="IntaRNA_target.bed2"
outBed3="IntaRNA_query_clean.bed1"
outBed4="IntaRNA_target_clean.bed2"


#df$hybridDPfull[1]

### E;E_hybrid;E_norm;
### Energy; Hybridization Energy; 



df <- read.csv(inFile,sep = ";")
message(nrow(df),"  ",ncol(df)) # 820 > 623
query <- unlist(sapply(strsplit(df$hybridDPfull, "&", fixed = T), "[", 1))
target <- unlist(sapply(strsplit(df$hybridDPfull, "&", fixed = T), "[", 2))

# for(i in 1:22){
#   print(lengths(regmatches(query[i], gregexpr("(", query[i], fixed = T))))
# }
#below seem to handle vec input
query.1 <- lengths(regmatches(query, gregexpr("(", query, fixed = T))) # query: (
query.2 <- lengths(regmatches(query, gregexpr(")", query, fixed = T))) # query: ), should be 0?
query.total <- nchar(query)
query.pair.frac <- (query.1-query.2)/query.total
#summary(query.pair.frac)

target.1 <- lengths(regmatches(target, gregexpr("(", target, fixed = T))) # target: (, should be 0?
target.2 <- lengths(regmatches(target, gregexpr(")", target, fixed = T))) # target: )
target.total <- nchar(target)
target.pair.frac <- (target.2-target.1)/target.total
#summary(target.pair.frac)

df$interPairLen1 <- query.1-query.2
df$interPairLen2 <- target.2-target.1
df$interPairFrac1 <- query.pair.frac
df$interPairFrac2 <- target.pair.frac


table(df$interPairFrac1 >= interPairFrac) 
table(df$interPairFrac2 >= interPairFrac)
table(df$interPairLen1 >= interPairLen)
table(df$interPairLen1 >= interPairLen)

sort(df$E)

# filter inter pair frac
filter <- (df$interPairFrac1 >= interPairFrac) & (df$interPairFrac2 >= interPairFrac) & 
  (df$interPairLen1 >= interPairLen) & (df$interPairLen2 >= interPairLen) # future consideration: df$E, df$Etotal

filter <- (df$E < -50)
table(filter)


filter <- (df$interPairLen1 >= interPairLen) & (df$interPairLen2 >= interPairLen) # future consideration: df$E, df$Etotal
table(filter)




df2 <- df[filter,c("id1","id2","E","Eall1","Eall2","Etotal","interPairFrac1","interPairFrac2","seq1","seq2","hybridDPfull")]

# write reduced table + Pairing Frac.
message(nrow(df)," > ",nrow(df2))


write.table(df2, paste0(outFile,""),quote=F,row.names=F,col.names=T,sep="\t")

## detailed columns explained at https://github.com/BackofenLab/IntaRNA
#E : overall interaction energy  = E_hybrid + ED1 + ED2
#Eall : ensemble energy of all considered interactions (-RT*log(Zall)) (Eall ~ E)
#Eall1 : ensemble energy of all considered intra-molecular structures of seq1 (given its accessibility constraints)
#Eall2 : ensemble energy of all considered intra-molecular structures of seq2 (given its accessibility constraints)
#Etotal : total energy of an interaction including the ensemble energies of intra-molecular structure formation (E+Eall1+Eall2)


# write paired gn bed6
coordTobed6 <- function(coord){
  bed1 <- data.frame(chr=unlist(sapply(strsplit(coord, ":", fixed = T), "[", 1)),
                    start=unlist(sapply(strsplit(coord, ":", fixed = T), "[", 2)),
                    strand=unlist(sapply(strsplit(coord, "_", fixed = T), "[", 2))
                    )
  bed1$end <- unlist(sapply(strsplit(bed1$start, "-", fixed = T), "[", 2))  
  bed1$start <- unlist(sapply(strsplit(bed1$start, "-", fixed = T), "[", 1))
  bed1$end <- unlist(sapply(strsplit(bed1$end, "_", fixed = T), "[", 1)) 
  bed1$score <- 1
  bed1$name <- paste0(bed1$chr,":",bed1$start,"-",bed1$end,"_",bed1$strand)
  bed1 <- bed1[,c("chr","start","end","name","score","strand")]
  #bed1[1:3,]
  return(bed1)
}

write.table(coordTobed6(df2$id1), paste0(outBed1,""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(coordTobed6(df2$id2), paste0(outBed2,""),quote=F,row.names=F,col.names=F,sep="\t")


 
