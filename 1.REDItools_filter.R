#!/usr/bin/env Rscript

options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser(description='filter Tab from REDItools')
parser$add_argument('--inTab', type='character',
                    help='input REDItools table path')
parser$add_argument('--outTab', type='character',
                    help='outnput REDItools table path')
parser$add_argument('--covQ30', type='integer', default=1,
                    help='filter in record >= Coverage-q30, default: 1')
parser$add_argument('--MeanQ', type='double', default=0.0,
                    help='filter in record >= MeanQ, default: 0')
parser$add_argument('--AllSubs', type='logical', default=TRUE,
                    help='whether to filter in AllSubs within c("AG","CT"), this will be as "-" if Pvalue>1*e-04.  default: TRUE (only keep A-->I/G and C-->U)')
parser$add_argument('--Frequency', type='double', default=0.0,
                    help='filter in record >= Frequency, default: 0')
#parser$add_argument('--Pvalue', type='logical', default=FALSE
#                    help='filter in record <= Pvalue, default: 0.05, (REDItoolsDenovo only)')

args <- parser$parse_args()
inTab <- args$inTab
outTab <- args$outTab
covQ30 <- args$covQ30
MeanQ <- args$MeanQ
AllSubs <- args$AllSubs
Frequency <- args$Frequency
Pvalue <- FALSE


# run filter
df <- read.table(inTab, header = T, check.names = F, sep = "\t", stringsAsFactors = F)
#colnames(df)[5] <- c("covQ30") # c("Coverage-q30","MeanQ","BaseCount[A,C,G,T]","AllSubs","Frequency","Pvalue")
if(ncol(df)==10){
  message("input is denovo mode, 10 cols")
}else if(ncol(df)==9){
  message("input is dnarna mode, 9 cols")
}else if(ncol(df)==14){
  message("input is denove? mode, 14 cols")
}else{
  message("input is ? mode, ",ncol(df)," cols")
}

# summary 
hist(df$`Coverage-q30`,xlim = c(0,199),breaks = 1000)
summary(df$MeanQ)
table(df$AllSubs %in% c("AG","CT")) 
table(grepl("AG|CT", df$AllSubs, perl = T))
summary(df$Frequency)


if(AllSubs){
  df2 <- df[(df$`Coverage-q30`>=covQ30) &
              (df$MeanQ>=MeanQ) &
              grepl("AG|CT", df$AllSubs, perl = T) &
              (df$Frequency>=Frequency) #(df$Pvalue<=Pvalue)
            ,]
}else{
  df2 <- df[(df$`Coverage-q30`>=covQ30) &
              (df$MeanQ>=MeanQ) &
              (df$Frequency>=Frequency) #(df$Pvalue<=Pvalue)
            ,]
}

message(nrow(df))
message(nrow(df2))

write.table(df2, file = outTab,quote = F,sep = "\t",row.names = F, col.names = T)

