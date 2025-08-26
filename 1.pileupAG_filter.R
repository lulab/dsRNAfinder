#!/usr/bin/env Rscript

options(stringsAsFactors = F)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

parser <- ArgumentParser(description='filter Tab from pileup A-to-G')
parser$add_argument('--inTab', type='character',
                    help='input pileup table path')
parser$add_argument('--outTab', type='character',
                    help='outnput pileup table path')
parser$add_argument('--cov', type='integer', default=1,
                    help='Coverage filter in record >= Coverage, default: 1')
parser$add_argument('--edited', type='integer', default=1,
                    help='Edited reads filter in record >= edited, default: 1')
parser$add_argument('--snp', type='logical', default=TRUE,
                    help='SNP filter, default: TRUE')

#cov = 3
#edited = 1


args <- parser$parse_args()
inTab <- args$inTab
outTab <- args$outTab
dbSNP <- args$snp
cov <- args$cov
edited <- args$edited

# run filter
df <- read.table(inTab, header = F, check.names = F, sep = "\t", stringsAsFactors = F)
colnames(df) <- c("Chrom","Position","Strand","Coverage","Editedreads","Editlevel",
                "Type","dbSNP","Repeat","Fun_V34","Gene_V34","AA_V34","Fun_Ref","Exonic_Ref","AA_Ref")
df2 <- df %>%
      filter(dbSNP=="-") %>% 
      filter(Coverage>=cov) %>%
      filter(Editedreads>=edited) %>%
      # filter(Editlevel > 0.05) %>%
      mutate(Strand=ifelse(Strand=="-",0,1)) %>%
      transmute(Region=Chrom,
                Position=Position,
                Reference="AG",
                Strand=Strand,
                Coverage=Coverage,
                Editedreads=Editedreads,
                Editlevel=Editlevel,
                Type=Type,
                Repeat=Repeat,
                Gene_func=Fun_V34,
                Gene=Gene_V34)

message("reads before: ",nrow(df),"\nreads filtered: ",nrow(df2))

write.table(df2, file = outTab,quote = F,sep = "\t",row.names = F, col.names = T)

