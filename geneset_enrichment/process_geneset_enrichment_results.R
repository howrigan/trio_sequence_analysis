#!bin/Rscript

# process_geneset_enrichment_results.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: February 2019


## ----- DESCRIPTION:

## Command line script to convert the default DNM gene set enrichment results into multiple .tsv tables
## Orders by mutation model p-value
## Adds prettier column names

# USAGE:

# Rscript process_geneset_enrichment_results.R \
# [results file] \
# [output path] \
# [order by pval] \


# example usage:

# Rscript process_geneset_enrichment_results.R \
# SCZ_n2772_CCDS_coverageQC_candidate_genesets.tsv \
# results/SCZ_n2772_CCDS_coverageQC/ \
# TRUE


## ----- NOTES:

## requires the default annotations output from geneset_enrichment.R 





## read in command line argument
args <- commandArgs(TRUE)
infile <- args[1]
outpath <- args[2]
ordering <- args[3]

## check / create results directory
system(paste0('mkdir -p ',outpath))

input <- read.table(infile,h=T,sep='\t',stringsAsFactors=F)

nms <- c('gene set','genes listed','genes not tested','genes tested','percent genes tested',
        'SCZ DNMs tested','Control DNMs tested','SCZ DNMs overlapping','Control DNMs overlapping',
        'SCZ observed proportion','Model expectation proportion','Model enrichment','Model pval','Model lower 95% CI','Model upper 95% CI',
        'Control observed proportion','CaseControl enrichment','CaseControl pval','CaseControl lower 95% CI','CaseControl upper 95% CI')

all <- input[,c(1:5,6,7,16,17,27,26,28:31,33:37)]
if (ordering==T) { all <- all[order(all$all_pval),] }
names(all) <- nms
write.table(all,paste0(outpath,'all.tsv'),col=T,row=F,quo=F,sep='\t')

lofmis <- input[,c(1:5,8,9,18,19,39,38,40:43,45:49)]
if (ordering==T) { lofmis <- lofmis[order(lofmis$lofmis_pval),] }
names(lofmis) <- nms
write.table(lofmis,paste0(outpath,'lofmis.tsv'),col=T,row=F,quo=F,sep='\t')

lof <- input[,c(1:5,10,11,20,21,51,50,52:55,57:61)]
if (ordering==T) { lof <- lof[order(lof$lof_pval),] }
names(lof) <- nms
write.table(lof,paste0(outpath,'lof.tsv'),col=T,row=F,quo=F,sep='\t')

mis <- input[,c(1:5,12,13,22,23,63,62,64:67,69:73)]
if (ordering==T) { mis <- mis[order(mis$mis_pval),] }
names(mis) <- nms
write.table(mis,paste0(outpath,'mis.tsv'),col=T,row=F,quo=F,sep='\t')

syn <- input[,c(1:5,14,15,24,25,75,74,76:79,81:85)]
if (ordering==T) { syn <- syn[order(syn$syn_pval),] }
names(syn) <- nms
write.table(syn,paste0(outpath,'syn.tsv'),col=T,row=F,quo=F,sep='\t')


## END of process_geneset_enrichment_results.R
