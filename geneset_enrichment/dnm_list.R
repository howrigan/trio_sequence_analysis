## dnm_list.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2019


## ----- DESCRIPTION:

## Interactive Script to generate simple DNM lists for gene set enrichment analysis

## ----- NOTES:

## Does not require trio count as gene set enrichment is relative to total DNM count


## read in DNM calls
dnm <- read.table('../files/combined_cohorts_DNM_list.tsv',h=T,sep='\t',stringsAsFactors=F)

## remove DNMs in genes outside of CCDS database
dnm <- subset(dnm,!is.na(dnm$CCDS))

## subet DNMs
scz <- subset(dnm,dnm2$DISEASE=='SCZ')
con <- subset(dnm,dnm2$DISEASE=='SIB_CONTROL')

## write out simple DNM gene lists (gene symbol and annotation)
write.table(scz[,c('gene_symbol_used','annotation_used')],'SCZ_n2772_CCDS.dnm',col=F,row=F,quo=F,sep='\t')
write.table(con[,c('gene_symbol_used','annotation_used')],'CON_n2216_CCDS.dnm',col=F,row=F,quo=F,sep='\t')
