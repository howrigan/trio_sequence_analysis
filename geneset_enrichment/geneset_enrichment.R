#!bin/Rscript

# geneset_enrichment.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: February 2019


## ----- DESCRIPTION:

## Command line script to generate DNM gene set enrichment results


# USAGE:

# Rscript geneset_enrichment.R \
# [affected DNM list] \
# [unaffected DNM list] \
# [mutation expection gene list] \
# [gene sets file] \
# [output results file] \
# [DNM overlap details path] \
# [annotation type] \
# [header status] \

## --- argument descriptions:

## affected DNM list: DNM list for affected samples from dnm_list.R containing gene symbol (column 1) and annotation (column 2)
## unaffected DNM list: DNM list for unaffected samples from dnm_list.R containing gene symbol (column 1) and annotation (column 2)
## mutation expection gene list: list of genes with mutation expectation for default annotations
## gene sets file: gene set list with gene symbols (column 1) and gene set name (column 2) 
## output results file: results filename
## DNM overlap details path: directory path for listing the overlapping genes in each gene set. Produces one file per gene set
## annotation type: use 'default_annotation' (5 tests: all,ptv,ptv+mis,mis,syn) or 'custom_annotation' (1 aff/unaff test)
## OPTIONS: [default_annotation] [custom_annotation]
## header status: does the gene set list have a column header or not? 
## OPTIONS: [header] [no_header]


# example usage:

# Rscript geneset_enrichment.R \
# SCZ_n2772_CCDS.dnm \
# CON_n2216_CCDS.dnm \
# ../files/gencode_pct75_gene17925.tsv \
# ../files/candidate_genesets.tsv \
# SCZ_n2772_CCDS_coverageQC_candidate_genesets.tsv \
# overlap/SCZ_n2772_CCDS_coverageQC/ \
# default_annotation \
# header


## ----- SCRIPT SET UP:

## PART 1: read in command line arguments
## PART 2: read in files / parse command line arguments
## PART 3: Running enrichment with default annotations (all_types==T)
## PART 4: Running enrichment with custom annotation (all_types==F)


## ----- NOTES:

## gene set enrichment does not test DNM rate, but DNM proportions, and thus doesn't require sample counts

## the script currently restricts to genes in the mutation expection gene list, even for custom annotations and case/control comparisons

## mutation expection gene list requires these four columns:
## 1) gene - gene symbol
## 2) p_all - all coding mutation probability
## 2) p_syn - synonymous mutation probability
## 2) p_mis - missense mutation probability
## 2) p_lof - LoF/PTV/LGD mutation probability (sum of nonsense, frameshift, and essential splice probabilities)








## PART 1: read in command line arguments
args <- commandArgs(TRUE)
dnmlist <- args[1] ## observed DNMs
unafflist <- args[2] ## observed DNMs
genelist <- args[3] ## gene file
setfile <- args[4] ## gene set directory
out <- args[5] ## output file name
overlap_path <- args[6] ## list output directory and name
annot_types <- args[7] ## run primary annotation types or only single type
setfile_header <- args[8] ## does the geneset file have a header

cat(paste('affected DNM file =',dnmlist),'\n')
cat(paste('unaffected DNM file =',unafflist),'\n')
cat(paste('mutation expection list =',genelist),'\n')
cat(paste('gene sets file =',setfile),'\n')
cat(paste('output results file =',out),'\n')
cat(paste('DNM overlap details path =',overlap_path),'\n')
cat(paste('running all annots =',annot_types),'\n')


## PART 2: read in files / parse command line arguments

## read in files
dnm <- read.table(dnmlist,stringsAsFactors=F)
unaff <- read.table(unafflist,stringsAsFactors=F)
genes <- read.table(genelist,h=T,stringsAsFactors=F)


## --- parsing command line input

## check / create overlap directory
system(paste0('mkdir -p ',overlap_path))

## selecting DNM annotation types
if (annot_types=='default_annotation') {all_types <- TRUE}
if (annot_types=='custom_annotation') {all_types <- FALSE}

## read in gene set file
if (setfile_header=='header') { setlist <- read.table(setfile,h=T,sep='\t',stringsAsFactors=F) }
if (setfile_header=='no_header') { setlist <- read.table(setfile,h=F,sep='\t',stringsAsFactors=F) }


## get gene sets
setlist_name <- names(table(setlist[,2]))


## PART 3: Running enrichment with default annotations (all_types==T)

if (all_types==T) {

## get DNM counts
all <- dnm$V1
lofmis <- subset(dnm$V1,dnm$V2=='ptv' | dnm$V2=='missense')
lof <- subset(dnm$V1,dnm$V2=='ptv')
mis <- subset(dnm$V1,dnm$V2=='missense')
syn <- subset(dnm$V1,dnm$V2=='synonymous')

all_dnm <- rep(length(dnm$V1),length(setlist_name))
lofmis_dnm <- rep(sum(dnm$V2=='ptv' | dnm$V2=='missense'),length(setlist_name))
lof_dnm <- rep(sum(dnm$V2=='ptv'),length(setlist_name))
mis_dnm <- rep(sum(dnm$V2=='missense'),length(setlist_name))
syn_dnm <- rep(sum(dnm$V2=='synonymous'),length(setlist_name))

## --- Gene matching
genes_listed <- NA
genes_missed <- NA
genes_tested <- NA
genes_pct <- NA
genenames_missed <- NA

## --- Binomal testing
all_expected <- NA; all_actual <- NA; all_enrichment <- NA; all_pval <- NA; all_count <- NA; all_low95 <- NA; all_hi95 <- NA
lofmis_expected <- NA; lofmis_actual <- NA; lofmis_enrichment <- NA; lofmis_pval <- NA; lofmis_count <- NA; lofmis_low95 <- NA; lofmis_hi95 <- NA
lof_expected <- NA; lof_actual <- NA; lof_enrichment <- NA; lof_pval <- NA; lof_count <- NA; lof_low95 <- NA; lof_hi95 <- NA
mis_expected <- NA; mis_actual <- NA; mis_enrichment <- NA; mis_pval <- NA; mis_count <- NA; mis_low95 <- NA; mis_hi95 <- NA
syn_expected <- NA; syn_actual <- NA; syn_enrichment <- NA; syn_pval <- NA; syn_count <- NA; syn_low95 <- NA; syn_hi95 <- NA


## -- unaffected counts
unaff_all <- unaff$V1
unaff_lofmis <- subset(unaff$V1,unaff$V2=='ptv' | unaff$V2=='missense')
unaff_lof <- subset(unaff$V1,unaff$V2=='ptv')
unaff_mis <- subset(unaff$V1,unaff$V2=='missense')
unaff_syn <- subset(unaff$V1,unaff$V2=='synonymous')

unaff_all_dnm <- rep(length(unaff$V1),length(setlist_name))
unaff_lofmis_dnm <- rep(sum(unaff$V2=='ptv' | unaff$V2=='missense'),length(setlist_name))
unaff_lof_dnm <- rep(sum(unaff$V2=='ptv'),length(setlist_name))
unaff_mis_dnm <- rep(sum(unaff$V2=='missense'),length(setlist_name))
unaff_syn_dnm <- rep(sum(unaff$V2=='synonymous'),length(setlist_name))

## --- Proportion testing
unaff_all_count <- NA; all_prop <- NA; unaff_all_prop <- NA; all_prop_enrichment <- NA; all_prop_pval <- NA; all_prop_low95 <- NA; all_prop_hi95 <- NA
unaff_lofmis_count <- NA; lofmis_prop <- NA; unaff_lofmis_prop <- NA; lofmis_prop_enrichment <- NA; lofmis_prop_pval <- NA; lofmis_prop_low95 <- NA; lofmis_prop_hi95 <- NA
unaff_lof_count <- NA; lof_prop <- NA; unaff_lof_prop <- NA; lof_prop_enrichment <- NA; lof_prop_pval <- NA; lof_prop_low95 <- NA; lof_prop_hi95 <- NA
unaff_mis_count <- NA; mis_prop <- NA; unaff_mis_prop <- NA; mis_prop_enrichment <- NA; mis_prop_pval <- NA; mis_prop_low95 <- NA; mis_prop_hi95 <- NA
unaff_syn_count <- NA; syn_prop <- NA; unaff_syn_prop <- NA; syn_prop_enrichment <- NA; syn_prop_pval <- NA; syn_prop_low95 <- NA; syn_prop_hi95 <- NA


## --- LooP through gene sets
for (i in 1:length(setlist_name)) {

    ll <- subset(setlist[,1],setlist[,2]==setlist_name[i])

    genes_listed[i] <- length(ll)
    genes_missed[i] <- sum(!(ll %in% genes$gene))
    genes_tested[i] <- sum(ll %in% genes$gene)
    genes_pct[i] <- sum(ll %in% genes$gene) / length(ll)   
    ll2 <- ll[!(ll %in% genes$gene)]
    genenames_missed[i] <- paste(ll2,collapse=';')

    ## get mutation expectation
    gene2 <- subset(genes,genes$gene %in% ll) 
    all_expected[i] <- sum(gene2$p_all) / sum(genes$p_all)
    lofmis_expected[i] <- (sum(gene2$p_lof)+sum(gene2$p_mis)) / (sum(genes$p_lof) + sum(genes$p_mis))
    lof_expected[i] <- sum(gene2$p_lof) / sum(genes$p_lof)
    mis_expected[i] <- sum(gene2$p_mis) / sum(genes$p_mis)
    syn_expected[i] <- sum(gene2$p_syn) / sum(genes$p_syn)

    ## get DNM count
    all_count[i] <- sum(all %in% gene2$gene)
    lofmis_count[i] <- sum(lofmis %in% gene2$gene)
    lof_count[i] <- sum(lof %in% gene2$gene)
    mis_count[i] <- sum(mis %in% gene2$gene)
    syn_count[i] <- sum(syn %in% gene2$gene)

    ## Run binomal test
    if (all_dnm[i] == 0) { all_pval[i] <- NA; all_low95[i] <- NA; all_hi95[i] <- NA; all_enrichment[i] <- NA }
    if (all_dnm[i] > 0) { 
    all_binom <- binom.test(x=all_count[i],n=all_dnm[i],p=all_expected[i],alternative ="two.sided")
    all_pval[i] <- all_binom$p.value ; all_actual[i] <- all_binom$estimate 
    all_low95[i] <- all_binom$conf.int[1] ; all_hi95[i] <- all_binom$conf.int[2] 
    all_enrichment[i] <- all_actual[i] / all_expected[i]
    }

    if (lofmis_dnm[i] == 0) { lofmis_pval[i] <- NA; lofmis_low95[i] <- NA; lofmis_hi95[i] <- NA; lofmis_enrichment[i] <- NA }
    if (lofmis_dnm[i] > 0) { 
    lofmis_binom <- binom.test(x=lofmis_count[i],n=lofmis_dnm[i],p=lofmis_expected[i],alternative ="two.sided")
    lofmis_pval[i] <- lofmis_binom$p.value ; lofmis_actual[i] <- lofmis_binom$estimate 
    lofmis_low95[i] <- lofmis_binom$conf.int[1] ; lofmis_hi95[i] <- lofmis_binom$conf.int[2] 
    lofmis_enrichment[i] <- lofmis_actual[i] / lofmis_expected[i]
    }

    if (lof_dnm[i] == 0) { lof_pval[i] <- NA; lof_low95[i] <- NA; lof_hi95[i] <- NA; lof_enrichment[i] <- NA }
    if (lof_dnm[i] > 0) { 
    lof_binom <- binom.test(x=lof_count[i],n=lof_dnm[i],p=lof_expected[i],alternative ="two.sided")
    lof_pval[i] <- lof_binom$p.value ; lof_actual[i] <- lof_binom$estimate 
    lof_low95[i] <- lof_binom$conf.int[1] ; lof_hi95[i] <- lof_binom$conf.int[2] 
    lof_enrichment[i] <- lof_actual[i] / lof_expected[i]
    }

    if (mis_dnm[i] == 0) { mis_pval[i] <- NA; mis_low95[i] <- NA; mis_hi95[i] <- NA; mis_enrichment[i] <- NA }
    if (mis_dnm[i] > 0) { 
    mis_binom <- binom.test(x=mis_count[i],n=mis_dnm[i],p=mis_expected[i],alternative ="two.sided")
    mis_pval[i] <- mis_binom$p.value ; mis_actual[i] <- mis_binom$estimate 
    mis_low95[i] <- mis_binom$conf.int[1] ; mis_hi95[i] <- mis_binom$conf.int[2] 
    mis_enrichment[i] <- mis_actual[i] / mis_expected[i]
    }

    if (syn_dnm[i] == 0) { syn_pval[i] <- NA; syn_low95[i] <- NA; syn_hi95[i] <- NA; syn_enrichment[i] <- NA }
    if (syn_dnm[i] > 0) { 
    syn_binom <- binom.test(x=syn_count[i],n=syn_dnm[i],p=syn_expected[i],alternative ="two.sided")
    syn_pval[i] <- syn_binom$p.value ; syn_actual[i] <- syn_binom$estimate 
    syn_low95[i] <- syn_binom$conf.int[1] ; syn_hi95[i] <- syn_binom$conf.int[2] 
    syn_enrichment[i] <- syn_actual[i] / syn_expected[i]
    }


    ## get unaffected DNM count
    unaff_all_count[i] <- sum(unaff_all %in% gene2$gene)
    unaff_lofmis_count[i] <- sum(unaff_lofmis %in% gene2$gene)
    unaff_lof_count[i] <- sum(unaff_lof %in% gene2$gene)
    unaff_mis_count[i] <- sum(unaff_mis %in% gene2$gene)
    unaff_syn_count[i] <- sum(unaff_syn %in% gene2$gene)

    ## Run proportion test
    if (all_dnm[i] == 0 | unaff_all_dnm[i] == 0) { all_prop[i] <- NA; all_prop_pval[i] <- NA; all_prop_low95[i] <- NA; all_prop_hi95[i] <- NA; all_prop_enrichment[i] <- NA }
    if (all_dnm[i] > 0 & unaff_all_dnm[i] > 0) { 
    all_proptest <- suppressWarnings(prop.test(c(all_count[i],unaff_all_count[i]),c(all_dnm[i],unaff_all_dnm[i]),alternative ="two.sided",correct=F))
    all_prop[i] <- all_proptest$estimate[1]; unaff_all_prop[i] <- all_proptest$estimate[2]
    all_prop_pval[i] <- all_proptest$p.value 
    all_prop_low95[i] <- all_proptest$conf.int[1] ; all_prop_hi95[i] <- all_proptest$conf.int[2] 
    all_prop_enrichment[i] <- all_prop[i] / unaff_all_prop[i]
    }

    if (lofmis_dnm[i] == 0 | unaff_lofmis_dnm[i] == 0) { lofmis_prop[i] <- NA; lofmis_prop_pval[i] <- NA; lofmis_prop_low95[i] <- NA; lofmis_prop_hi95[i] <- NA; lofmis_prop_enrichment[i] <- NA }
    if (lofmis_dnm[i] > 0 & unaff_lofmis_dnm[i] > 0) { 
    lofmis_proptest <- suppressWarnings(prop.test(c(lofmis_count[i],unaff_lofmis_count[i]),c(lofmis_dnm[i],unaff_lofmis_dnm[i]),alternative ="two.sided",correct=F))
    lofmis_prop[i] <- lofmis_proptest$estimate[1]; unaff_lofmis_prop[i] <- lofmis_proptest$estimate[2]
    lofmis_prop_pval[i] <- lofmis_proptest$p.value 
    lofmis_prop_low95[i] <- lofmis_proptest$conf.int[1] ; lofmis_prop_hi95[i] <- lofmis_proptest$conf.int[2] 
    lofmis_prop_enrichment[i] <- lofmis_prop[i] / unaff_lofmis_prop[i]
    }

    if (lof_dnm[i] == 0 | unaff_lof_dnm[i] == 0) { lof_prop[i] <- NA; lof_prop_pval[i] <- NA; lof_prop_low95[i] <- NA; lof_prop_hi95[i] <- NA; lof_prop_enrichment[i] <- NA }
    if (lof_dnm[i] > 0 & unaff_lof_dnm[i] > 0) { 
    lof_proptest <- suppressWarnings(prop.test(c(lof_count[i],unaff_lof_count[i]),c(lof_dnm[i],unaff_lof_dnm[i]),alternative ="two.sided",correct=F))
    lof_prop[i] <- lof_proptest$estimate[1]; unaff_lof_prop[i] <- lof_proptest$estimate[2]
    lof_prop_pval[i] <- lof_proptest$p.value 
    lof_prop_low95[i] <- lof_proptest$conf.int[1] ; lof_prop_hi95[i] <- lof_proptest$conf.int[2] 
    lof_prop_enrichment[i] <- lof_prop[i] / unaff_lof_prop[i]
    }

    if (mis_dnm[i] == 0 | unaff_mis_dnm[i] == 0) { mis_prop[i] <- NA; mis_prop_pval[i] <- NA; mis_prop_low95[i] <- NA; mis_prop_hi95[i] <- NA; mis_prop_enrichment[i] <- NA }
    if (mis_dnm[i] > 0 & unaff_mis_dnm[i] > 0) { 
    mis_proptest <- suppressWarnings(prop.test(c(mis_count[i],unaff_mis_count[i]),c(mis_dnm[i],unaff_mis_dnm[i]),alternative ="two.sided",correct=F))
    mis_prop[i] <- mis_proptest$estimate[1]; unaff_mis_prop[i] <- mis_proptest$estimate[2]
    mis_prop_pval[i] <- mis_proptest$p.value 
    mis_prop_low95[i] <- mis_proptest$conf.int[1] ; mis_prop_hi95[i] <- mis_proptest$conf.int[2] 
    mis_prop_enrichment[i] <- mis_prop[i] / unaff_mis_prop[i]
    }

    if (syn_dnm[i] == 0 | unaff_syn_dnm[i] == 0) { syn_prop[i] <- NA; syn_prop_pval[i] <- NA; syn_prop_low95[i] <- NA; syn_prop_hi95[i] <- NA; syn_prop_enrichment[i] <- NA }
    if (syn_dnm[i] > 0 & unaff_syn_dnm[i] > 0) {   
    syn_proptest <- suppressWarnings(prop.test(c(syn_count[i],unaff_syn_count[i]),c(syn_dnm[i],unaff_syn_dnm[i]),alternative ="two.sided",correct=F))
    syn_prop[i] <- syn_proptest$estimate[1]; unaff_syn_prop[i] <- syn_proptest$estimate[2]
    syn_prop_pval[i] <- syn_proptest$p.value 
    syn_prop_low95[i] <- syn_proptest$conf.int[1] ; syn_prop_hi95[i] <- syn_proptest$conf.int[2] 
    syn_prop_enrichment[i] <- syn_prop[i] / unaff_syn_prop[i]
    }

    ## get DNM list
    lof_list <- lof[lof %in% gene2$gene]
    mis_list <- mis[mis %in% gene2$gene]
    syn_list <- syn[syn %in% gene2$gene]

    ## combine and write to file
    full_list <- rbind(c('PTV overlap = ',paste(lof_list,collapse=';')),
    	      c('Missense overlap = ',paste(mis_list,collapse=';')),
    	      c('Synonymous overlap = ',paste(syn_list,collapse=';')))

    ## write to list 
    write.table(full_list,paste0(overlap_path,setlist_name[i],'.overlap'),col=F,row=F,quo=F,sep='\t')

    print(i) 

} ## END of i LooP


## combine results
set_data <- cbind.data.frame(setlist_name,genes_listed,genes_missed,genes_tested,genes_pct,
	 all_dnm,unaff_all_dnm,lofmis_dnm,unaff_lofmis_dnm,lof_dnm,unaff_lof_dnm,mis_dnm,unaff_mis_dnm,syn_dnm,unaff_syn_dnm,
	 all_count,unaff_all_count,lofmis_count,unaff_lofmis_count,lof_count,unaff_lof_count,mis_count,unaff_mis_count,syn_count,unaff_syn_count,
    all_expected,all_actual,all_enrichment,all_pval,all_low95,all_hi95,
    all_prop,unaff_all_prop,all_prop_enrichment,all_prop_pval,all_prop_low95,all_prop_hi95,
    lofmis_expected,lofmis_actual,lofmis_enrichment,lofmis_pval,lofmis_low95,lofmis_hi95,
    lofmis_prop,unaff_lofmis_prop,lofmis_prop_enrichment,lofmis_prop_pval,lofmis_prop_low95,lofmis_prop_hi95,
    lof_expected,lof_actual,lof_enrichment,lof_pval,lof_low95,lof_hi95,
    lof_prop,unaff_lof_prop,lof_prop_enrichment,lof_prop_pval,lof_prop_low95,lof_prop_hi95,
    mis_expected,mis_actual,mis_enrichment,mis_pval,mis_low95,mis_hi95,
    mis_prop,unaff_mis_prop,mis_prop_enrichment,mis_prop_pval,mis_prop_low95,mis_prop_hi95,
    syn_expected,syn_actual,syn_enrichment,syn_pval,syn_low95,syn_hi95,
    syn_prop,unaff_syn_prop,syn_prop_enrichment,syn_prop_pval,syn_prop_low95,syn_prop_hi95,
    genenames_missed)

## write to file
write.table(set_data,out,col=T,row=F,quo=F,sep='\t')


} ## END of all_types==T




## PART 4: Running enrichment with custom annotation (all_types==F)

if (all_types==F) {

## get DNM counts
all <- dnm$V1
all_dnm <- rep(length(dnm$V1),length(setlist_name))

## --- Gene matching
genes_listed <- NA
genes_missed <- NA
genes_tested <- NA
genes_pct <- NA
genenames_missed <- NA

## --- Binomal testing
all_expected <- NA; all_actual <- NA; all_enrichment <- NA; all_pval <- NA; all_count <- NA; all_low95 <- NA; all_hi95 <- NA


## -- unaffected counts
unaff_all <- unaff$V1
unaff_all_dnm <- rep(length(unaff$V1),length(setlist_name))

## --- Proportion testing
unaff_all_count <- NA; all_prop <- NA; unaff_all_prop <- NA; all_prop_enrichment <- NA; all_prop_pval <- NA; all_prop_low95 <- NA; all_prop_hi95 <- NA


## --- LooP through gene sets
for (i in 1:length(setlist_name)) {

    ll <- subset(setlist[,1],setlist[,2]==setlist_name[i])

    genes_listed[i] <- length(ll)
    genes_missed[i] <- sum(!(ll %in% genes$gene))
    genes_tested[i] <- sum(ll %in% genes$gene)
    genes_pct[i] <- sum(ll %in% genes$gene) / length(ll)   
    ll2 <- ll[!(ll %in% genes$gene)]
    genenames_missed[i] <- paste(ll2,collapse=';')

    ## get mutation expectation
    gene2 <- subset(genes,genes$gene %in% ll) 
    all_expected[i] <- sum(gene2$p_all) / sum(genes$p_all)

    ## get DNM count
    all_count[i] <- sum(all %in% gene2$gene)

    ## Run binomal test
    if (all_dnm[i] == 0) { all_pval[i] <- NA; all_low95[i] <- NA; all_hi95[i] <- NA; all_enrichment[i] <- NA }
    if (all_dnm[i] > 0) { 
    all_binom <- binom.test(x=all_count[i],n=all_dnm[i],p=all_expected[i],alternative ="two.sided")
    all_pval[i] <- all_binom$p.value ; all_actual[i] <- all_binom$estimate 
    all_low95[i] <- all_binom$conf.int[1] ; all_hi95[i] <- all_binom$conf.int[2] 
    all_enrichment[i] <- all_actual[i] / all_expected[i]
    }


    ## get unaffected DNM count
    unaff_all_count[i] <- sum(unaff_all %in% gene2$gene)

    ## Run proportion test
    if (all_dnm[i] == 0 | unaff_all_dnm[i] == 0) { all_prop[i] <- NA; all_prop_pval[i] <- NA; all_prop_low95[i] <- NA; all_prop_hi95[i] <- NA; all_prop_enrichment[i] <- NA }
    if (all_dnm[i] > 0 & unaff_all_dnm[i] > 0) { 
    all_proptest <- suppressWarnings(prop.test(c(all_count[i],unaff_all_count[i]),c(all_dnm[i],unaff_all_dnm[i]),alternative ="two.sided",correct=F))
    all_prop[i] <- all_proptest$estimate[1]; unaff_all_prop[i] <- all_proptest$estimate[2]
    all_prop_pval[i] <- all_proptest$p.value 
    all_prop_low95[i] <- all_proptest$conf.int[1] ; all_prop_hi95[i] <- all_proptest$conf.int[2] 
    all_prop_enrichment[i] <- all_prop[i] / unaff_all_prop[i]
    }

    ## get DNM list
    overlap_list <- all[all %in% gene2$gene]

    # write to list 
    write.table(overlap_list,paste0(overlap_path,setlist_name[i],'.overlap'),col=F,row=F,quo=F,sep='\t')

    print(i) 

} ## END of i LooP


## combine results
set_data <- cbind.data.frame(setlist_name,genes_listed,genes_missed,genes_tested,genes_pct,
	 all_dnm,unaff_all_dnm,all_count,unaff_all_count,
    all_expected,all_actual,all_enrichment,all_pval,all_low95,all_hi95,
    all_prop,unaff_all_prop,all_prop_enrichment,all_prop_pval,all_prop_low95,all_prop_hi95,
    genenames_missed)

## write to file
write.table(set_data,out,col=T,row=F,quo=F,sep='\t')


} ## END of all_types==F



# END of geneset_enrichment.R
