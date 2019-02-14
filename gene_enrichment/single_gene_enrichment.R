## single_gene_enrichment.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: February 2019

## ----- DESCRIPTION:

## Run a poisson test on each gene

## account for X chromosome genes (males + 2*females)

## tests:
# - all DNM
# - ptv
# - ptv+mis
# - mis
# - syn

## poisson.test against mutation expectation

## poisson.test([DNM count],[Inherited chromosome count],[mutation expectation])

# example:
# > poisson.test(2,(2772*2),(8.287213e-05*2))

#   Exact Poisson test

# data:  2 time base: (2772 * 2)
# number of events = 2, time base = 5544, p-value = 0.2344
# alternative hypothesis: true event rate is not equal to 0.0001657443
# 95 percent confidence interval:
#  4.368854e-05 1.303154e-03
# sample estimates:
#   event rate 
# 0.0003607504 



## ----- SCRIPT SET UP:

## PART 1: get chromosome counts
## PART 2: read in relevant files
## PART 3: Run poisson test in each gene
## PART 4: combine results and write to file


## ----- NOTES:

## restricting to genes passing coverage QC 
## Using gencode_pct75_gene17925_haploid.tsv, where the expectation is per-inherited chromosome



## PART 1: get chromosome counts

## inherited autosomal chromosomes (numer of trios * 2)
auto_nchrobs <- 2772*2
## inherited sex chromosomes (number of trios + number of female probands)
x_nchrobs <- 2772+1118



## PART 2: read in relevant files

## get DNM list
dnm <- read.table('../files/combined_cohorts_DNM_list.tsv',h=T,stringsAsFactors=F)
scz <- subset(dnm,dnm$DISEASE=='SCZ')

scz_lof <- subset(scz,scz$annotation_used=='ptv')
scz_mis <- subset(scz,scz$annotation_used=='missense')
scz_lofmis <- subset(scz,scz$annotation_used=='ptv' | scz$annotation_used=='missense')
scz_syn <- subset(scz,scz$annotation_used=='synonymous')

## read in per-gene adjusted mutation rates
gene_list <- read.table('../files/gencode_pct75_gene17925_haploid.tsv',h=T,sep='\t',stringsAsFactors=F)

## get percentiles
gene_list <- gene_list[order(gene_list$bp),]
gene_list$bp_prop <- seq(1,nrow(gene_list),1) / nrow(gene_list)




## PART 3: Run poisson test in each gene

## ======= output variables


mut_type <- c('all','lof','lofmis','mis','syn')

gene <- NA
mut <- NA
count <- NA
exp_rate <- NA
obs_rate <- NA
pval <- NA


X <- 1
## LooP through genes
for (a in 1:nrow(gene_list)) {

    # subset to gene
    gg <- gene_list[a,]

    ## LooP through mut type
    for (b in 1:length(mut_type)) {

    	## collect gene level info
	   gene[X] <- gg$gene
	   	   mut[X] <- mut_type[b]
		   	  if (mut[X]=='all') { dd <- scz ; gm <- as.numeric(gg$p_all) }
			     if (mut[X]=='lof') { dd <- scz_lof ; gm <- as.numeric(gg$p_lof) }
			     	if (mut[X]=='lofmis') { dd <- scz_lofmis ; gm <- as.numeric(gg$p_lof) + as.numeric(gg$p_mis) }
				   if (mut[X]=='mis') { dd <- scz_mis ; gm <- as.numeric(gg$p_mis) }
				      if (mut[X]=='syn') { dd <- scz_syn ; gm <- as.numeric(gg$p_syn) }

				      	 ## subset dd to gene
					    count[X] <- sum(dd$gene_symbol_used %in% gg$gene)

					    	     ## autosomal test
						     	if (gg$chr != 'X') {
							   	   mod <- poisson.test(count[X],auto_nchrobs,gm,alternative='greater')
								       	  exp_rate[X] <- mod$null.value
									  	      	 obs_rate[X] <- mod$estimate
											 	     	pval[X] <- mod$p.value
														} ## END if

														  ## Xchr test
														     if (gg$chr == 'X') {
														     		mod <- poisson.test(count[X],x_nchrobs,gm,alternative='greater')
																       exp_rate[X] <- mod$null.value
																       		      obs_rate[X] <- mod$estimate
																		      		     pval[X] <- mod$p.value
																				     	     } ## END if

																					       print(X)
																					        X <- X+1

																						} ## END of b LooP

} ## END of a LooP




## PART 4: combine results and write to file

res <- cbind.data.frame(gene,
mut,
count,
exp_rate,
obs_rate,
pval)

## combine with gene information
res2 <- merge(res,gene_list,by='gene',all=T)

## order by lowest p-value
res2 <- res2[order(res2$pval),]

write.table(res2,'SCZ_poisson_test_gencode_pct75_gene17925.txt',col=T,row=F,quo=F,sep='\t')



## END of single_gene_enrichment.R
