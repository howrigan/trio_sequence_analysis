## DNM_burden.R

## Author: Daniel P. Howrigan (daniel.howrigan@gmail.com)
## Last Modified: January 2019


## ----- DESCRIPTION:

## Interactive Script to generate DNM burden results

## ----- SCRIPT SET UP:

## PART 1: read in relevant files
## PART 2: aggregate DNMs by cohort
## PART 3: test DNM burden in SCZ vs control
## PART 4: test DNM burden in SCZ vs mutation model


## ----- NOTES:

## The mutation model used has already adjusted X chromosome DNM rates to the proportion of male probands among bothe SCZ and control cohorts


## PART 1: read in relevant files

## de novo mutation list
dnm <- read.table('../files/combined_cohorts_DNM_list.tsv',h=T,stringsAsFactors=F,sep='\t')

## list of cohorts and trio counts
grps <- read.table('../files/DNM_studies.tsv',h=T,stringsAsFactors=F,sep='\t')

## mutation model
model <- read.table('../files/gencode_pct75_gene17925.tsv',h=T,stringsAsFactors=F,sep='\t')




## PART 2: aggregate DNMs by cohort

## Remove DNMs outside of Agilent Interval
dnm <- subset(dnm,dnm$agilent_interval==T)

## Remove DNMs outside of genes with low 10x coverage
dnm <- subset(dnm,dnm$gene_symbol_used %in% model$gene)

## create per-trio count for various measures
all <- NA
ptv <- NA
mis <- NA
syn <- NA

## aggregate DNMs by cohort
for (i in 1:nrow(grps)) {

    ## subset to cohort
    dnm2 <- subset(dnm,dnm$STUDY %in% grps$STUDY[i])
    dnm2 <- subset(dnm2,dnm2$DISEASE %in% grps$DISEASE[i])

    all[i] <- nrow(dnm2)
    ptv[i] <- sum(dnm2$annotation_used=='ptv')
    mis[i] <- sum(dnm2$annotation_used=='missense')
    syn[i] <- sum(dnm2$annotation_used=='synonymous')

    print(i)
} ## End of i LooP

## add counts to cohort list
grps <- cbind(grps,all,ptv,mis,syn)







## PART 3: test DNM burden in SCZ vs control


## data <- read.table('burden_scz/published_ttrio_dnmCounts_agilentIntervals.txt',h=T,stringsAsFactors=F)
NOW grps

## remove outlier cohorts
grps <- subset(grps,grps$STUDY != 'Xu')

## split by disease
study <- names(table(grps$DISEASE))

## set up variables
trios <- NA
dnm <- NA; dnm_rate <- NA; dnm_lowci <- NA; dnm_highci <- NA; dnm_enrich <- NA; dnm_pval <- NA
ptv <- NA; ptv_rate <- NA; ptv_lowci <- NA; ptv_highci <- NA; ptv_enrich <- NA; ptv_pval <- NA
mis <- NA; mis_rate <- NA; mis_lowci <- NA; mis_highci <- NA; mis_enrich <- NA; mis_pval <- NA
syn <- NA; syn_rate <- NA; syn_lowci <- NA; syn_highci <- NA; syn_enrich <- NA; syn_pval <- NA



for (i in 1:length(study)) {

    grp1 <- subset(grps,grps$DISEASE==study[i])
    grp2 <- subset(grps,grps$DISEASE!=study[i])

    trios[i] <- sum(grp1$TRIOS)

    ## Overall DNM
    mod <- poisson.test(c(sum(grp1$all),sum(grp2$all)),
    	c(sum(grp1$TRIOS),sum(grp2$TRIOS)))
    dnm[i] <- sum(grp1$all)
    dnm_rate[i] <- sum(grp1$all) / sum(grp1$TRIOS)
    grp2_rate <- sum(grp2$all) / sum(grp2$TRIOS)
    dnm_lowci[i] <- as.numeric(mod$conf.int[1]*grp2_rate)
    dnm_highci[i] <- as.numeric(mod$conf.int[2]*grp2_rate)
    dnm_pval[i] <- mod$p.value 
    dnm_enrich[i] <- mod$estimate

    ## PTV
    mod <- poisson.test(c(sum(grp1$ptv),sum(grp2$ptv)),
    	c(sum(grp1$TRIOS),sum(grp2$TRIOS)))
    ptv[i] <- sum(grp1$ptv)
    ptv_rate[i] <- sum(grp1$ptv) / sum(grp1$TRIOS)
    grp2_rate <- sum(grp2$ptv) / sum(grp2$TRIOS)
    ptv_lowci[i] <- as.numeric(mod$conf.int[1]*grp2_rate)
    ptv_highci[i] <- as.numeric(mod$conf.int[2]*grp2_rate)
    ptv_pval[i] <- mod$p.value 
    ptv_enrich[i] <- mod$estimate

    ## MIS
    mod <- poisson.test(c(sum(grp1$mis),sum(grp2$mis)),
    	c(sum(grp1$TRIOS),sum(grp2$TRIOS)))
    mis[i] <- sum(grp1$mis)
    mis_rate[i] <- sum(grp1$mis) / sum(grp1$TRIOS)
    grp2_rate <- sum(grp2$mis) / sum(grp2$TRIOS)
    mis_lowci[i] <- as.numeric(mod$conf.int[1]*grp2_rate)
    mis_highci[i] <- as.numeric(mod$conf.int[2]*grp2_rate)
    mis_pval[i] <- mod$p.value 
    mis_enrich[i] <- mod$estimate

    ## SYN
    mod <- poisson.test(c(sum(grp1$syn),sum(grp2$syn)),
    	c(sum(grp1$TRIOS),sum(grp2$TRIOS)))
    syn[i] <- sum(grp1$syn)
    syn_rate[i] <- sum(grp1$syn) / sum(grp1$TRIOS)
    grp2_rate <- sum(grp2$syn) / sum(grp2$TRIOS)
    syn_lowci[i] <- as.numeric(mod$conf.int[1]*grp2_rate)
    syn_highci[i] <- as.numeric(mod$conf.int[2]*grp2_rate)
    syn_pval[i] <- mod$p.value 
    syn_enrich[i] <- mod$estimate

} ## END of i LooP

(burden_casecon <- cbind.data.frame(study,trios,
dnm, dnm_rate, dnm_lowci, dnm_highci, dnm_enrich, dnm_pval,
ptv, ptv_rate, ptv_lowci, ptv_highci, ptv_enrich, ptv_pval,
mis, mis_rate, mis_lowci, mis_highci, mis_enrich, mis_pval,
syn, syn_rate, syn_lowci, syn_highci, syn_enrich, syn_pval))

## write out to file
## write.table(burden_casecon,'burden_casecon_results.txt',col=T,row=F,quo=F,sep='\t')





## PART 4: test DNM burden in SCZ vs mutation model

## aggregate exome wide mutation model rate (NOTE: X chromosome rates already accounted for in mutation model file)
p_all <- sum(model$p_all)
p_mis <- sum(model$p_mis)
p_ptv <- sum(model$p_lof)
p_syn <- sum(model$p_syn)

exp_rate <- c(p_all,p_mis,p_ptv,p_syn)

## split by disease
study <- names(table(grps$DISEASE))

## set up variables
trios <- NA
dnm <- NA; dnm_rate <- NA; dnm_lowci <- NA; dnm_highci <- NA; dnm_enrich <- NA; dnm_pval <- NA
ptv <- NA; ptv_rate <- NA; ptv_lowci <- NA; ptv_highci <- NA; ptv_enrich <- NA; ptv_pval <- NA
mis <- NA; mis_rate <- NA; mis_lowci <- NA; mis_highci <- NA; mis_enrich <- NA; mis_pval <- NA
syn <- NA; syn_rate <- NA; syn_lowci <- NA; syn_highci <- NA; syn_enrich <- NA; syn_pval <- NA


for (i in 1:length(study)) {

    grp1 <- subset(grps,grps$DISEASE==study[i])

    trios[i] <- sum(grp1$TRIOS)

    ## Overall DNM
    mod <- poisson.test(sum(grp1$all),sum(grp1$TRIOS),p_all)
    dnm[i] <- sum(grp1$all)
    dnm_rate[i] <- sum(grp1$all) / sum(grp1$TRIOS)
    dnm_lowci[i] <- as.numeric(mod$conf.int[1])
    dnm_highci[i] <- as.numeric(mod$conf.int[2])
    dnm_pval[i] <- mod$p.value 
    dnm_enrich[i] <- mod$estimate / p_all

    ## PTV
    mod <- poisson.test(sum(grp1$ptv),sum(grp1$TRIOS),p_ptv)
    ptv[i] <- sum(grp1$ptv)
    ptv_rate[i] <- sum(grp1$ptv) / sum(grp1$TRIOS)
    ptv_lowci[i] <- as.numeric(mod$conf.int[1])
    ptv_highci[i] <- as.numeric(mod$conf.int[2])
    ptv_pval[i] <- mod$p.value 
    ptv_enrich[i] <- mod$estimate / p_ptv

    ## MIS
    mod <- poisson.test(sum(grp1$mis),sum(grp1$TRIOS),p_mis)
    mis[i] <- sum(grp1$mis)
    mis_rate[i] <- sum(grp1$mis) / sum(grp1$TRIOS)
    mis_lowci[i] <- as.numeric(mod$conf.int[1])
    mis_highci[i] <- as.numeric(mod$conf.int[2])
    mis_pval[i] <- mod$p.value 
    mis_enrich[i] <- mod$estimate / p_mis

    ## SYN
    mod <- poisson.test(sum(grp1$syn),sum(grp1$TRIOS),p_syn)
    syn[i] <- sum(grp1$syn)
    syn_rate[i] <- sum(grp1$syn) / sum(grp1$TRIOS)
    syn_lowci[i] <- as.numeric(mod$conf.int[1])
    syn_highci[i] <- as.numeric(mod$conf.int[2])
    syn_pval[i] <- mod$p.value 
    syn_enrich[i] <- mod$estimate / p_syn

} ## END of i LooP

(burden_model <- cbind.data.frame(study,trios,
dnm, dnm_rate, dnm_lowci, dnm_highci, dnm_enrich, dnm_pval,
ptv, ptv_rate, ptv_lowci, ptv_highci, ptv_enrich, ptv_pval,
mis, mis_rate, mis_lowci, mis_highci, mis_enrich, mis_pval,
syn, syn_rate, syn_lowci, syn_highci, syn_enrich, syn_pval))

## write out to file
## write.table(burden_model,'burden_model_results.txt',col=T,row=F,quo=F,sep='\t')
