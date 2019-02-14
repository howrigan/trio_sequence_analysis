# Table of Contents
* [Overview](#overview) 
* [dnm_list.R](#dnm_list.R)
* [geneset_enrichment.R](#geneset_enrichment.R)
* [process_geneset_enrichment_results.R](#process_geneset_enrichment_results.R)

### Overview

Gene set enrichment tests for relative enrichment of DNM hits in a given gene set relative to the whole exome 
  * Details of the enrichment test can be found in the supplementary material here: https://personal.broadinstitute.org/howrigan/taiwanese_trios_SCZ_supplement/Taiwanese-trio_Supplement_BioArxiv.pdf 
  * Supplementary Section 7 - Statistical tests used to analyze DNM rates and patterns
  * Supplementary Section 14 - DNM burden in candidate gene sets


### dnm_list.R

Extracts simple lists of the gene symbol and primary annotation (PTV / missense / synonymous) from exome-wide coding DNM lists 
  * primary goal - split DNMs by disease or cohort
  * secondary goal - split by custom annotation (not in ExAC or damaging annotation)


### geneset_enrichment.R

Test for enrichment of DNM overlapping genes in a given gene set against:
  * Mutation model expectations using a binomial exact test
  * Another DNM gene list (preferrably unaffected individuals) using a two-sample proportion test

```
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
```

### process_geneset_enrichment_results.R

When looking at default annotations, convert larger results file into multiple annotation-specific files sorted by p-value

```
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
```
