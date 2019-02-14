# Table of Contents
* [Overview](#overview) 
* [single_gene_enrichment.R](#single_gene_enrichment.R)

### Overview

Single gene tests for rate of DNMs given your trio count
  * Details of the test can be found in the supplementary material here: https://personal.broadinstitute.org/howrigan/taiwanese_trios_SCZ_supplement/Taiwanese-trio_Supplement_BioArxiv.pdf 
  * Supplementary Section 7 - Statistical tests used to analyze DNM rates and patterns
  * Supplementary Section 16 - Single gene DNM enrichment


### single_gene_enrichment.R

Runs a one-sided poisson rate test on individual genes relative to mutation expectations 
  * Requires per-gene mutation expectations for PTV, missense, and synonymous coding variation
  * Requires total proband count and female proband count (for X-linked genes) 

