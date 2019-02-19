# Table of Contents
* [Overview](#overview) 
* [DNM_burden.R](#DNM_burden.R)

### Overview

Exome-wide tests for rate of DNMs given your trio count
  * Details of the test can be found in the supplementary material here: https://personal.broadinstitute.org/howrigan/taiwanese_trios_SCZ_supplement/Taiwanese-trio_Supplement_BioArxiv.pdf 
  * Supplementary Section 7 - Statistical tests used to analyze DNM rates and patterns
  * Supplementary Section 10 - DNM burden in combined SCZ cohorts


### DNM_burden.R

Interactive script to perform a two-sided poisson rate test across the exome relative to controls and mutation model expectations 
  * Requires exome-wide mutation expectations for PTV, missense, and synonymous coding variation

```
## ----- DESCRIPTION:

## Interactive Script to generate DNM burden results

## ----- SCRIPT SET UP:

## PART 1: read in relevant files
## PART 2: aggregate DNMs by cohort
## PART 3: test DNM burden in SCZ vs control
## PART 4: test DNM burden in SCZ vs mutation model


## ----- NOTES:

## The mutation model used has already adjusted X chromosome DNM rates to the proportion of male probands among bothe SCZ and control cohorts

```

