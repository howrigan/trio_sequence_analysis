# Overview
The files are used to run DNM exome-wide burden, gene set, and single gene enrichment analyses

# Files
DNM_studies.tsv
 * list of cohorts, disease designation, and number of trios reported

combined_cohorts_DNM_list.tsv
 * full list of QC-passing (or published) DNM calls from the trios listed in DNM_studies.tsv
 * ALl DNMs aligned to hg19 (GRCh37) reference
 * contains annotations used in subsequent analyses

gencode_pct75_gene17925.tsv
 * Per-gene mutation model probabilities
 * Probabilities split by primary annotation (PTV, missense, synonymous)
 * Using mutation model probabilities from Gencode v19 transcripts (See supplementary section 8 - Mutation rate model testing)
 * Restricted to 17925 'well-covered' genes in the SCZ Taiwanese trio cohort (at least 75% of capture target in gene meets 10x coverage)
 * Probabilities additionally adjusted to coverage and proportion of females in full trio set

gencode_pct75_gene17925_haploid.tsv
 * Similar file to above but using per-chromosome probabilities (as opposed to per-proband probabilities)
 * For autosomal genes = probability * 2
 * For X-linked genes = female probands (prob. * 2) / male probands (prob.)

candidate_genesets.tsv
 * HGNC gene symbols for 85 gene sets used in the main analysis
