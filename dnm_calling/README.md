# Overview
This subdirectory contains the python scripts that were used to parse the VCF for results and generate parimary annotations, along with the shell scripts that call the python scripts with specific files and filters

de_novo_finder_3.py (Author: Kaitlin Samocha) - script to find de novo mutations in a VCF file

TDT_CC.py (Author: Jack Kosmicki) - script to collect variant level transmission results
 * filters.py = helper script to set genotype filters and identify PAR regions
 * FamilyPed.py = helper script to split trios and case-control individuals from pedigree file
 * processTDT.py = helper script to run TDT test on results

run_lof_annotation.py (Author: Konrad Karczewski) - script to annotation the VCF with VEP (for use on the LSF job system)
 * loftee_utils.py = helper script to parse VEP annotation
