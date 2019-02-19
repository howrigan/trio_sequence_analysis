# Overview
This subdirectory contains python scripts that were parse the VCF for results and generate parimary annotations, along with the shell scripts simply run the python scripts with specific files and filters

de_novo_finder_3.py = Kaitlin Samocha's Script to find de novo mutations in a VCF file

TDT_CC.py = Jack Kosmicki's script to collect variant level transmission results
filters.py = helper script to set genotype filters and identify PAR regions
FamilyPed.py = helper script to split trios and case-control individuals from pedigree file
processTDT.py = helper script to run TDT test on results 

run_lof_annotation.py = Konrad Karczewski's script to annotation the VCF with VEP (for use on the LSF job system)
loftee_utils.py = helper script to parse VEP annotation
