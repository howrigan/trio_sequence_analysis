# trio_sequence_analysis

this is a repository to exchange scripts that analysis trio-based exome sequence data

In general, python scripts parse the VCF for results, and the shell scripts simply run the python scripts with specific files and filters

de_novo_finder.py = Kaitlin Samocha's Script to find de novo mutations
TDT_CC.py = Jack Kosmicki's script to collect variant level transmission results
vep_annotation.py = Konrad K's script to annotation the VCF with VEP (for use on the LSF job system)
