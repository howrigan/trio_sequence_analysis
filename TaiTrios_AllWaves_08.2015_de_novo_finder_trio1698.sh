#!/bin/bash
python \
/psych/genetics_data/howrigan/projects/taiwan_trio/scripts/de_novo_finder_3.py \
/seq/dax/TaiTrios_AllWaves_08.2015/Exome/v1/TaiTrios_AllWaves_08.2015.vcf.gz \
/psych/genetics_data/howrigan/projects/taiwan_trio/data/Taiwanese_Trio_AllWaves_08.2015_confirmedTrios_wholeBloodProbands.fam \
/humgen/atgu1/fs03/wip/kaitlin/all_ESP_counts_5.28.13.txt > \
/psych/genetics_data/howrigan/projects/taiwan_trio/denovo_data/TaiTrios_AllWaves_08.2015.dnm
