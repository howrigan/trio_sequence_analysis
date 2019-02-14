#!/bin/bash

python TDT_CC.py doTDT \
/seq/dax/TaiTrios_AllWaves_08.2015/Exome/v1/TaiTrios_AllWaves_08.2015.vcf.gz \
/psych/genetics_data/howrigan/projects/taiwan_trio/data/Taiwanese_Trio_AllWaves_08.2015_confirmedTrios_uncontaminatedTrios.fam \
/psych/genetics_data/howrigan/projects/taiwan_trio/transmission_data/variant_transmission/TaiTrios_AllWaves_08.2015_n1142.txt \
--pl 25 --ab_Ref 0.05 --ab_Het 0.2 --ab_Alt 0.95 --gq_Par 25 --gq_Kid 25

