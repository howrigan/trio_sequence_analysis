#!/usr/bin/env python
"""
:File:         VCF_VEP.py
:OriginalAuthor:   Kamil Slowikowski <kslowikowski@fas.harvard.edu>
:Edited:       Jack A. Kosmicki
:Last updated: 2015-08-10

Read in a line from a VCF and return a hash table with the columns of the VCF as keys
    and the records for each variant as the values.  The info field is broken up into
    key value pairs.
Hence, vcfLine['CHROM'] will return the chromosome the variant is on.

Note:
vcf lines that do not pass GATK filters will not be processed nor are multi-allelic variants.
"""


import gzip
import sys
from sets import Set


def parse(line, hashTable, PASS):
    """Parse a single VCF line and return a dictionary.

    Parameters
    ----------
    line: a line in the VCF (string)
    hashTables: hash tables of individuals to look up in the VCF
    PASS: Flag indicating whether only PASS sites should be examined
    """

    vcfLine = {}

    # the 1st 8 defined Columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    FIELDS = line.split('\t')[:9]

    vcfLine['FILTER'] = FIELDS[6]

    if vcfLine['FILTER'].startswith('InbreedingCoeff'):
	return None
    
    # if PASS is true, examine only sites whose FILTER is PASS
    if PASS:
        if not vcfLine['FILTER'].startswith('PASS'):
            return None

    # ignore multi-allelics
    if ',' in FIELDS[4]:
        return None

    # load dictionary
    vcfLine['CHROM'] = FIELDS[0]
    vcfLine['POS'] = int(FIELDS[1])
    vcfLine['ID'] = FIELDS[2]
    vcfLine['REF'] = FIELDS[3]
    vcfLine['ALT'] = FIELDS[4]
    vcfLine['QUAL'] = FIELDS[5]

    # VCF file format must be GT:AD:DP:GQ:PL
    if not checkFORMAT(FIELDS[8]):
        return None
    else:
        format = FIELDS[8].split(':')
        format_len = len(format)

    # INFO field consists of "key1=value;key2=value;...".
    infos = FIELDS[7].split(';')

    for i, info in enumerate(infos, 1):
        # It should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        if not value:
            value = None
        vcfLine[key] = value

    if not PASS and float(vcfLine['VQSLOD']) < -2.632:
        return None

    fields = line.split('\t')[9:]       # fields are all of the individuals in the vcf.

    for indivID in hashTable.keys():
        # individuals values based on the format field
        stats = {}

        # indivAttr (individual attributes) are GT:AD:DP:GQ:PL
        indivAttr = fields[hashTable[indivID]].strip('"').split(':')

        # If the sizes don't match, the individual has missing values in the format field.
        # i.e., format = GT:AD:DP and indivAttr = ./.
        if len(indivAttr) != format_len:
            vcfLine[indivID] = None
        else:
            # create a dictionary of GT:AD:DP:GQ:PL
            for i in range(format_len):
                if ',' in indivAttr[i]:
                    try:
                        indivAttr[i] = map(int, indivAttr[i].split(','))
                    except ValueError:
                        vcfLine[indivID] = None
                stats[format[i]] = indivAttr[i]

            # sometimes AD is none
            if stats['AD'] != '.':
                stats['GT'] = find_gtype(stats)
                vcfLine[indivID] = stats
            else:
                vcfLine[indivID] = None

    return vcfLine


def find_gtype(stats):
    """ Determine the genotype of the individual by checking both
        the genotype and AD.

    Parameters
    ----------
    stats = format field [GT,[AD_ref, AD_alt],DP,GQ,[PL_ref, PL_het, PL_alt]]
    """

    # Throw away individuals with 0 reference and alternate reads.
    if stats['AD'][0] == 0 and stats['AD'][1] == 0:
        return None
    elif stats['GT'] in ('0/0', '0|0'):
        return 'homoRef'
    elif stats['GT'] in ('0/1', '0|1', '1/0', '1|0'):
        return 'het'
    elif stats['GT'] in ('1/1', '1|1'):
        return 'homoAlt'
    else:
        return None


def checkFORMAT(formatField):
    """ This program requires the following values in the FORMAT field:
        GT: genotype
        AD: allelic depth
        DP: total read depth
        GQ: Genotype Quality
        PL: Phred-scaled likelihoods

        If any of these fields are missing then the line cannot be read.

        Parameters
        ----------
        formatField: (: delimited string) 9th column of the vcf.
    """

    format = Set(formatField.split(':'))
    if 'GT' not in format:
        sys.stderr.write('WARNING: GT not in format field.')
        return False
    elif 'AD' not in format:
        sys.stderr.write('WARNING: AD not in format field.')
        return False
    elif 'DP' not in format:
        sys.stderr.write('WARNING: DP not in format field.')
        return False
    elif 'GQ' not in format:
        sys.stderr.write('WARNING: GQ not in format field.')
        return False
    elif 'PL' not in format:
        sys.stderr.write('WARNING: PL not in format field.')
        return False
    else:
        return True
