"""
:File: TDT_CC.py
:Author: Jack A. Kosmicki & Kamil Slowikowski
:Last updated: 2016-04-15

for each variant count the number of transmitted, untransmitted variants in families
find variants in cases and controls

Steps:
1)  Read in the family relations file; save the family relations in a global hash table
2)  Read in the individuals in the VCF and save them in a global hash table
3)  Assign the optional parameters as global variables
4)  Read in the VCF file one line at a time
5)  For each line in the VCF, apply filters to determine if it should be examined
6)  Upon passing the filters, examine each individual in the VCF
7)  For doTDT, only look at the probands (ignore cases/controls/parents)
8)  Apply filters to the proband
9)  If the proband passes the filters, apply filters to
    the parents and determine transmission
10) After determining transmission, update the number of
    transmitted and untransmitted variants
11) Print out the number of transmissions and untransmissions

Usage:
    TDT_CC.py doTDT <vcf_File> <ped_File> <outputFile_Name> [options]
    TDT_CC.py doCaseControl <vcf_File> <ped_File> <outputFile_Name> [options]

Options:
    --pl=VALUE           Specify PL threshold [default: 30]
    --dp=DP_VAL          Specify minimum depth [default: 10]
    --ab_Ref=AB_Ref_VAL  Specify allele balance threshold for homoRef [default: 0.1]
    --ab_Het=AB_Het_VAL  Specify allele balance threshold for het [default: 0.3]
    --ab_Alt=AB_Alt_VAL  Specify allele balance threshold for homoAlt [default: 0.9]
    --gq_Par=GQ_PAR_VAL  GQ threshold for parents [default: 30]
    --gq_Kid=GQ_KID_VAL  GQ threshold for child [default: 30]
    --gq_CC=GQ_CC_VAL    GQ threshold for cases and controls [default: 30]
    --pass               Use only variants with FILTER == PASS. Default is False.  -- Note: not implemented yet.
    --vep                Pull VEP annotations. Default is False.
    --unaff              Examine unaffected probands (doTDT only). Default is False.
    -h, --help           Show this message and exit.
    -v, --version        Show the version and exit.
"""


from __future__ import division
import gzip
import sys
import os
import itertools
import cProfile
import pstats
import StringIO
import VCF_VEP
import FamilyPed
import filters
import VEP_Annotation as vepA
import numpy as np
from sets import Set
from docopt import docopt


__version__ = 1.511
__author__ = 'Jack A. Kosmicki <jkosmicki@fas.harvard.edu>'
__date__ = '04/15/2016'


def doTDT(v, family, thresh):
    """ Perform the Transmission Disequilibrium Test (TDT).

    Parameters
    ----------
    v: line of the VCF
    family is a hash table with:
        Key: individual id      Value: [father id, mother id, sex]
    thresh: hash table of thresholds
    """

    # If filters on the line failed move on.
    if not v:
        return None

    TU = [0, 0]         # Array of [transmissions, untransmissions]
    TU_m = [0, 0]       # same as TU but for males
    TU_f = [0, 0]       # same as TU but for females
    mErr = 0            # Count mendelian error: Parents: ref, child: het
    mErr_o = 0          # Count other mendelian errors
    Nproband_alt = 0    # Number of homozygous alt probands that passed all thresholds
    AN = 0              # Number of families that passed all thresholds
    indivs_T = [v['CHROM'], v['POS'], v['ID'], v['REF'], v['ALT']]       # Array of individuals who were transmitted the variant
    indivs_U = [v['CHROM'], v['POS'], v['ID'], v['REF'], v['ALT']]       # Array of individuals who did not receive the variant

    for indiv_id in family.keys():         #loop through all the probands
        # indiv_data is their GT:AD:DP:GQ:PL stats
        indiv_data = v[indiv_id]

        if indiv_data == None:
            continue

        # Apply quality control filters on proband.
        if not filters.passFilters(indiv_data, thresh, thresh['GQ_Kid_Thresh']):
            continue

        # Apply PL filter to child.
        if not filters.PhredScaleFilter(indiv_data, thresh['PL_Thresh']):
            continue

        father = v[family[indiv_id][0]]
        mother = v[family[indiv_id][1]]

        # Check if the parents have the alternate allele
        # so they can pass it on AND apply quality control filters.
        if filters.TDT_Parent_Filters(father, mother, thresh):
            AN += 1         # all individuals in the nuclear family passed the filters

            # TDT operates differently in the hemizygous chromosomes
            # PAR regions defined from
            # http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/
            # in this case we are in the Par region so transmission is normal
            if filters.check_Hemizgyous(v['CHROM'], family[indiv_id][2],
                                        filters.inPar(v['POS'])):
                TU, TU_m, TU_f, mErr, mErr_o, transFlag = numberTransmissions(indiv_data['GT'], father['GT'], mother['GT'], TU, TU_m, TU_f, family[indiv_id][2], False, mErr, mErr_o)

            else:
                TU, TU_m, TU_f, mErr, mErr_o, transFlag = numberTransmissions(indiv_data['GT'], father['GT'], mother['GT'], TU, TU_m, TU_f, family[indiv_id][2], True, mErr, mErr_o)

            if indiv_data['GT'] == 'homoAlt':
                Nproband_alt += 1

            if transFlag == True:               # if the variant was transmitted
                indivs_T.extend((indiv_id, family[indiv_id][2], indiv_data,
                                 family[indiv_id][0], father, family[indiv_id][1],
                                 mother))
            elif transFlag == False:
                indivs_U.extend((indiv_id, family[indiv_id][2], indiv_data,
                                 family[indiv_id][0], father, family[indiv_id][1],
                                 mother))

    # Ignore the cases in which we have 0 transmissions and 0 untransmissions.
    if TU[0] + TU[1] == 0:
        return None

    # Calculate percentage of mendelian errors.
    mendErrorPercent = (mErr + mErr_o) / (TU[0] + TU[1] + mErr + mErr_o)

    if vepFieldNames:
        gene, anno, pph2, sift, lof = vepA.findVariantAnnotation(v, args, vepFieldNames)
    else:
        gene, anno, pph2, sift, lof = ('', '', '', '', '')

    return [v['CHROM'], v['POS'], v['ID'], v['REF'], v['ALT'], v['FILTER'], v['VQSLOD'],
            gene, anno, pph2, sift, lof, v['AF'], v['AC'], AN,
            Nproband_alt, TU[0], TU[1], TU_m[0], TU_m[1], TU_f[0], TU_f[1],
            mErr, mErr_o, mendErrorPercent], indivs_T, indivs_U


def numberTransmissions(kid, dad, mom, TU, TU_m, TU_f, sex, xFlag, mErr, mErr_o):
    """
    Determine the number of transmissions and nontransmissions.
    The X chromosome is different for males so xFlag indicates
         if we should examine those unique cases (which otherwise
         are Mendelian errors).

    Series of cases:
    Kid Dad Mom Transmissions   Untransmissions
    Ref Het Het 0               2
    Ref Ref Het 0               1
    Het Het Het 1               1
    Het Ref Het 1               0
    Het Alt Het 0               1
    Alt Het Het 2               0
    Alt Het Alt 1               0

    - - X CHROM specific cases - -
    Ref Het Het 0               1
    Ref Ref Het 0               1
    Ref Alt Het 0               1
    Alt Het Het 1               0
    Alt Ref Het 1               0
    Alt Alt Het 1               0
    """

    # Child Father Mother: [transmissions, untransmissions]
    transmissions = {
        'homoRef het het': [0, 2],
        'homoRef homoRef het': [0, 1],
        'homoRef het homoRef': [0, 1],
        'het het het': [1, 1],
        'het homoRef het': [1, 0],
        'het het homoRef': [1, 0],
        'het homoAlt het': [0, 1],
        'het het homoAlt': [0, 1],
        'homoAlt het het': [2, 0],
        'homoAlt het homoAlt': [1, 0],
        'homoAlt homoAlt het': [1, 0]
    }

    xTransmissions = {
        'homoRef het het': [0, 1],
        'homoRef homoRef het': [0, 1],
        'homoRef homoAlt het': [0, 1],
        'homoAlt het het': [1, 0],
        'homoAlt homoRef het': [1, 0],
        'homoAlt homoAlt het': [1, 0]
    }

    # mendelian errors we don't care about
    mendelErrors = Set([
        'homoRef homoRef homoAlt',
        'homoRef homoAlt homoRef',
        'homoRef het homoAlt',
        'homoRef homoAlt het',
        'homoRef homoAlt homoAlt',
        'het homoAlt homoAlt',
        'homoAlt homoRef homoRef',
        'homoAlt homoRef het',
        'homoAlt het homoRef',
        'homoAlt homoAlt homoRef',
        'homoAlt homoRef homoAlt'
    ])

    key = '{} {} {}'.format(kid, dad, mom)
    transFlag = False  # flag to indicate whether the variant was transmitted or not

    counts = transmissions.get(key)

    # If we are in a unique x chromosome case.
    if xFlag:
        counts = xTransmissions.get(key)

    # If not in counts, it's a Mendelian Error.
    if key == 'het homoRef homoRef':
        mErr += 1

    if key in mendelErrors and not xFlag:
        mErr_o +=1

    if not counts:
        return TU, TU_m, TU_f, mErr, mErr_o, None

    if counts[0] > 0:    # if the variant was transmitted, indicate so
        transFlag = True

    # Update the total number of transmissions and untransmissions using vector addition.
    TU = [x+y for x,y in zip(TU, counts)]

    if sex == 'male':
        TU_m = [x+y for x,y in zip(TU_m, counts)]
    elif sex == 'female':
        TU_f = [x+y for x,y in zip(TU_f, counts)]

    return TU, TU_m, TU_f, mErr, mErr_o, transFlag


def doCaseControl(v, cases, controls, thresh):
    """ Analyzes case/control data, counting the
        number of reference and alternate alleles.

    Parameters
    ----------
    v: line of the VCF
    case and control are hash tables with:
        Key: individual id      Value: gender
    thresh is a hash table with:
        Key: name of threshold  Value: threshold
    """

    # If filters on the line failed move on.
    if not v:
        return None

    # Count the number of ref and alt alleles in cases and controls.
    caseRefs_m = 0
    caseAlts_m = 0
    controlRefs_m = 0
    controlAlts_m = 0
    caseRefs_f = 0
    caseAlts_f = 0
    controlRefs_f = 0
    controlAlts_f = 0


    # Loop through all the individuals in the case hash table.
    for indiv_id in cases:

        # indiv_data is their GT:AD:DP:GQ:PL stats
        indiv_data = v[indiv_id]

        if indiv_data == None:
            continue

        # Apply filters and update counts   Note: cases[indiv_id] is gender
        if ProcessCC(indiv_data, thresh):
            parFlag = filters.check_Hemizgyous(v['CHROM'], cases[indiv_id],
                                               filters.inPar(v['POS']) )

            if cases[indiv_id] == 'male':
                caseRefs_m, caseAlts_m = Counts(indiv_data, caseRefs_m, caseAlts_m,
                                                parFlag)
            elif cases[indiv_id] == 'female':
                caseRefs_f, caseAlts_f = Counts(indiv_data, caseRefs_f, caseAlts_f,
                                                parFlag)

    # Loop through all indivs in the control hash table.
    for indiv_id in controls:

        # indiv_data is their GT:AD:DP:GQ:PL stats
        indiv_data = v.get(indiv_id)

        if indiv_data == None:
            continue

        # Apply filters and update counts   Note: controls[indiv_id] is gender
        if ProcessCC(indiv_data, thresh):
            parFlag = filters.check_Hemizgyous(v['CHROM'], controls[indiv_id],
                                               filters.inPar(v['POS']) )

            if controls[indiv_id] == 'male':
                controlRefs_m, controlAlts_m = Counts(indiv_data, controlRefs_m,
                                                      controlAlts_m, parFlag)
            elif controls[indiv_id] == 'female':
                controlRefs_f, controlAlts_f = Counts(indiv_data, controlRefs_f,
                                                      controlAlts_f, parFlag)

    caseAlts = caseAlts_f + caseAlts_m
    caseRefs = caseRefs_f + caseRefs_m
    controlAlts = controlAlts_f + controlAlts_m
    controlRefs = controlRefs_f + controlRefs_m

    if (caseAlts + controlAlts == 0) | (caseRefs + controlRefs == 0):
        return None

    AC = caseAlts + controlAlts
    AN = AC + caseRefs + controlRefs
    AF = AC / AN

    if vepFieldNames:
        gene, anno, pph2, sift, lof = vepA.findVariantAnnotation(v, args, vepFieldNames)
    else:
        gene, anno, pph2, sift, lof = ('', '', '', '', '')

    return [v['CHROM'], v['POS'], v['ID'], v['REF'], v['ALT'], v['FILTER'], v['VQSLOD'],
            gene, anno, pph2, sift, lof, AF, AC, AN,
            caseRefs, caseAlts, controlRefs, controlAlts,
            caseRefs_m, caseAlts_m, controlRefs_m, controlAlts_m,
            caseRefs_f, caseAlts_f, controlRefs_f, controlAlts_f]


def ProcessCC(indivAttr, thresh):
    """ Apply filters to the individual: True = pass, False = fail

    Parameters
    ----------
    indivAttr: [GT:AD:DP:GQ:PL]
    thresh: hash table of thresholds
    """

    # Apply quality filters.
    if not filters.passFilters(indivAttr, thresh, thresh['GQ_CC_Thresh']):
        return False

    # Apply PL filter to individual.
    if not filters.PhredScaleFilter(indivAttr, thresh['PL_Thresh']):
        return False

    return True


def Counts(indivAttr, Refcount, Altcount, parFlag):
    """ Update the number of reference and alternate variant counts.

    Parameters
    ----------
    indivAttr: [GT:AD:DP:GQ:PL]
    parFlag indicates if we are in a hemizgyous case or not
        True: not in hemizygous case      False: in hemizygous
    """

    # homoRef has 2 reference alleles
    if indivAttr['GT'] == 'homoRef':
        if parFlag:
            Refcount += 2
            return Refcount, Altcount
        else:
            Refcount += 1
            return Refcount, Altcount

    # het has 1 reference allele and 1 alternate allele
    elif indivAttr['GT'] == 'het':
        if parFlag:
            Refcount += 1
            Altcount += 1
            return Refcount, Altcount
        else:
            return Refcount, Altcount

    # homoAlt has 2 alternate alleles
    elif indivAttr['GT'] == 'homoAlt':
        if parFlag:
            Altcount += 2
            return Refcount, Altcount
        else:
            Altcount += 1
            return Refcount, Altcount

    # If the individual has another genotype then return the input.
    return Refcount, Altcount


if __name__ == "__main__":
    args = docopt(__doc__, version='0.6.1')
    print(args)

    # GLOBAL VARIABLES
    # Case and control are hash tables with:
    #   Key: individual id      Value: gender
    # Family is a hash table with:
    #   Key: individual id      Value: [father id, mother id, sex]

    # Hash table of thresholds
    thresholds = {
        'PL_Thresh': float(args['--pl']),
        'DP_Thresh': float(args['--dp']),
        'AB_Ref_Thresh': float(args['--ab_Ref']),
        'AB_Het_Thresh': float(args['--ab_Het']),
        'AB_Alt_Thresh': float(args['--ab_Alt']),
        'GQ_Kid_Thresh': float(args['--gq_Kid']),
        'GQ_Parent_Thresh': float(args['--gq_Par']),
        'GQ_CC_Thresh': float(args['--gq_CC']) }

    PASS = args['--pass']

    # Set VEP field names line to None unless the VCF is annotated.
    vepFieldNames = None

    fname, fextension = os.path.splitext(args['<outputFile_Name>'])

    writer = open(args['<outputFile_Name>'], 'wb')
    writer.write('\t'.join(['CHROM', 'POSITION', 'ID', 'REF', 'ALT', 'FILTER', 'VQSLOD',
                            'GENE_NAME', 'FUNCTIONAL_CLASS', 'SIFT', 'PolyPhen2',
                            'loftee', 'AF', 'AC', 'AN'])+ '\t')

    fn_open = gzip.open if args['<vcf_File>'].endswith('.gz') else open

    if args['doTDT']:
        pr = cProfile.Profile()
        pr.enable()

        #writer2 = open(fname + '_indivs_T' + fextension, 'wb')
        #writer3 = open(fname + '_indivs_U' + fextension, 'wb')

        # Write out the header.
        writer.write('\t'.join(['Nproband_HomoAlt',
                            'transmitted','untransmitted',
                            'transmitted_m','untransmitted_m',
                            'transmitted_f','untransmitted_f',
                            'De-Novo_Mendel_Errors','Other_Mendel_Errors',
                            'Percent_of_Mendel_Errors']) + '\n')
        #writer2.write('\t'.join(['CHROM','POSITION','ID','REF','ALT']) + '\n')
        #writer3.write('\t'.join(['CHROM','POSITION','ID','REF','ALT']) + '\n')

        with fn_open(args['<vcf_File>']) as fh:
            for line in fh:
                line = line.rstrip('\r\n').rstrip('\n').rstrip('\t')
                if line.startswith('##INFO=<ID=CSQ'):
                    vepFieldNames = line.split('Format: ')[-1].strip('">').split('|')
                if line.startswith('#CHROM'):
                    indivs = line.split('\t')[9:]  # individual IDs in the VCF
                    pedigree, vcfIndivs = FamilyPed.readFamily(args['<ped_File>'],
                                                               indivs, args['--unaff'])
                if line.startswith('#'):
                    continue
                # the following lines contain the variants
                else:
                    result = doTDT(VCF_VEP.parse(line, vcfIndivs, PASS), pedigree,
                                                 thresholds)
                    if result:
                        writer.write('\t'.join(map(str,result[0])) + '\n')
                        #writer2.write('\t'.join(map(str,result[1])) + '\n')
                        #writer3.write('\t'.join(map(str,result[2])) + '\n')

        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print s.getvalue()

    elif args['doCaseControl']:
        # Write out the header.
        writer.write('\t'.join(['caseRefs','caseAlts','controlRefs','controlAlts',
                        'caseRefs_m','caseAlts_m','controlRefs_m','controlAlts_m',
                        'caseRefs_f','caseAlts_f','controlRefs_f','controlAlts_f'])+'\n')

        with fn_open(args['<vcf_File>']) as fh:
            for line in fh:
                line = line.rstrip('\r\n').rstrip('\n').rstrip('\t')
                if line.startswith('##INFO=<ID=CSQ'):
                    vepFieldNames = line.split('Format: ')[-1].strip('">').split('|')
                if line.startswith('#CHROM'):
                    indivs = line.split('\t')[9:]  # individual IDs in the VCF
                    case, control, caseControl = FamilyPed.readCC(args['<ped_File>'],
                                                                  indivs)
                if line.startswith('#'):
                    continue
                else:
                    result = doCaseControl(VCF_VEP.parse(line, caseControl, PASS), case,
                                           control, thresholds)
                    if result:
                        writer.write('\t'.join(map(str,result)) + '\n')
