"""
:File: FamilyPed.py
:Author: Jack A. Kosmicki
:Last updated: 2015-04-23

FamilyPed.py reads in the Pedigree file (.ped file).

The .ped file must be tab delimited ('\t') and have the following 5 columns in this order:
    family_ID   indiv_ID    father_ID   mother_ID   sex     affected_Status
See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml for more information.

FamilyPed.py stores the family relationships into three hash tables, one for trios,
    one for cases, and one for controls.  An additional final hash table
    contains all of the case/control or trio individuals and their index in the vcf.
"""


import sys
from sets import Set


def readFamily(Ped_File, vcfIndivs, unaff_Flag):
    """ Read in the pedigree file for the trios (.ped file).

    Parameters
    ----------
    Ped_File: Pedigree file (.ped) of the individuals in the vcf.
        columns: FamilyID, IndividualID, dadID, momID, sex, phenotype
        sex: 1=male, 2=female, 3=unknown
        phenotype: 1=unaffected, 2=affected

    vcfIndivs: Array of individuals in the vcf file.
    unaff_Flag: Boolean indicating whether TDT is run on affected/unaffected individuals
            True: unaffected    False: affected
    """

    family = {}     # family hash table    Key = ID   Value = (Father ID, Mother ID, Sex)
    indivs = {}     # individuals in VCF   Key = ID   Value = index in the vcf

    indivSet = Set(vcfIndivs)  # Convert array to Set to decrease lookup time.

    with open(Ped_File) as file:
        for line in file:
            field = line.strip().split('\t')

            family_ID = field[0]
            indiv_ID = field[1]
            father_ID = field[2]
            mother_ID = field[3]

            if indiv_ID not in indivSet:
                sys.stderr.write('Individual {} is not in vcf.\n'.format(indiv_ID))
                continue

            if field[4] == '1':
                sex = 'male'
            elif field[4] == '2':
                sex = 'female'
            else:
                sex = 'NA'

            # Parents, cases, and controls will not have parental IDs.
            if(father_ID == '0' or mother_ID == '0'):
                continue

            # Check to see if the parents are in the vcf.
            if father_ID not in indivSet or mother_ID not in indivSet:
                sys.stderr.write('Family {} is incomplete.\n'.format(family_ID))
                continue

            # If we only want affected probands.
            if not unaff_Flag:
                if field[5] != '2':
                    continue
            # If we are only looking at unaffected probands.
            else:
                if field[5] != '1':
                    continue

            # Family dictionary is in the form: {child_ID} = [Dad_ID, Mom_ID, Sex]
            family[indiv_ID] = (father_ID, mother_ID, sex)
            indivs[indiv_ID] = vcfIndivs.index(indiv_ID)
            indivs[father_ID] = vcfIndivs.index(father_ID)
            indivs[mother_ID] = vcfIndivs.index(mother_ID)

    print 'Number of families in hash table = {}.'.format(len(family))
    return family, indivs


def readCC(Ped_File, vcfIndivs):
    """ Read in the pedigree file for cases and controls (.ped file).

    Parameters
    ----------
    Ped_File: pedigree file (.ped) of the individuals in the vcf.
        columns = FamilyID, IndividualID, dadID, momID, sex, phenotype
        sex: 1=male, 2=female, 3=unknown
        phenotype: 1=unaffected, 2=affected

    vcfIndivs: Array of individuals in the vcf file.
    """

    case = {}         # case hash table:     Key = ID   Value = Sex
    control = {}      # control hash table:  Key = ID   Value = Sex
    caseControl = {}  # cases and controls hash table:   Key = ID   Value = index in vcf

    indivSet = Set(vcfIndivs)  # convert array to Set to decrease lookup time.

    with open(Ped_File) as file:
        for line in file:
            field = line.strip().split('\t')

            indiv_ID = field[1]
            father_ID = field[2]
            mother_ID = field[3]
            ptype = field[5]      # case/control status: 1=control, 2=case

            if indiv_ID not in indivSet:
                sys.stderr.write('Individual {} is not in vcf.\n'.format(indiv_ID))
                continue

            if field[4] == '1':
                sex = 'male'
            elif field[4] == '2':
                sex = 'female'
            else:
                sex = 'NA'

            if(father_ID != '0' or mother_ID != '0'):
                continue

            elif(ptype == '2'):
                case[indiv_ID] = sex
                caseControl[indiv_ID] = vcfIndivs.index(indiv_ID)

            elif(ptype == '1'):
                control[indiv_ID] = sex
                caseControl[indiv_ID] = vcfIndivs.index(indiv_ID)

    print 'Number of cases in hash table = {}.'.format(len(case))
    print 'Number of controls in hash table = {}.'.format(len(control))
    return case, control, caseControl


def getGender(gender):
    """ """
    Pass