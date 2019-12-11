from collections import OrderedDict

CODON_TYPES = {
    '2f': 'FYHQNKDECs',
    '4fcG': 'VSPTA',
    '4f': 'VSPTAG'
}
BASES = 'TCAG'
CODONS = [a+b+c for b in BASES for a in BASES for c in BASES]
STOP_CODONS = ['TGA', 'TAG', 'TAA']
AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'
transl_code_i = 'F2L6I3M1V4S4P4T4A4Y2*2H2Q2N2K2D2E2C2*1W1R4s2R2G4'
transl_i = ''.join([
    str(a*int(b)) for a, b in zip(transl_code_i[::2], transl_code_i[1::2])])
GENETIC_CODE_i = OrderedDict(zip(CODONS, transl_i))

TC_END_2f = 'FYHNDCs'
AG_END_2f = 'QKE'

TC = list('TC')
TCs = [TC for _ in range(len(TC_END_2f))]
AG = list('AG')
AGs = [AG for _ in range(len(AG_END_2f))]
TCAG = TC + AG
TCAGs = [TCAG for _ in range(len(CODON_TYPES['4f']))]

AA_CODON_3rdpos_BASES = dict(zip(
    TC_END_2f+AG_END_2f+CODON_TYPES['4f'], 
    TCs+AGs+TCAGs
))

# For SFS
FREQUENCY_POOL_CLASS = {
    'Dm' :{
        1: [1], 2: [2], 3: [3], 4:[4, 5], 5:[6, 7],
        6: [8, 9, 10], 7: [11, 12, 13], 8: [14]
    },
    'Ds':{
        1: [1], 2: [2], 3: [3], 4: [4, 5], 5: [6, 7], 6: [8, 9, 10],
        7: [11, 12, 13, 14], 8: [15, 16, 17, 18, 19, 20], 9:[21]
    }
}

MUTATION_POOL_CLASS = {
    'WS': ['AC', 'AG', 'TC', 'TG'],
    'SW': ['CA', 'CT', 'GA', 'GT'],
    'WW': ['AT', 'TA'],'SS': ['CG', 'GC'],
    'neu': ['AT', 'TA', 'CG', 'GC']
}

MUTATIONS = [b1+b2 for b1 in list(BASES) for b2 in list(BASES) if b1!=b2] \
    + list(MUTATION_POOL_CLASS.keys())

