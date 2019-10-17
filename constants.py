from collections import OrderedDict

CODON_TYPES = {
    '2f': 'FYHQNKDECs',
    '4fcG': 'VSPTA',
    '4f': 'VSPTA'
}
BASES = 'TCAG'
CODONS = [a+b+c for b in BASES for a in BASES for c in BASES]
STOP_CODONS = ['TGA', 'TAG', 'TAA']
AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'
transl_code_i = 'F2L6I3M1V4S4P4T4A4Y2*2H2Q2N2K2D2E2C2*1W1R4s2R2G4'
transl_i = ''.join([
    str(a*int(b)) for a, b in zip(transl_code_i[::2], transl_code_i[1::2])])
GENETIC_CODE_i = OrderedDict(zip(CODONS, transl_i))
