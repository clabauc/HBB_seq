# HBB_seq
f = open('hbb_mutated_seq_hs.txt', 'r')
hbb_dna = f.readlines()
del hbb_dna[0]
hbb_dna = "".join(i for i in hbb_dna)
hbb_dna = hbb_dna.replace("\n", "")


def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    hbb_protein = []
    if len(seq) % 3 == 0:
        for t_3 in range(0, len(seq), 3):
            codon = seq[t_3:t_3 + 3]
            hbb_protein += (table[codon])
        if hbb_protein[6] == "V":
            return "There is signs of a sickle cell mutation"
        if hbb_protein[6] == "K":
            return "There is signs of a Hemoglobin C mutation"
        if hbb_protein[26] == "K":
            return "There is signs of a Hemoglobin E mutation"
        if hbb_protein[121] == "Q":
            return "There is signs of a Hemoglobin D mutation"
        else:
            return "There is no signs of HBB gene mutations."
    return hbb_protein




#print(hbb_dna)
print(translate(hbb_dna))
