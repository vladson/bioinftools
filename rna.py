import dna
import protein

__author__ = 'vladson'

class Rna(dna.Dna):

    conversion_table = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C'}
    alfabet = 'A', 'C', 'U', 'G'
    codon_table = {}
    table = open('RNA_codon_table_1.txt', 'r')
    for line in table:
        data = line.strip().split(' ')
        codon_table[data[0]] = data[1] if len(data) > 1 else False

    def __init__(self, genome=""):
        self.genome = genome


    def to_protein(self):
        """
        >>> rna = Rna('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
        >>> rna.to_protein()
        MAMAPRTEINSTRING
        """
        return protein.Protein(''.join(self.amynoacids()))

    def to_dna(self):
        """
        >>> r = Rna('AUGGCA')
        >>> r.to_dna()
        ATGGCA
        """
        return dna.Dna(''.join(map(lambda c: 'T' if c == 'U' else c, self.genome)))

    def amynoacids(self):
        """
        >>> r = Rna('AUGGCC')
        >>> print list(r.amynoacids())
        ['M', 'A']
        """
        for codon in self.triplets():
            if self.codon_table[codon]:
                yield self.codon_table[codon]

    def triplets(self):
        """
        >>> r = Rna('AUGGCC')
        >>> print list(r.triplets())
        ['AUG', 'GCC']
        """
        for i in xrange(0, len(self.genome), 3):
            yield self.genome[i:i+3]