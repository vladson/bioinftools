__author__ = 'vladson'
from protein import Protein

class Rna:

    def __init__(self, rna = ""):
        self.rna = rna
        self.codon_table = {}
        table = open('RNA_codon_table_1.txt', 'r')
        for line in table:
            data = line.strip().split(' ')
            self.codon_table[data[0]] = data[1] if len(data) > 1 else False


    def __repr__(self):
        return self.rna

    def to_protein(self):
        """
        >>> rna = Rna('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
        >>> rna.to_protein()
        MAMAPRTEINSTRING
        """
        return Protein(''.join(self.amynoacids()))

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
        for i in xrange(0, len(self.rna), 3):
            yield self.rna[i:i+3]