#! /usr/bin/env python
__author__ = 'vladson'

class Dna:

    def __init__(self, str=''):
        self.genome = str.upper()
        self.conversion_table = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

    def complementary(self):
        """
        >>> dna = Dna('AAAACCCGGT')
        >>> dna.complementary()
        'ACCGGGTTTT'
        """
        return ''.join(map(lambda c: self.conversion_table[c], self.genome))[::-1]


    def substr_finder(self, length = 3):
        """
        >>> dna = Dna("ACGTTGCATGTCGCATGATGCATGAGAGCT")
        >>> dna.substr_finder(4)
        'CATG GCAT'
        """
        results = {}
        for i in range(0, len(self.genome) - length):
            substr = self.genome[i:i + length]
            if results.has_key(substr):
                results[substr] += 1
            else:
                results[substr] = 1
        maximum = max(results.values())
        return ' '.join(map(lambda (key,_): key, filter(lambda (_,v): v == maximum, results.iteritems())))

    def starting_positions(self, fragment):
        """
        >>> dna = Dna("GATATATGCATATACTT")
        >>> dna.starting_positions('ATAT')
        [1, 3, 9]
        """
        positions = []
        fragment_length = len(fragment)
        for i in range(0, len(self.genome) - fragment_length):
           if self.genome[i:fragment_length + i] == fragment:
               positions.append(i)
        return positions

    def h_starting_positions(self, fragment):
        


