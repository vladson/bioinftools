#! /usr/bin/env python
__author__ = 'vladson'

class Dna:

    def __init__(self, str=''):
        self.genome = str.upper()
        self.conversion_table = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

    def complementary(self):
        return ''.join(map(lambda c: self.conversion_table[c], self.genome))[::-1]


    def substr_finder(self, length = 3):
        results = {}
        for i in range(0, len(self.genome) - length):
            substr = self.genome[i:i + length]
            if results.has_key(substr):
                results[substr] += 1
            else:
                results[substr] = 1
        maximum = max(results.values())
        return ' '.join(map(lambda (key,_): key, filter(lambda (_,v): v == maximum, results.iteritems())))

if __name__ == '__main__':
    dna = Dna("ACGTTGCATGTCGCATGATGCATGAGAGCT")
    print(dna.substr_finder(4))
    raw_input()

