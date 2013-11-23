#! /usr/bin/env python
import __builtin__
import itertools
import rna

__author__ = 'vladson'


class Dna:
    def __init__(self, str=''):
        self.genome = str.upper()
        self.conversion_table = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        self.alfabet = 'A', 'C', 'T', 'G'

    def __repr__(self):
        return self.genome

    def __eq__(self, other):
        return self.genome == other.genome

    def complementary(self):
        """
        >>> dna = Dna('AAAACCCGGT')
        >>> dna.complementary()
        ACCGGGTTTT
        """
        return self.__class__(''.join(map(lambda c: self.conversion_table[c], self.genome))[::-1])

    def to_rna(self):
        """
        >>> r = Dna('ATGGCA')
        >>> r.to_rna()
        AUGGCA
        """
        return rna.Rna(''.join(map(lambda c: 'U' if c == 'T' else c, self.genome)))

    def substr_finder(self, length=3):
        results = {}
        for substr in self.n_substr_generator(length):
            if results.has_key(substr):
                results[substr] += 1
            else:
                results[substr] = 1
        return results

    def n_substr_finder(self, length, treshold):
        """
        >>> dna = Dna("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA")
        >>> dna.n_substr_finder(5, 4)
        ['CGACA', 'GAAGA']
        """
        results = self.substr_finder(length)
        return map(lambda (key, _): key, filter(lambda (_, val): val >= treshold, results.iteritems()))

    def n_mismatch_substr_finder(self, k, d):
        """
        k is length of k-mer
        d is number of possible mismatches
        >>> dna = Dna('ACGTTGCATGTCGCATGATGCATGAGAGCT')
        >>> dna.n_mismatch_substr_finder(4, 1)
        ['GATG', 'ATGC', 'ATGT']
        #>>> dna = Dna('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC')
        #>>> dna.n_mismatch_substr_finder(10, 2)
        #['GCACACAGAC', 'GCGCACACAC']
        """
        # Generating raw possible kmers
        #results = {}
        #for substr in self.n_substr_generator(k):
        #    comparator, count = self.get_num_mismatch_comparator(substr, d), 1
        #    for kmer, d in results.iteritems():
        #        if d['comparator'](substr):
        #            results[kmer]['count'] += 1
        #        if comparator(kmer):
        #            count += 1
        #    if not results.has_key(substr):
        #        results[substr] = {'comparator': comparator, 'count': count}
        #
        ## selecting candidats
        #counts = map(lambda (a, b): b['count'], results.iteritems())
        #treshold = min(sorted(counts)[len(counts) / 2:])
        #kmers = map(lambda (key, _): key, filter(lambda (_, val): val['count'] >= treshold, results.iteritems()))
        #kmers += results.keys()
        #kmers = set(kmers)
        kmers = map(lambda t: ''.join(t), map(lambda i: itertools.permutations(i, k), itertools.combinations_with_replacement(self.alfabet, k)))
        candidates = dict(zip((kmer for kmer in kmers), map(lambda kmer: dict(
            comparator=self.get_num_mismatch_comparator(kmer, d), count=0), kmers)))
        print len(candidates)
        #counting
        for substr in self.n_substr_generator(k):
             for kmer, d in candidates.iteritems():
                if d['comparator'](substr):
                    candidates[kmer]['count'] += 1

        #select winners
        treshold = max(map(lambda (_, v): v['count'], candidates.iteritems()))
        print treshold
        return map(lambda (key, _): key, filter(lambda (_, val): val['count'] == treshold, candidates.iteritems()))

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

    def apprx_starting_positions(self, pattern, max_mismatches):
        """
        >>> dna = Dna('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT')
        >>> dna.apprx_starting_positions('ATTCTGGA', 3)
        [6, 7, 26, 27]
        """
        positions = []
        pattern_length = len(pattern)
        comparator = self.get_num_mismatch_comparator(pattern, max_mismatches)
        for i in xrange(len(self.genome) - pattern_length + 1):
            if comparator(self.genome[i:pattern_length + i]):
                positions.append(i)
        return positions

    def batches(self, batch):
        for i in xrange(0, len(self.genome)):
            yield self.__class__(self.genome[i: i + batch])

    def clump_finder(self, batch, length, treshold):
        """
        >>> dna = Dna("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA")
        >>> list(dna.clump_finder(50, 5, 4))
        ['CGACA', 'GAAGA']
        """
        clump = set()
        for batch in self.batches(batch):
            for kmer in batch.n_substr_finder(length, treshold):
                clump.add(kmer)
        return clump

    def skew(self):
        """
        >>> dna = Dna('CATGGGCATCGGCCATACGCC')
        >>> print ' '.join(map(lambda x: str(x), dna.skew()))
        0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
        """
        skew = [0]
        for i in self.genome:
            if i == 'C':
                inc = -1
            elif i == 'G':
                inc = 1
            else:
                inc = 0
            skew.append(skew[-1] + inc)
        return skew

    def min_skew(self):
        """
        >>> dna = Dna('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
        >>> print ' '.join(map(lambda x: str(x), dna.min_skew()))
        11 24
        """
        skew = self.skew()
        return [i for i, x in enumerate(skew) if x == min(skew)]

    # Service functions

    def get_num_mismatch_comparator(self, pattern, max_mismathces):
        def comparator(test_string):
            errors = 0
            for i in xrange(len(pattern)):
                if test_string[i] != pattern[i]:
                    errors += 1
                    if errors > max_mismathces:
                        return False
            return True

        func = comparator
        func.pattern = pattern
        return func

    def n_substr_generator(self, length):
        for i in xrange(0, len(self.genome) - length + 1):
            yield self.genome[i:i + length]

    def permutate_substr(self, substr, level=1):
        """
        >>> d = Dna('ACTG')
        >>> list(d.permutate_substr(d.genome))
        ['CCTG', 'TCTG', 'GCTG', 'AATG', 'ATTG', 'AGTG', 'ACAG', 'ACCG', 'ACGG', 'ACTA', 'ACTC', 'ACTT']
        """
        seen = set()
        for i in xrange(len(substr)):
            for c in self.conversion_table.keys():
                if i == 0:
                    new = c + substr[1:]
                else:
                    new = substr[:i] + c + substr[i + 1:]
                if new != substr and new not in seen:
                    seen.add(new)
                    yield new

    # Human output

    def h_max_substr_finder(self, length=3):
        """
        >>> dna = Dna("ACGTTGCATGTCGCATGATGCATGAGAGCT")
        >>> dna.h_max_substr_finder(4)
        'CATG GCAT'
        """
        results = self.substr_finder(length)
        maximum = max(results.values())
        return ' '.join(map(lambda (key, _): key, filter(lambda (_, v): v == maximum, results.iteritems())))

    def h_starting_positions(self, fragment):
        """
        >>> dna = Dna("GATATATGCATATACTT")
        >>> dna.h_starting_positions('ATAT')
        '1 3 9'
        """
        return ' '.join(map(lambda c: str(c), self.starting_positions(fragment)))

    # File initiated calls
    @classmethod
    def h_file_starting_positions(cls, path):
        data = open(path, 'r')
        fragment = data.readline().strip()
        dna = cls(data.readline().strip())
        data.close()
        return dna.h_starting_positions(fragment)

    @classmethod
    def h_file_clump_finder(cls, path):
        data = open(path, 'r')
        dna = Dna(data.readline().strip())
        length, batch, treshold = map(lambda x: int(x), data.readline().strip().split())
        print dna.clump_finder(batch, length, treshold)

    @classmethod
    def h_file_min_skew(cls, path):
        data = open(path, 'r')
        dna = Dna(data.readline().strip())
        data.close()
        print ' '.join(map(lambda x: str(x), dna.min_skew()))

    @classmethod
    def h_file_apprx_start_pos(cls, path):
        data = open(path, 'r')
        fragment = data.readline().strip()
        dna = Dna(data.readline().strip())
        max_ed = int(data.readline().strip())
        data.close()
        print ' '.join(map(lambda x: str(x), dna.apprx_starting_positions(fragment, max_ed)))