#! /usr/bin/env python
import itertools

__author__ = 'vladson'


class Dna:
    conversion_table = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    alfabet = 'A', 'C', 'T', 'G'

    def __init__(self, str=''):
        self.genome = str.upper()

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
        import rna
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

    def n_mismatch_substr_finder(self, k, d, cmrtr_type='n_mismatch'):
        """
        k is length of k-mer
        d is number of possible mismatches
        #>>> dna = Dna('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC')
        #>>> dna.n_mismatch_substr_finder(10, 2)
        #['GCACACAGAC', 'GCGCACACAC']
        >>> Dna('CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT').n_mismatch_substr_finder(9, 3, 'n_mismatch_compl')
        ['AGCGCCGCT', 'AGCGGCGCT']
        >>> dna = Dna('ACGTTGCATGTCGCATGATGCATGAGAGCT').n_mismatch_substr_finder(4, 1, 'n_mismatch')
        ['GATG', 'ATGC', 'ATGT']
        >>> Dna('ACGTTGCATGTCGCATGATGCATGAGAGCT').n_mismatch_substr_finder(4, 1, 'n_mismatch_compl')
        ['ATGT', 'ACAT']
        """
        # Generating raw possible kmers
        results = {}

        get_comparator = {
          'n_mismatch': self.get_num_mismatch_comparator,
          'n_mismatch_compl': self.get_num_mismatch_complementary_comparator
        }.get(cmrtr_type, self.get_num_mismatch_comparator)

        for substr in self.n_substr_generator(k):
            comparator, count = get_comparator(substr, d), 1
            for kmer, data in results.iteritems():
                if data['comparator'](substr):
                    results[kmer]['count'] += 1
                if comparator(kmer):
                    count += 1
            if not results.has_key(substr):
                results[substr] = {'comparator': comparator, 'count': count}

        # selecting candidats
        counts = sorted(map(lambda (a, b): b['count'], results.iteritems()))
        res_treshold = int((counts[2]+counts[-2])/2)
        #res_treshold = 0
        prekmers = map(lambda (key, _): key, filter(lambda (_, val): val['count'] >= res_treshold, results.iteritems()))

        #kmers = set(self.complementary_filter(possible for prekmer in prekmers for possible in self.n_mismatch_generator(prekmer, d)))
        kmers = set(possible for prekmer in prekmers for possible in self.n_mismatch_generator(prekmer, d))
        candidates = dict(zip((kmer for kmer in kmers), map(lambda kmer: dict(
            comparator=get_comparator(kmer, d), count=0), kmers)))
        #counting
        for substr in self.n_substr_generator(k):
             for kmer, d in candidates.iteritems():
                if d['comparator'](substr):
                    candidates[kmer]['count'] += 1

        #select winners
        treshold = max(map(lambda (_, v): v['count'], candidates.iteritems()))
        print "Raw results: %i, min: %i, :max: %i, median: %i, candidates: %i, treshold_count: %i" % (len(counts), min(counts), max(counts), res_treshold, len(candidates), treshold)
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

    def get_num_mismatch_complementary_comparator(self, pattern, max_mismatches):
        """
        >>> c = Dna().get_num_mismatch_complementary_comparator('ATGT', 1)
        >>> [c(i) for i in ['ATGT', 'ACAT', 'AGGT', 'ATGA', 'TCAT']]
        [True, True, True, True, True]
        """
        comp_pattern = Dna(pattern).complementary().genome
        def comparator(test_string):
            errors = 0
            c_errors = 0
            for i in xrange(len(pattern)):
                if errors <= max_mismatches and test_string[i] != pattern[i]:
                    errors += 1
                if c_errors <= max_mismatches and test_string[i] != comp_pattern[i]:
                    c_errors += 1
                if errors > max_mismatches and c_errors > max_mismatches:
                    return False
            return True

        func = comparator
        func.pattern = pattern
        return func


    def n_substr_generator(self, length):
        for i in xrange(0, len(self.genome) - length + 1):
            yield self.genome[i:i + length]

    def n_mismatch_generator(self, substring, num_mismatches=1):
        """
        >>> d = Dna('ACTG')
        >>> list(d.n_mismatch_generator(d.genome, 2))
        ['CATG', 'CTTG', 'CGTG', 'TATG', 'TTTG', 'TGTG', 'GATG', 'GTTG', 'GGTG', 'CCAG', 'CCCG', 'CCGG', 'TCAG', 'TCCG', 'TCGG', 'GCAG', 'GCCG', 'GCGG', 'CCTA', 'CCTC', 'CCTT', 'TCTA', 'TCTC', 'TCTT', 'GCTA', 'GCTC', 'GCTT', 'AAAG', 'AACG', 'AAGG', 'ATAG', 'ATCG', 'ATGG', 'AGAG', 'AGCG', 'AGGG', 'AATA', 'AATC', 'AATT', 'ATTA', 'ATTC', 'ATTT', 'AGTA', 'AGTC', 'AGTT', 'ACAA', 'ACAC', 'ACAT', 'ACCA', 'ACCC', 'ACCT', 'ACGA', 'ACGC', 'ACGT']
        >>> list(d.n_mismatch_generator(d.genome, 1))
        ['CCTG', 'TCTG', 'GCTG', 'AATG', 'ATTG', 'AGTG', 'ACAG', 'ACCG', 'ACGG', 'ACTA', 'ACTC', 'ACTT']
        """
        for locs in itertools.combinations(range(len(substring)), num_mismatches):
            current_substr = [[char] for char in substring]
            for loc in locs:
                orig_char = substring[loc]
                current_substr[loc] = [l for l in self.alfabet if l != orig_char]
            for poss in itertools.product(*current_substr):
                yield ''.join(poss)

    def complementary_filter(self, sequence):
        """
        >>> list(Dna().complementary_filter(['CCTG', 'TCTG', 'GCTG', 'CAGG', 'CAGC']))
        ['CCTG', 'TCTG', 'GCTG']
        """
        heap = set()
        for genome in sequence:
            if genome in heap or Dna(genome).complementary().genome in heap:
                continue
            heap.add(genome)
            yield genome

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

    @classmethod
    def h_file_n_mismatch_substr(cls, path):
        data = open(path, 'r')
        fragment, k, d = data.readline().strip().split(' ')
        k = int(k)
        d = int(d)
        dna = Dna(fragment)
        data.close()
        print 'DNA: %s, k-length: %i, num_mismatch: %i' % (fragment, k, d)
        print ' '.join(map(lambda x: str(x), dna.n_mismatch_substr_finder(k, d)))

    @classmethod
    def h_file_n_mismatch_substr_compl(cls, path):
        data = open(path, 'r')
        fragment = Dna(data.readline().strip())
        k, d = map(lambda x: int(x), data.readline().strip().split())
        data.close()
        print 'DNA: %s, k-length: %i, num_mismatch: %i' % (fragment.genome, k, d)
        print ' '.join(map(lambda x: str(x), fragment.n_mismatch_substr_finder(k, d, 'n_mismatch_compl')))


