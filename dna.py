#! /usr/bin/env python
import itertools
import random
import operator
from scipy import stats

__author__ = 'vladson'


class Dna:
    conversion_table = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    alfabet = 'A', 'C', 'T', 'G'
    test = True

    def __init__(self, str=''):
        self.genome = str.upper()
        self.hamming_score = len(self.genome) * 100

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

    #
    #   Motiff Problems
    #

    @classmethod
    def motiff_generator(cls, dnas=[], k=3, d=1):
        """
        dnas = list of DNA strings
        k = kmer length
        d = max mismatch number
        >>> print list(Dna.motiff_generator(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]))
        Existing: 18, seen: 64
        ['ATA', 'ATT', 'TTT', 'GTT']
        """
        motiff_checker = Dna.get_dnas_d_motiff_checker(dnas, d)
        existing = set(kmer for genome in dnas for kmer in Dna(genome).kmer_generator(k))

        ## Strategy one
        #possibles = set(possible for kmer in existing for possible in Dna().n_mismatch_generator(kmer, d))
        #for possible in possibles:
        #   if motiff_checker(possible):
        #        yield possible

        # Strategy two
        seen = set()
        for kmer in existing:
            for possible in Dna().n_mismatch_generator(kmer, k):
                if possible in seen:
                    continue
                seen.add(possible)
                if motiff_checker(possible):
                    yield possible
        if Dna.test:
            print "Existing: %i, seen: %i" % (len(existing), len(seen))

    @staticmethod
    def median_string(dnas=[], length=3):
        """
        #>>> dnas = ["TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT", "CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA", "TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT", "TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA", "ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG", "TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA", "TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC", "GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA", "CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG", "CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG"]
        #>>> Dna.median_string(dnas, 6)
        CGCCGA
        >>> l = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTACGGGACAG"]
        >>> Dna.median_string(l)
        ACG
        """
        leader = Dna('A' * length)
        for kmer in Dna.all_kmer_generator(length):
            kmer = Dna(kmer)
            if Dna.hd_list(kmer, dnas) < leader.hamming_score:
                leader = kmer
        return leader

    @staticmethod
    def greedy_motiff_search(motiffs, k, t, pseudocounts=True):
        """
        >>> Dna.greedy_motiff_search(["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5, False)
        ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
        >>> Dna.greedy_motiff_search(["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"], 3, 5)
        ['TTC', 'ATC', 'TTC', 'ATC', 'TTC']
        """
        best = map(lambda kmer: kmer[0:k], motiffs)
        best_score = Dna.motiff_scorer(best)
        for offset in xrange(1, len(motiffs[0]) - k):
            candidates = [motiffs[0][offset:offset + k]]
            for datum in motiffs[1:]:
                candidates.append(Dna(datum).most_probable_motiff(Profile.from_motiffs(candidates, not pseudocounts)))
            candidates_score = Dna.motiff_scorer(candidates)
            if candidates_score < best_score:
                best, best_score = candidates, candidates_score
        return best

    @staticmethod
    def iterative_randomized_motiff_search(dnas, k, t, treshold=10, out=False):
        """
        This algorithm is Monte-Carlo so output is not guarantied
        #>>> dnas = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        #>>> Dna.iterative_randomized_motiff_search(dnas, 8, 5, treshold=9)
        ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
        """
        try:
            itercount = 0
            while itercount < 1000:
                itercount += 1
                best_motifs, inneritercount, best_score = Dna.randomized_motiff_search(dnas, k, t)
                if best_score <= treshold:
                    if out:
                        print "Reached %i on iteration %i with score %i (on iter %i)" % (
                            treshold, itercount, best_score, inneritercount)
                    return best_motifs
        except KeyboardInterrupt:
            print "Interrupted on iteration %i with score %i (on iter %i)" % (itercount, best_score, inneritercount)
            for motiff in best_motifs:
                print motiff
            raise

    @staticmethod
    def randomized_motiff_search(dnas, k, t, out=False):
        try:
            best_motifs = motifs = [Dna(line).rand_kmer(k) for line in dnas]
            best_score = Dna.motiff_scorer(best_motifs, False)
            itercount = 0
            while True:
                itercount += 1
                profile = Profile.from_motiffs(motifs, False)
                motifs = Dna.profile_motiffs(dnas, profile)
                score = Dna.motiff_scorer(motifs, False)
                if score < best_score:
                    best_motifs, best_score = motifs, score
                else:
                    if out:
                        print "Finished on iteration %i with score %i" % (itercount, best_score)
                    return best_motifs, itercount, best_score
        except KeyboardInterrupt:
            print "Interrupted on iteration %i with score %i" % (itercount, best_score)
            raise

    @staticmethod
    def iterative_gibbs_sampler(dnas, k, t, n, treshold=10, out=False):
        """
        Monte-Carlo hard to test
        #>>> dnas = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
        #>>> Dna.iterative_gibbs_sampler(dnas, 8, 5, 100, treshold=9, out=False)
        ['TGTTCAGT', 'GCCAAGGT', 'AGTACCGA', 'ATCAAGTT', 'ACCAGCTC']
        """
        try:
            itercount = 0
            best_motifs, best_score = [], 10000
            while itercount < 1000:
                itercount += 1
                motifs, score = Dna.gibbs_sampler(dnas, k, t, n)
                if score < best_score:
                    best_motifs, best_score = motifs, score
                if best_score <= treshold:
                    if out:
                        print "Reached %i on iteration %i with score %i " % (treshold, itercount, best_score)
                    return best_motifs
        except KeyboardInterrupt:
            print "Interrupted on iteration %i with min score %i" % (itercount, best_score)
            for motiff in best_motifs:
                print motiff
            raise

    @staticmethod
    def gibbs_sampler(dnas, k, t, n, out=False):
        try:
            best_motifs = motifs = [Dna(line).rand_kmer(k) for line in dnas]
            best_score = Dna.motiff_scorer(best_motifs, False)
            itercount = 0
            for i in xrange(n):
                itercount += 1
                i = random.randint(0, t - 1)
                motifs.pop(i)
                profile = Profile.from_motiffs(motifs, False)
                motifs.insert(i, profile.draw_kmer(Dna(dnas[i]).kmer_generator(k)))
                score = Dna.motiff_scorer(motifs, False)
                if score < best_score:
                    best_motifs, best_score = motifs, score
            if out:
                print "Finished with score %i" % best_score
            return best_motifs, best_score
        except KeyboardInterrupt:
            print "Interrupted on iteration %i with score %i" % (itercount, best_score)
            raise

    @staticmethod
    def profile_motiffs(dnas, profile):
        """
        >>> pr = Profile.from_lines_legend(['A C G T', "0.2 0.4 0.3 0.1", "0.2 0.3 0.3 0.2", "0.3 0.1 0.5 0.1", "0.2 0.5 0.2 0.1", "0.3 0.1 0.4 0.2"])
        >>> dnas = ['ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGTGGGCCT']
        >>> Dna.profile_motiffs(dnas, pr)
        ['CCGAG', 'CCGAA']
        """
        return [Dna(line).most_probable_motiff(profile) for line in dnas]

    def most_probable_motiff(self, profile):
        """
        >>> prof = Profile([{'A': 0.357, 'C': 0.357, 'T': 0.107, 'G': 0.179}, {'A': 0.393, 'C': 0.179, 'T': 0.286, 'G': 0.143}, {'A': 0.179, 'C': 0.179, 'T': 0.321, 'G': 0.321}, {'A': 0.214, 'C': 0.179, 'T': 0.286, 'G': 0.321}, {'A': 0.286, 'C': 0.25, 'T': 0.214, 'G': 0.25}, {'A': 0.286, 'C': 0.25, 'T': 0.286, 'G': 0.179}, {'A': 0.393, 'C': 0.143, 'T': 0.25, 'G': 0.214}])
        >>> d = Dna('GGTATGCGCACTTCCGAAGAAGGATGCTCAATCATACAAGACACATTCCATCGAGGTAGTTTGACTGGCGAAGTCCCGACTCGCTCACAACTAGTATCCTGTGAAGTCCAGCGTTGAACGACGTGTTGGCTTTAAGCGCCCTGCTTTTCACCAGTTTCTCTCCTAAGTTCGTTCCAGGTCCAAACTGTGGCACTGCAAAT')
        >>> d.most_probable_motiff(prof)
        'CATTCCA'
        >>> prof.pr('CATTCCA')
        0.000316376632947375
        >>> d = Dna('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT')
        >>> pr = Profile.from_lines_legend(['A C G T', "0.2 0.4 0.3 0.1", "0.2 0.3 0.3 0.2", "0.3 0.1 0.5 0.1", "0.2 0.5 0.2 0.1", "0.3 0.1 0.4 0.2"])
        >>> d.most_probable_motiff(pr)
        'CCGAG'
        >>> pr.pr(pr.consensus())
        0.012
        """
        generator = self.kmer_generator(profile.length())
        leader = generator.next()
        probability = profile.pr(leader)
        for kmer in generator:
            if profile.pr(kmer) > probability:
                leader, probability = kmer, profile.pr(kmer)
        return leader

    def substr_finder(self, length=3):
        results = {}
        for substr in self.kmer_generator(length):
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
        Frequent words with mismatches and complementary
        k is length of k-mer
        d is number of possible mismatches
        #>>> dna = Dna('CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC')
        #>>> dna.n_mismatch_substr_finder(10, 2)
        #['GCACACAGAC', 'GCGCACACAC']
        #>>> Dna('CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT').n_mismatch_substr_finder(9, 3, 'n_mismatch_compl')
        ['AGCGCCGCT', 'AGCGGCGCT']
        #>>> dna = Dna('ACGTTGCATGTCGCATGATGCATGAGAGCT').n_mismatch_substr_finder(4, 1, 'n_mismatch')
        ['GATG', 'ATGC', 'ATGT']
        #>>> Dna('ACGTTGCATGTCGCATGATGCATGAGAGCT').n_mismatch_substr_finder(4, 1, 'n_mismatch_compl')
        ['ATGT', 'ACAT']
        """
        # Generating raw possible kmers
        results = {}

        get_comparator = {
            'n_mismatch': self.get_num_mismatch_comparator,
            'n_mismatch_compl': self.get_num_mismatch_complementary_comparator
        }.get(cmrtr_type, self.get_num_mismatch_comparator)

        for substr in self.kmer_generator(k):
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
        res_treshold = int((counts[2] + counts[-2]) / 2)
        #res_treshold = 0
        prekmers = map(lambda (key, _): key, filter(lambda (_, val): val['count'] >= res_treshold, results.iteritems()))

        #kmers = set(self.complementary_filter(possible for prekmer in prekmers for possible in self.n_mismatch_generator(prekmer, d)))
        kmers = set(possible for prekmer in prekmers for possible in self.n_mismatch_generator(prekmer, d))
        candidates = dict(zip((kmer for kmer in kmers), map(lambda kmer: dict(
            comparator=get_comparator(kmer, d), count=0), kmers)))
        #counting
        for substr in self.kmer_generator(k):
            for kmer, d in candidates.iteritems():
                if d['comparator'](substr):
                    candidates[kmer]['count'] += 1

        #select winners
        treshold = max(map(lambda (_, v): v['count'], candidates.iteritems()))
        print "Raw results: %i, min: %i, :max: %i, median: %i, candidates: %i, treshold_count: %i" % (
            len(counts), min(counts), max(counts), res_treshold, len(candidates), treshold)
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

    def at_position(self, offset=0, length=3):
        """
        >>> print Dna('CGCCCGAATCCAGAACGCATTCCCATATTTCGG').at_position(5, 8)
        GAATCCAG
        """
        return self.genome[offset:offset + length]

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

    #
    # Service functions
    #


    @classmethod
    def get_random_string(cls, length):
        """
        >>> len(Dna.get_random_string(5).genome)
        5
        """
        return Dna(''.join(random.choice(cls.alfabet) for x in range(length)))

    # motiff scorer
    @staticmethod
    def motiff_scorer(motiffs, real=True):
        """
        >>> Dna.motiff_scorer(['TCTGAGCTTGCGTTATTTTTAGACC', 'GTTTGACGGGAACCCGACGCCTATA', 'TTTTAGATTTCCTCAGTCCACTATA', 'CTTACAATTTCGTTATTTATCTAAT', 'CAGTAGGAATAGCCACTTTGTTGTA', 'AAATCCATTAAGGAAAGACGACCGT'])
        69
        """
        profile = Profile.from_motiffs(motiffs, real)
        return Dna.hd_list(profile.consensus(), motiffs)

    # Hamming distance related
    @staticmethod
    def hamming(seq1, seq2):
        """
        Hamming distance
        >>> Dna.hamming('ACCCTG', 'ACCTTG')
        1
        """
        if isinstance(seq1, Dna):
            seq1 = seq1.genome
        if isinstance(seq2, Dna):
            seq2 = seq2.genome
        assert len(seq1) == len(seq2)
        return sum(itertools.imap(operator.ne, seq1, seq2))

    @staticmethod
    def hd_long(pattern, genome, indx=False):
        """
        >>> Dna.hd_long('AAA', 'TTACCTTAAC')
        1
        """
        if not isinstance(genome, Dna):
            genome = Dna(genome)
        hamming_offsets = []
        for possible in genome.kmer_generator(len(pattern)):
            hamming_offsets.append(Dna.hamming(pattern, possible))
        if indx:
            return hamming_offsets.index(min(hamming_offsets))
        else:
            return min(hamming_offsets)

    @staticmethod
    def hd_list(pattern, dnas):
        """
        >>> l =["TTACCTTAAC", 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
        >>> Dna.hd_list('AAA', l)
        5
        >>> p = Dna('AAA')
        >>> print p.hamming_score
        300
        >>> Dna.hd_list(p, l)
        5
        >>> print p.hamming_score
        5
        """
        if isinstance(pattern, Dna):
            pattern.hamming_score = sum(map(Dna.hd_long, [pattern.genome] * len(dnas), dnas))
            return pattern.hamming_score
        else:
            return sum(map(Dna.hd_long, [pattern] * len(dnas), dnas))

    @staticmethod
    def motiff(pattern, genome):
        """
        >>> Dna.motiff('GATTCTCA', 'GCAAAGACGCTGACCAA')
        'GACGCTGA'
        """
        if not isinstance(genome, Dna):
            genome = Dna(genome)
        return genome.at_position(Dna.hd_long(pattern, genome, True), len(pattern))


    @classmethod
    def get_dnas_d_motiff_checker(cls, dnas, d):
        """
        list checker if kmer is present in all genomes
        #>>> ck =Dna.get_dnas_d_motiff_checker(['TCTGAGCTTGCGTTATTTTTAGACC', 'GTTTGACGGGAACCCGACGCCTATA', 'TTTTAGATTTCCTCAGTCCACTATA', 'CTTACAATTTCGTTATTTATCTAAT', 'CAGTAGGAATAGCCACTTTGTTGTA', 'AAATCCATTAAGGAAAGACGACCGT'], 2)
        #>>> ck('GCCTT')
        True
        >>> c = Dna.get_dnas_d_motiff_checker(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 1)
        >>> c('ATA')
        True
        >>> c('GGG')
        False
        """

        def checker(kmer):
            dnas_objs = dict(zip(dnas, [False] * len(dnas)))
            for mutation in Dna().n_mismatch_generator(kmer, d, True):
                for (genome, val) in dnas_objs.iteritems():
                    if not val:
                        if mutation in genome:
                            dnas_objs[genome] = True
                if reduce(lambda a, b: a and b, dnas_objs.values(), True):
                    return True
            return False

        return checker

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

    #
    # Assemble related
    #

    def kmer_composition(self, k):
        """
        >>> Dna('CAATCCAAC').kmer_composition(5)
        ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
        """
        return sorted([kmer for kmer in self.kmer_generator(k)])

    def kdmer_composition(self, k, d):
        """
        >>> Dna('TAATGCCATGGGATGTT').kdmer_composition(3,1)
        ['CAT|ATG', 'CCA|AAT', 'GAT|ATG', 'GGA|CAT', 'GGG|CCA', 'GTT|GGA', 'TAA|GCC', 'TGC|ATG', 'TGG|ATG', 'TGG|GCC', 'TGT|GGG']
        """
        return sorted(map(lambda kdmer: '|'.join(kdmer), self.kdmer_generator(k, d)))

    @staticmethod
    def kdmer_splitter(kdmer):
        """
        >>> Dna.kdmer_splitter('GAGA|TTGA')
        ['GAG|TTG', 'AGA|TGA']
        """
        return map(lambda kd: '|'.join(kd), zip(*map(lambda kmer: [kmer[:-1], kmer[1:]], kdmer.split('|'))))

    @staticmethod
    def get_kdmer_nodes_overlap_assembler(d, real_k=False):
        """
        #>>> Dna.get_kdmer_nodes_overlap_assembler(1)(["AG|AG", "GC|GC", "CA|CT", "AG|TG", "GC|GC", "CT|CT", "TG|TG", "GC|GC", "CT|CA"])
        AGCAGCTGCTGCA
        """
        def spl(kdmer):
            return kdmer.split('|')

        def assembler(iterable):
            first = True
            for node in iter(iterable):
                a, b = spl(node.key)
                #a, b = spl(node)
                if first:
                    if real_k:
                        k = len(a)
                    else:
                        k = len(a)+1
                    a_string, b_string = a, b
                    first = False
                else:
                    a_string += a[-1]
                    b_string += b[-1]
            return a_string + b_string[-(k+d):]
        return assembler

    @staticmethod
    def nodes_overlap_assembler(iterable):
        result = ""
        first = True
        for node in iter(iterable):
            if first:
                result += node.key
                first = False
            else:
                result += node.key[-1]
        return str(result)

    @staticmethod
    def nodes_cyrcular_assembler(iterable):
        result = ""
        first = True
        for node in iter(iterable):
            if first:
                first = False
                continue
            else:
                result += node.key[-1]
        return str(result)
    # End of Assemble related

    def kdmer_generator(self, k, d):
        """
        >>> list(Dna('TAATGCCATGGGATGTT').kdmer_generator(3,1))
        [['TAA', 'GCC'], ['CCA', 'AAT'], ['CAT', 'ATG'], ['TGC', 'ATG'], ['TGG', 'GCC'], ['GGG', 'CCA'], ['GGA', 'CAT'], ['GAT', 'ATG'], ['TGG', 'ATG'], ['TGT', 'GGG'], ['GTT', 'GGA']]
        """
        for kmer in self.kmer_generator(2*k+d):
            yield sorted([kmer[:k], kmer[-k:]], reverse = True)

    def kmer_generator(self, k):
        """
        Important method: make all the kmer strings
        """
        for i in xrange(0, len(self.genome) - k + 1):
            yield self.genome[i:i + k]

    def rand_kmer(self, k):
        offset = random.randint(0, len(self.genome) - k)
        return self.genome[offset:offset + k]

    @classmethod
    def all_kmer_generator(cls, length=4):
        """
        >>> len(list(Dna.all_kmer_generator(3)))
        64
        """
        for kmer in itertools.product(cls.alfabet, repeat=length):
            yield ''.join(kmer)

    @staticmethod
    def n_mismatch_generator(substring, num_mismatches=1, eager=False):
        """
        >>> d = Dna('ACTG')
        >>> list(d.n_mismatch_generator(d.genome, 2))
        ['CATG', 'CTTG', 'CGTG', 'TATG', 'TTTG', 'TGTG', 'GATG', 'GTTG', 'GGTG', 'CCAG', 'CCCG', 'CCGG', 'TCAG', 'TCCG', 'TCGG', 'GCAG', 'GCCG', 'GCGG', 'CCTA', 'CCTC', 'CCTT', 'TCTA', 'TCTC', 'TCTT', 'GCTA', 'GCTC', 'GCTT', 'AAAG', 'AACG', 'AAGG', 'ATAG', 'ATCG', 'ATGG', 'AGAG', 'AGCG', 'AGGG', 'AATA', 'AATC', 'AATT', 'ATTA', 'ATTC', 'ATTT', 'AGTA', 'AGTC', 'AGTT', 'ACAA', 'ACAC', 'ACAT', 'ACCA', 'ACCC', 'ACCT', 'ACGA', 'ACGC', 'ACGT']
        >>> list(d.n_mismatch_generator(d.genome, 1))
        ['CCTG', 'TCTG', 'GCTG', 'AATG', 'ATTG', 'AGTG', 'ACAG', 'ACCG', 'ACGG', 'ACTA', 'ACTC', 'ACTT']
        >>> list(d.n_mismatch_generator(d.genome, 1, True))
        ['ACTG', 'CCTG', 'TCTG', 'GCTG', 'AATG', 'ATTG', 'AGTG', 'ACAG', 'ACCG', 'ACGG', 'ACTA', 'ACTC', 'ACTT']
        """
        if eager:
            yield substring
        for locs in itertools.combinations(range(len(substring)), num_mismatches):
            current_substr = [[char] for char in substring]
            for loc in locs:
                orig_char = substring[loc]
                current_substr[loc] = [l for l in Dna.alfabet if l != orig_char]
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

    #
    # Human output
    # TODO: remove or move to Stepic
    #
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

    @classmethod
    def h_file_motiff_generator(cls, path, test=False):
        data = open(path, 'r')
        if test:
            data.readline()
        k, d = map(lambda i: int(i), data.readline().strip().split())
        dnas = []
        if test:
            for i in range(6):
                dnas.append(data.readline().strip())
            data.readline()
            answer = set(data.readline().strip().split())
        else:
            for line in data.readlines():
                dnas.append(line.strip())
            data.close()
        print "Dnas:"
        print dnas
        print 'Begin'
        results = set(Dna.motiff_generator(dnas, k, d))
        if test:
            print "Done needed: %i, generated: %i" % (len(answer), len(results))
            print answer - results
        else:
            print ' '.join(results)

    @staticmethod
    def h_file_median_string(path):
        data = open(path, 'r')
        k = int(data.readline().strip())
        dnas = []
        for line in data.readlines():
            dnas.append(line.strip())
        data.close()
        median_string = Dna.median_string(dnas, k)
        print median_string, median_string.hamming_score


class Profile:
    legend = ['A', 'C', 'T', 'G']

    def __init__(self, profile=[], legend=['A', 'C', 'T', 'G']):
        self.matrix = profile
        self.legend = legend

    def __repr__(self):
        return self.matrix.__repr__()

    def __len__(self):
        return len(self.matrix)

    @classmethod
    def from_lines_legend(cls, iterable):
        """
        >>> lines = ['A C G T', "0.2 0.4 0.3 0.1", "0.2 0.3 0.3 0.2", "0.3 0.1 0.5 0.1", "0.2 0.5 0.2 0.1", "0.3 0.1 0.4 0.2"]
        >>> Profile.from_lines_legend(lines)
        [{'A': 0.2, 'C': 0.4, 'T': 0.1, 'G': 0.3}, {'A': 0.2, 'C': 0.3, 'T': 0.2, 'G': 0.3}, {'A': 0.3, 'C': 0.1, 'T': 0.1, 'G': 0.5}, {'A': 0.2, 'C': 0.5, 'T': 0.1, 'G': 0.2}, {'A': 0.3, 'C': 0.1, 'T': 0.2, 'G': 0.4}]
        >>> len(Profile.from_lines_legend(lines))
        5
        """
        iterable = iter(iterable)
        legend = iterable.next().strip().split()
        profile = []
        for line in iterable:
            profile.append(dict(zip(legend, map(lambda x: float(x), line.strip().split()))))
        return cls(profile, legend)

    @classmethod
    def from_motiffs(cls, motiffs, real=True):
        """
        >>> pr = Profile.from_motiffs(["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"])
        >>> pr.consensus()
        'TCGGGGATTTCC'
        """
        return cls(cls.countmatrix_to_probability(cls.count_matrix(motiffs, real)))

    def length(self):
        return len(self.matrix)

    def pr(self, kmer):
        """
        >>> lines = ['A C G T', "0.2 0.1 0 0.7", '0.2 0.6 0 0.2', '0 0 1 0', '0 0 1 0', '0 0 0.9 0.1', '0 0 0.9 0.1', '0.9 0 0.1 0', '0.1 0.4 0 0.5', '0.1 0.1 0 0.8', '0.1 0.2 0 0.7', '0.3 0.4 0 0.3', '0 0.6 0 0.4']
        >>> pr = Profile.from_lines_legend(lines)
        >>> len(pr)
        12
        >>> pr.pr('TCGGGGATTTCC')
        0.020575296
        >>> pr.pr('TCGCGGATTTCC')
        0.0
        """
        assert len(kmer) == len(self)
        return reduce(lambda s, p: s * p, map(lambda (i, na): self.matrix[i][na], enumerate(kmer)))

    def draw_kmer(self, kmer_generator):
        """
        Believe me - it works
        #>>> pr = Profile.from_motiffs(["TCGCCTGTTTTT", "CCGGTGACTTAC", "ACGGGGATTCTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"], False)
        #>>> pr.draw_kmer(Dna('TGACTTCGCTGGAAAACCGTA').kmer_generator(12))
        'TCGCTGGAAAAC'
        """
        kmers = [kmer for kmer in iter(kmer_generator)]
        raw_probabilities = [self.pr(kmer) for kmer in kmers]
        # normalise probabilities
        sum_probailities = sum(raw_probabilities)
        probabilities = map(lambda p: p / sum_probailities, raw_probabilities)
        distr = stats.rv_discrete(values=(range(len(kmers)), probabilities))
        draw_index = distr.rvs()
        return kmers[draw_index]

    def consensus(self):
        """
        >>> lines = ['A C G T', "0.2 0.1 0 0.7", '0.2 0.6 0 0.2', '0 0 1 0', '0 0 1 0', '0 0 0.9 0.1', '0 0 0.9 0.1', '0.9 0 0.1 0', '0.1 0.4 0 0.5', '0.1 0.1 0 0.8', '0.1 0.2 0 0.7', '0.3 0.4 0 0.3', '0 0.6 0 0.4']
        >>> pr = Profile.from_lines_legend(lines)
        >>> pr.consensus()
        'TCGGGGATTTCC'
        """
        return ''.join(map(lambda pd: max(pd.iteritems(), key=operator.itemgetter(1))[0], self.matrix))

    @staticmethod
    def count_matrix(motiffs, real=True):
        """
        >>> Profile.count_matrix( ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"])
        [{'A': 2, 'C': 1, 'T': 7, 'G': 0}, {'A': 2, 'C': 6, 'T': 2, 'G': 0}, {'A': 0, 'C': 0, 'T': 0, 'G': 10}, {'A': 0, 'C': 0, 'T': 0, 'G': 10}, {'A': 0, 'C': 0, 'T': 1, 'G': 9}, {'A': 0, 'C': 0, 'T': 1, 'G': 9}, {'A': 9, 'C': 0, 'T': 0, 'G': 1}, {'A': 1, 'C': 4, 'T': 5, 'G': 0}, {'A': 1, 'C': 1, 'T': 8, 'G': 0}, {'A': 1, 'C': 2, 'T': 7, 'G': 0}, {'A': 3, 'C': 4, 'T': 3, 'G': 0}, {'A': 0, 'C': 6, 'T': 4, 'G': 0}]
        """
        matrix = []
        scaffold = {'A': 0, 'C': 0, 'T': 0, 'G': 0} if real else {'A': 1, 'C': 1, 'T': 1, 'G': 1}
        for motiff in motiffs:
            for i, na in enumerate(motiff):
                if len(matrix) == i:
                    matrix.append(dict(scaffold))
                matrix[i][na] += 1
        return matrix

    @staticmethod
    def countmatrix_to_probability(matrix):
        """
         >>> counts = Profile.count_matrix( ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC", "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"])
         >>> Profile.countmatrix_to_probability(counts)
         [{'A': 0.2, 'C': 0.1, 'T': 0.7, 'G': 0.0}, {'A': 0.2, 'C': 0.6, 'T': 0.2, 'G': 0.0}, {'A': 0.0, 'C': 0.0, 'T': 0.0, 'G': 1.0}, {'A': 0.0, 'C': 0.0, 'T': 0.0, 'G': 1.0}, {'A': 0.0, 'C': 0.0, 'T': 0.1, 'G': 0.9}, {'A': 0.0, 'C': 0.0, 'T': 0.1, 'G': 0.9}, {'A': 0.9, 'C': 0.0, 'T': 0.0, 'G': 0.1}, {'A': 0.1, 'C': 0.4, 'T': 0.5, 'G': 0.0}, {'A': 0.1, 'C': 0.1, 'T': 0.8, 'G': 0.0}, {'A': 0.1, 'C': 0.2, 'T': 0.7, 'G': 0.0}, {'A': 0.3, 'C': 0.4, 'T': 0.3, 'G': 0.0}, {'A': 0.0, 'C': 0.6, 'T': 0.4, 'G': 0.0}]
        """
        probabilities = list(matrix)
        counts = map(lambda x: sum(x.values()), matrix)
        for i, count in enumerate(counts):
            for na, nacount in probabilities[i].iteritems():
                probabilities[i][na] = float(nacount) / count
        return probabilities
