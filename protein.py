import numpy
import rna
import dna


class Protein:
    reverse_codon_table = {}
    table = open('RNA_codon_table_1.txt', 'r')
    for line in table:
        data = line.strip().split(' ')
        if len(data) > 1:
            if reverse_codon_table.has_key(data[1]):
                reverse_codon_table[data[1]].append(data[0])
            else:
                reverse_codon_table[data[1]] = [data[0]]
    table.close()
    integer_masses = {}
    integer_masses_uniq = {}
    integer_masses_reverse = {}
    table = open('integer_mass_table.txt', 'r')
    for line in table:
        aa, mass = line.strip().split(' ')
        integer_masses[aa] = int(mass)
        if not int(mass) in integer_masses_uniq.values():
            integer_masses_uniq[aa] = int(mass)
            integer_masses_reverse[int(mass)] = aa
    table.close()

    def __init__(self, seq="", seq_score=0):
        self.seq = seq
        self.seq_score = seq_score

    def __repr__(self):
        return self.seq

    def mass_repr(self):
        """
        >>> Protein('KWI').mass_repr()
        '128-186-113'
        """
        return '-'.join(map(lambda x: str(self.integer_masses[x]), self.seq))

    def nalength(self):
        # length of corresponding nucleic acid
        return len(self.seq) * 3

    def mass(self):
        """
        >>> Protein('').mass()
        0
        >>> Protein('LN').mass()
        227
        """
        mass = 0
        for aa in self.seq:
            mass += self.integer_masses[aa]
        return mass

    def spectrum_leaderboard_seq(self, spectrum, leaders=5, convolute=0):
        """
        >>> print Protein().spectrum_leaderboard_seq([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460], 10).mass_repr()
        71-147-113-129
        >>> print Protein().spectrum_leaderboard_seq([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 60, 20).mass_repr()
        57-137-71-99-129
        >>> Protein().restore_alfabet()
        """
        dbg = False
        if dbg:
            print "Spectrum len: %i, n: %i, max_mass: %i, convolution num: %i" % (
                len(spectrum), leaders, max(spectrum), convolute)
        if convolute:
            self.set_alfabet(self.convolution_alfabet(spectrum, convolute, 'existing_t'))
            if dbg:
                print Protein.integer_masses
        parentMasses = set(spectrum)
        scorer = self.get_spectrum_scorer(spectrum)
        leader = Protein()
        leader.seq_score = 1
        leaderboard = [Protein()]
        while len(leaderboard):
            new_leaderboard = []
            for candidate in self.expanded(leaderboard):
                if candidate.mass() <= max(parentMasses):
                    if scorer(candidate) > leader.seq_score:
                        leader = candidate
                    new_leaderboard.append(candidate)
            leaderboard = self.score_cutted(new_leaderboard, leaders)
            if dbg:
                print "score: %i, size: %i, trimmed_size: %i" % (
                    leader.seq_score, len(new_leaderboard), len(leaderboard))
        return leader


    def get_spectrum_scorer(self, pattern):
        """
        >>> scorer = Protein().get_spectrum_scorer([0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])
        >>> peptide = Protein('NQEL')
        >>> print scorer(peptide)
        14
        """

        def scorer(candidate):
            score_list = list(pattern)
            for mass in list(candidate.theoretical_ms()):
                if mass in score_list:
                    score_list.remove(mass)
            score = len(pattern) - len(score_list)
            candidate.seq_score = score
            return score

        return scorer


    def score_cutted(self, new_leaderboard=[], leaders=5):
        """
        >>> l = map(lambda i: Protein('K'*i,i), xrange(5))
        >>> Protein().score_cutted(l,3)
        [KKKK, KKK, KK]
        """
        if len(new_leaderboard) <= leaders:
            return new_leaderboard
        new_leaderboard = sorted(new_leaderboard, key=lambda c: c.seq_score, reverse=True)
        min_score = new_leaderboard[-1].seq_score if len(new_leaderboard) < leaders else new_leaderboard[
            leaders - 1].seq_score
        return filter(lambda p: p.seq_score >= (min_score), new_leaderboard)

    def spectrum_seq(self, spectrum):
        """
        >>> print ' '.join(map(lambda p: p.mass_repr(), Protein().spectrum_seq([0, 113, 128, 186, 241, 299, 314, 427]))[::-1])
        186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
        """
        checker = self.get_spectrum_consistence_checker(spectrum)
        candidates = [Protein()]
        while len(candidates):
            new_candidates = []
            for peptide in self.expanded(candidates):
                if checker(peptide.theoretical_ms(True)):
                    if sorted(list(peptide.theoretical_ms())) == sorted(spectrum):
                        yield peptide
                    else:
                        new_candidates.append(peptide)
            candidates = new_candidates

    def expanded(self, candidates):
        """
        >>> list(Protein().expanded([Protein()]))
        [A, C, E, D, G, F, I, H, K, M, N, P, S, R, T, W, V, Y]
        >>> map(lambda p: p.mass(), Protein().expanded([Protein()]))
        [71, 103, 129, 115, 57, 147, 113, 137, 128, 131, 114, 97, 87, 156, 101, 186, 99, 163]
        """
        for candidate in candidates:
            for aa in self.integer_masses_uniq:
                yield candidate.appended(aa)

    def convolute(self, spectrum=[]):
        """
        >>> list(Protein().convolute([0, 137, 186, 323]))
        [137, 186, 323, 49, 186, 137]
        """
        for traverse in xrange(0, len(spectrum)):
            i = spectrum[traverse]
            for j in spectrum[traverse:]:
                diff = abs(i - j)
                if diff > 0:
                    yield diff

    def convolution_alfabet(self, spectrum, num=20, type='concrete'):
        """
        >>> Protein().convolution_alfabet([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 20)
        [129, 137, 71, 99, 57, 194, 170, 186, 79, 91, 58, 95, 113, 115, 128, 136, 148, 151, 156, 157]
        >>> Protein().convolution_alfabet([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 20, 'tails')
        [129, 137, 71, 99, 57, 194, 170, 186, 79, 91, 58, 95, 113, 115, 128, 136, 148, 151, 156, 157, 162, 166, 171, 178, 65, 66, 72, 80, 87, 109, 123]
        >>> Protein().convolution_alfabet([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 20, 'existing_c')
        [129, 137, 71, 99, 57, 186, 113, 115, 128, 156]
        >>> Protein().convolution_alfabet([57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493], 20, 'existing_t')
        [129, 137, 71, 99, 57, 186, 113, 115, 128, 156, 87]
        """
        masset = {}
        for mass in self.convolute(spectrum):
            if 56 < mass < 200:
                if masset.has_key(mass):
                    masset[mass] += 1
                else:
                    masset[mass] = 1
        sorted_masses = sorted(masset.iteritems(), key=lambda (k, v): v, reverse=True)
        _, min_freq = sorted_masses[num - 1]
        return {
            'concrete': lambda: map(lambda (mass, freq): mass, sorted_masses[:num]),
            'tails': lambda: map(lambda (mass, freq): mass, filter(lambda (_, freq): freq >= min_freq, sorted_masses)),
            'existing_c': lambda: map(lambda (m, _): m, filter(lambda (m, _): m in self.integer_masses_reverse.keys(),
                                                               sorted_masses[:num])),
            'existing_t': lambda: map(lambda (m, _): m,
                                      filter(lambda (m, f): m in self.integer_masses_reverse.keys() and f >= min_freq,
                                             sorted_masses))
        }.get(type, lambda: map(lambda (mass, freq): mass, filter(lambda (_, freq): freq >= min_freq, sorted_masses)))()

    def get_spectrum_consistence_checker(self, pattern_spectrum):
        """
        >>> p = Protein('VPCHAMNIID')
        >>> c = p.get_spectrum_consistence_checker([0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093])
        >>> c(p.theoretical_ms())
        True
        """

        def check_spectrum(spectrum):
            for mass in spectrum:
                if mass not in set(pattern_spectrum):
                    return False
            return True

        return check_spectrum

    def theoretical_ms(self, cyclic=True):
        """
        >>> sorted(list(Protein('NQEL').theoretical_ms()))
        [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        >>> sorted(list(Protein('VPCHAMNIID').theoretical_ms()))
        [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
        """
        for pept in self.possible_subpeptides(cyclic):
            yield pept.mass()

    def possible_subpeptides(self, cyclic=True):
        """
        >>> list(Protein('NQEL').possible_subpeptides())
        [, N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, LNQ, ELN, NQEL]
        """
        yield Protein()
        for fraglen in xrange(1, len(self.seq) + 1):
            for offset in xrange(len(self.seq) - fraglen + 1):
                yield Protein(self.seq[offset:(offset + fraglen)])
            if cyclic and len(self.seq) > fraglen > 1:
                for neg in xrange(1, fraglen):
                    yield Protein(self.seq[-neg:] + self.seq[:(fraglen - neg)])

    def possible_rnas(self, position=0, assembly=''):
        """
        >>> prt = Protein('MA')
        >>> list(prt.possible_rnas())
        ['AUGGCA', 'AUGGCC', 'AUGGCG', 'AUGGCU']
        """
        for codon in self.reverse_codon_table[self.seq[position]]:
            if position + 1 == len(self.seq):
                yield assembly + codon
            else:
                for fragment in self.possible_rnas(position + 1, assembly + codon):
                    yield fragment

    def possible_dnas(self):
        for fragment in self.possible_rnas():
            fragment = rna.Rna(fragment)
            yield fragment.to_dna()
            yield fragment.complementary().to_dna()

    def set_alfabet(self, masset=[]):
        """
        Tests go well but ruines all
        >>> Protein().set_alfabet([10, 20, 40])
        >>> print Protein.integer_masses
        {'A': 10, 'C': 40, 'B': 20}
        >>> Protein().restore_alfabet()
        """
        self.__class__.old_integer_masses = self.integer_masses
        self.__class__.old_integer_masses_uniq = self.integer_masses_uniq
        self.__class__.old_integer_masses_reverse = self.integer_masses_reverse
        alfabet = {}
        alfabet_reverse = {}
        for i, mass in enumerate(sorted(masset)):
            char = chr((65 if i < 26 else 71) + i)
            alfabet[char] = mass
            alfabet_reverse[mass] = char
        self.integer_masses = alfabet
        self.__class__.integer_masses = alfabet
        self.integer_masses_reverse = alfabet_reverse
        self.__class__.integer_masses_reverse = alfabet_reverse
        self.integer_masses_uniq = alfabet
        self.__class__.integer_masses_uniq = alfabet

    def restore_alfabet(self):
        if self.old_integer_masses:
            self.integer_masses = self.old_integer_masses
            self.__class__.integer_masses = self.old_integer_masses
            self.integer_masses_reverse = self.old_integer_masses_reverse
            self.__class__.integer_masses_reverse = self.old_integer_masses_reverse
            self.integer_masses_uniq = self.old_integer_masses_uniq
            self.__class__.integer_masses_uniq = self.old_integer_masses_uniq


    def append(self, aa=''):
        """
        >>> p = Protein('L')
        >>> p.append('N')
        >>> print p
        LN
        """
        if aa in self.integer_masses.keys():
            self.seq += aa
        else:
            raise AttributeError

    def appended(self, aa=''):
        """
        >>> Protein('L').appended('N')
        LN
        """
        if aa in self.integer_masses.keys():
            return Protein(self.seq + aa)
        raise AttributeError


    @classmethod
    def from_mass_repr(cls, repr=''):
        """
        >>> Protein.from_mass_repr('128-186-113')
        KWI
        """
        return cls(''.join(map(lambda m: cls.integer_masses_reverse[int(m)], repr.split('-'))))


    @classmethod
    def check_file_getdnas(cls, path):
        file = open(path, 'r')
        dn = dna.Dna(file.readline().strip())
        pr = Protein(file.readline().strip())
        file.close()
        dnas = list(pr.possible_dnas())
        fragments = set()
        for fragment in dn.kmer_generator(pr.nalength()):
            if dna.Dna(fragment) in dnas:
                fragments.add(fragment)

        for fragment in fragments:
            print fragment


    @classmethod
    def leaderboard_seq(cls, path, test=False):
        file = open(path, 'r')
        if test:
            file.readline()
        n = int(file.readline().strip())
        l = map(lambda i: int(i), file.readline().strip().split(' '))
        print "Spectrum len: %i, n: %i, max_mass: %i" % (len(l), n, max(l))
        peptide = Protein().spectrum_leaderboard_seq(l, n)
        print peptide, peptide.mass_repr()
        if test:
            file.readline()
            test = Protein.from_mass_repr(file.readline())
            print test, test.mass_repr()

    @classmethod
    def leaderboard_seq_convol(cls, path, test=False):
        file = open(path, 'r')
        if test:
            file.readline()
        m = int(file.readline().strip())
        n = int(file.readline().strip())
        l = map(lambda i: int(i), file.readline().strip().split(' '))
        print "Spectrum len: %i, n: %i, max_mass: %i, convolve to: %i" % (len(l), n, max(l), m)
        peptide = Protein().spectrum_leaderboard_seq(l, n, m)
        print peptide, peptide.mass_repr()
        if test:
            file.readline()
            test = Protein.from_mass_repr(file.readline())
            print test, test.mass_repr()
        Protein().restore_alfabet()