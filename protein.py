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

    def __init__(self, seq="",seq_score=0):
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

    def spectrum_leaderboard_seq(self, spectrum, leaders=5):
        """
        >>> print Protein().spectrum_leaderboard_seq([0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460], 10).mass_repr()
        71-147-113-129
        """
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
            #print "score: %i, size: %i, trimmed_size: %i" % (leader.seq_score, len(new_leaderboard), len(leaderboard))
        return leader


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

    def get_spectrum_scorer(self, pattern):
        def scorer(candidate):
            score_list = list(pattern)
            #score = 0
            #score_list = set(pattern)
            for mass in set(candidate.theoretical_ms()):
                if mass in score_list:
                    score_list.remove(mass)
                    #score += 1
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
        min_score = new_leaderboard[-1].seq_score if len(new_leaderboard) < leaders else new_leaderboard[leaders-1].seq_score
        return filter(lambda p: p.seq_score >= (min_score), new_leaderboard)

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
        for fragment in dn.n_substr_generator(pr.nalength()):
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