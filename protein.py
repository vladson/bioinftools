import rna
import dna


class Protein:
    def __init__(self, seq=""):
        self.seq = seq
        self.reverse_codon_table = {}
        table = open('RNA_codon_table_1.txt', 'r')
        for line in table:
            data = line.strip().split(' ')
            if len(data) > 1:
                if self.reverse_codon_table.has_key(data[1]):
                    self.reverse_codon_table[data[1]].append(data[0])
                else:
                    self.reverse_codon_table[data[1]] = [data[0]]
        table.close()
        self.integer_masses = {}
        table = open('integer_mass_table.txt', 'r')
        for line in table:
            aa, mass = line.strip().split(' ')
            self.integer_masses[aa] = int(mass)
        table.close()

    def __repr__(self):
        return self.seq

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

    def theoretical_ms(self):
        """
        >>> sorted(list(Protein('NQEL').theoretical_ms()))
        [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        >>> sorted(list(Protein('IAQMLFYCKVATN').theoretical_ms()))
        [0, 71, 71, 99, 101, 103, 113, 113, 114, 128, 128, 131, 147, 163, 170, 172, 184, 199, 215, 227, 227, 231, 244, 259, 260, 266, 271, 286, 298, 298, 310, 312, 328, 330, 330, 372, 385, 391, 394, 399, 399, 399, 401, 413, 423, 426, 443, 443, 470, 493, 498, 502, 513, 519, 526, 527, 541, 554, 556, 557, 564, 569, 590, 598, 616, 626, 640, 654, 657, 658, 665, 670, 682, 697, 697, 703, 711, 729, 729, 753, 753, 771, 779, 785, 785, 800, 812, 817, 824, 825, 828, 842, 856, 866, 884, 892, 913, 918, 925, 926, 928, 941, 955, 956, 963, 969, 980, 984, 989, 1012, 1039, 1039, 1056, 1059, 1069, 1081, 1083, 1083, 1083, 1088, 1091, 1097, 1110, 1152, 1152, 1154, 1170, 1172, 1184, 1184, 1196, 1211, 1216, 1222, 1223, 1238, 1251, 1255, 1255, 1267, 1283, 1298, 1310, 1312, 1319, 1335, 1351, 1354, 1354, 1368, 1369, 1369, 1379, 1381, 1383, 1411, 1411, 1482]
        """
        for pept in self.possible_subpeptides():
            yield pept.mass()

    def possible_subpeptides(self, cyclic=True):
        """
        >>> list(Protein('NQEL').possible_subpeptides())
        [, L, N, Q, E, LN, NQ, EL, QE, LNQ, ELN, QEL, NQE, NQEL]
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
