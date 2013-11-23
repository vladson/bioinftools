import rna
import dna

class Protein:

    def __init__(self, seq = ""):
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

    def __repr__(self):
        return self.seq

    def nalen(self):
        return len(self.seq)*3

    def possible_rnas(self, position=0, assembly=''):
        """
        >>> prt = Protein('MA')
        >>> list(prt.possible_rnas())
        [AUGGCA, AUGGCC, AUGGCG, AUGGCU]
        """
        for codon in self.reverse_codon_table[self.seq[position]]:
            if position + 1 == len(self.seq):
                yield rna.Rna(assembly + codon)
            else:
                for fragment in self.possible_rnas(position + 1, assembly + codon):
                    yield rna.Rna(fragment)

    def possible_dnas(self):
        for fragment in self.possible_rnas():
            yield fragment.to_dna()
            yield fragment.complementary().to_dna()


    @classmethod
    def check_file(cls, path):
        file = open(path, 'r')
        file.readline()
        dn = dna.Dna(file.readline().strip())
        pr = Protein(file.readline().strip())
        return dn, pr
        file.close()
        dnas = list(pr.possible_dnas())
        fragments = set()
        for fragment in dn.n_substr_generator(pr.nalen()):
            if dna.Dna(fragment) in dnas:
                fragments.add(fragment)

        for fragment in fragments:
            print fragment
