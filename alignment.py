import re
import dynamic

class SequenceAlignment(dynamic.Dag):

    """
    # Correct but long
    #>>> ddg = SequenceAlignment('GTGTTTATGTTACCTTGTAGGGATGCAGCCATCGAAGCAATAGCGCCCTGAGCAACCCTACTAGCGAACTGTCTAGGTTCATAATTACCACTCTAGTCTACCCCACGCTTCGGTTGCACTCAATATCAGCGCTCGTTTAACTAAGCCCTCCATTACAATCACCGAATAATACAGGAACTCGTCGTTGCCGCCCTATGTTTGGAGACGCAAATGGCTACTGCGAAGACGGAAGTAGCGTTGTGTACCGGCTACTTCTGAACAGTCTACAGGTCATCCTGTTTCTCTATTTATAGTCGCCAGAGTTCCAATACATCAGCCACATCACAAATCCGGGGGCCCTCTTACGAACCATGACAAACTGAATCACTTGCGAAACCTCTACCCCACAGGACCCCGCCGATCAGCGTAGATTTACGTCGACCTGGCCCGACCCTTGACCTGCCGCCTCCCTTCCGGACCGATGGGTGGACCCCGGAGCTGAGGCATAGACACACATATCGCTCCTTGGTCACACGATTCGTTTAATTGTCTGGCGCAGGTCTAGGCTAGTAAGATGGTCAGAAGAAACACGAGACGTTACTAAGTGAGCCTCCCTTACGCTTCGTTCGCGTTCTTCCATCGTAGGGCCACAGTTCTAGTTTGTGCCCGATCTTGGGTAGCCCTTCTACTTTCAGCTAAGACAAATCACTGCGGGAACACAGGACTTCCGTGTCCCAGTGGCCTAGTGCGGATAACATTTTACTGGAAGCGGATCTGGGTCACGGGCCAATGCTCTAGACAGCATCGCTCGGTATGGTAGCCAGTAACGATGACGAACTAGGAGATATCCTCTTTAATACTTCTGCCACCCTTTCAACACATCTTTTGGGAACCCATCGGTGCATCCCTTGCGCTGTTTTGATTGCAAAACCGCTACAAGCAGATCTTAACGTAGAGGTCGTCTGGGGGCGTCTCCCTCCC', 'ACAGGGGACTACATGGCGATCTGTGGGGTTATACGGGTGTACTACACACTTTCTTTACGATAGGAGGTACCCTTATGGCACCTGGGTCTAATGTTCCTACACATGAACTAGCAACCCACCCTACCCCCCGATTGTTCTAAAGTATGACACACGAATATAGCGTAGCTCACCGTGTGATTTAACAGGCCGCATCACGGTGTTCGTTCGAAACTACCACCTCCCTCATGCATCACAGCACCACCGGAAGGTATTAAAAGCCGCCCTTAGCCACTTCTGAAGCAATCTATTTTGACCGAACGTAAATATTGCGATCTTAGCTCGACACGCAGACGGGTATTACCGCCTCGGCGTAGTAATGCGGGTGGGAGGAACGGCGTCGAGACTCCTCCCTATGACACGCACAAGGCCCCACAGAATCTCTGAGATATGGGTGCCTGACAGTCCGTGAAAACGCCGCGGTAGTTACTCTCAAAGTTATTAACATGCAGATGTAAACCGTGAGCCGGTTTCCTCGGGCCACGTACGGGGGCTGGGTACACTCTACTCAATGATATCCCAATGATCCTGTCCGCACTGTCGCAGTCAGCGATGCTCTAGGTCCAAAGAAATGATGGGAATCATGTTAGCCGACTTTGTTAACTGCCGGTATATCCTGGGAGGCTAACCGACTGCATAACATTAGGATCGCTAGGCGAAGCTCACCCAATGACGACAAGGTAGAACAGAAGTAAAGAGGAAACCAACCGATCTTACTGCTTACCACCGTTCCCTTGGCTTTCCATAACTGAACTTATTAACAGAGTCTATTTCATGTACTTCGCGCAATCCTGGTATTTGTCCGGGTGGTGCCGGCT')
    #>>> ddg.backtrack(ddg.longest_path())
    >>> dag = SequenceAlignment('AACCTTGG', 'ACACTGTGA')
    >>> dag.backtrack(dag.longest_path())
    'AACTTG'
    """

    def __init__(self, v, w):
        self.v = v
        self.w = w
        self.n = len(v)
        self.m = len(w)
        dynamic.Dag.__init__(self)

    def prepare_sides(self):
        pass

    def inbounds(self, i, j):
        yield self.vertices[i][j + 1] + 0
        yield self.vertices[i + 1][j] + 0
        if self.v[i] == self.w[j]:
            yield self.vertices[i][j] + 1

    @staticmethod
    def resolve_vertex(dag, vrtx, backtr):
        if vrtx.real:
            backtr.append(dag.w[vrtx.coords[0]-1])

class BlosumAlign(dynamic.Dag):

    """
    # Blosum 62 scoring matrix
    >>> g = BlosumAlign('PLEASANTLY', 'MEANLY')
    >>> p = g.longest_path()
    >>> print p
    8
    >>> print g.backtrack(p)
    -MEA--N-LY
    """

    blosum_62 = {}
    blosum_dat = open('./BLOSUM62.txt', 'r')
    line_splitter = re.compile("[\w-]+")
    aa_keys = line_splitter.findall(blosum_dat.readline().strip())
    for line in blosum_dat.readlines():
        dat = line_splitter.findall(line.strip())
        blosum_62[dat[0]] = dict(zip(aa_keys, map(lambda s: int(s), dat[1:])))
    blosum_dat.close()

    def __init__(self, v, w, sigma=-5):
        self.v = v
        self.w = w
        self.n = len(v)
        self.m = len(w)
        self.sigma = sigma
        dynamic.Dag.__init__(self, vrtx_real_dtrm=lambda x: not x == sigma)

    def prepare_sides(self):
        for i in range(self.n):
            self.vertices[i + 1][0].update(self.vertices[i][0] + self.sigma)
        for i in range(self.m):
            self.vertices[0][i + 1].update(self.vertices[0][i] + self.sigma)

    def inbounds(self, i, j):
        yield self.vertices[i][j + 1] + self.sigma
        yield self.vertices[i + 1][j] + self.sigma
        yield self.vertices[i][j] + self.blosum_62[self.v[i]][self.w[j]]

    @staticmethod
    def resolve_vertex(dag, vrtx, backtr):
        if vrtx.real:
            backtr.append(dag.w[vrtx.coords[0]-1])
        else:
            backtr.append('-')

