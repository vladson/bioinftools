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
    >>> print g.backtrack(p, BlosumAlign.resolve_vertex_1)
    PLEASANTLY
    >>> print g.backtrack(p, BlosumAlign.resolve_vertex_2)
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
            self.vertices[i + 1][0].update(self.vertices[i][0].__add__(self.sigma, -2))
        for i in range(self.m):
            self.vertices[0][i + 1].update(self.vertices[0][i].__add__(self.sigma, -1))

    def inbounds(self, i, j):
        yield self.vertices[i][j + 1].__add__(self.sigma, -2)
        yield self.vertices[i + 1][j].__add__(self.sigma, -1)
        yield self.vertices[i][j].__add__(self.blosum_62[self.v[i]][self.w[j]], 1)


    @staticmethod
    def resolve_vertex_1(dag, vrtx, backtr):
        if vrtx.direction == 1 or vrtx.direction == -2:
            backtr.append(dag.v[vrtx.coords[1]-1])
        elif vrtx.direction == -1:
            backtr.append('-')


    @staticmethod
    def resolve_vertex_2(dag, vrtx, backtr):
        if vrtx.direction == 1 or vrtx.direction == -1:
            backtr.append(dag.w[vrtx.coords[0]-1])
        elif vrtx.direction == -2:
            backtr.append('-')

