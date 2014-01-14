import re
import dna


class Permutation:

    def __init__(self, representation, dst = False):
        """
        @param repr: list or str
        @param dst: list (opt)
        @return: Permutaion
        >>> Permutation([-3, 4, 1, 5, -2])
        (-3 +4 +1 +5 -2)
        >>> Permutation(' (-3 +4 +1 +5 -2)  ')
        (-3 +4 +1 +5 -2)
        >>> Permutation(' (-3 +4 +1 +5 -2)  ').out_destination()
        '(+1 +2 +3 +4 +5)'
        """
        if isinstance(representation, list):
            perm = representation
        else:
            perm =Permutation.parse_repr(representation)[0]
        self.perm = perm
        self.initital_perm = perm
        if dst:
            if isinstance(dst, list):
                self.dst = dst
            else:
                self.dst = map(lambda s: int(s), dst.strip()[1:-1].split(' '))
        else:
            self.dst = self.default_destination()
        self.path = []

    def greedy_sorting(self):
        """
        >>> len(Permutation([-3, 4, 1, 5, -2]).greedy_sorting())
        7
        >>> map(lambda x: Permutation.out_perm(x), Permutation([-3, 4, 1, 5, -2]).greedy_sorting())
        ['(-1 -4 +3 +5 -2)', '(+1 -4 +3 +5 -2)', '(+1 +2 -5 -3 +4)', '(+1 +2 +3 +5 +4)', '(+1 +2 +3 -4 -5)', '(+1 +2 +3 +4 -5)', '(+1 +2 +3 +4 +5)']
        """
        for i in xrange(len(self.perm)):
            if self.perm[i] == self.dst[i]:
                continue
            else:
                self.greedy_sorting_iter(i)
        return self.path

    def greedy_sorting_iter(self, i):
        try:
            j = self.perm.index(self.dst[i])
        except ValueError:
            j = self.perm.index(-self.dst[i])
        self.perm = self.reversal(i, j)
        self.path.append(self.perm[:]) #record
        if -self.perm[i] == self.dst[i]:
            self.greedy_sorting_iter(i)

    def reversal(self, i, j):
        """
        @param i: int
        @param j: int
        @return: list
        >>> Permutation([-3, 4, 1, 5, -2]).reversal(0,0)
        [3, 4, 1, 5, -2]
        >>> Permutation([-3, 4, 1, 5, -2]).reversal(0,2)
        [-1, -4, 3, 5, -2]
        """
        return self.perm[:i] + map(lambda s: s*-1, self.perm[i:j+1][::-1]) + self.perm[j+1:]

    def breakpoints(self):
        """
        @return: int
        >>> len(Permutation('(+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14)').breakpoints())
        8
        """
        breakpoints = []
        perm = [0] + self.perm + [self.dst[-1]+1]
        last_j = 0
        for i in xrange(1, len(perm)):
            if not self.adjacent(perm[last_j:i+1]):
                breakpoints.append(i-1)
                last_j = i
        return breakpoints

    @staticmethod
    def adjacent(lst):
        """
        @param lst: list
        @return: bool
        >>> Permutation.adjacent([3,4,5])
        True
        >>> Permutation.adjacent([-8, -7, -6])
        True
        >>> Permutation.adjacent([10, 9])
        False
        >>> Permutation.adjacent([5, -12])
        False
        """
        for i in range(1, len(lst)):
            if not lst[i-1] + 1 == lst[i]:
                return False
        return True

    def default_destination(self):
        """
        @return: list
        >>> Permutation([-3, 4, 1, 5, -2]).default_destination()
        [1, 2, 3, 4, 5]
        """
        return sorted(map(lambda x: abs(x), self.perm))

    def __repr__(self):
        return self.__class__.out_perm(self.perm)

    def out_destination(self):
        return self.__class__.out_perm(self.dst)

    @staticmethod
    def out_perm(perm):
        """
        @param perm: []
        @return: str
        >>> Permutation.out_perm([1, 2, 3, 4, 5])
        '(+1 +2 +3 +4 +5)'
        >>> Permutation.out_perm([-3, 4, 1, 5, -2])
        '(-3 +4 +1 +5 -2)'
        """
        repr = []
        for i in perm:
            if i > 0:
                repr.append("+%i" % i)
            if i < 0:
                repr.append("%i" % i)
        return '(' + " ".join(repr) + ')'

    @staticmethod
    def parse_repr(repr):
        """
        >>> Permutation.parse_repr('(+1 -3 -6 -5)(+2 -4)')
        [[1, -3, -6, -5], [2, -4]]
        >>> Permutation.parse_repr(' (-3 +4 +1 +5 -2)  ')[0]
        [-3, 4, 1, 5, -2]
        """
        return map(lambda literal: map(lambda item: int(item), literal.split()) ,re.findall('\(([^()]*)\)+', repr))

class CircularGenome:

    def __init__(self, genome):
        self.genome = genome
        self.genome.append(genome[0])

    def pairs(self):
        for i in xrange(len(self.genome) - 1):
            yield self.genome[i:i+2]

class BreakpointGraph:

    """
    >>> BreakpointGraph.from_genome('(+1 +2 +3 +4 +5 +6)')
    (+1 +2 +3 +4 +5 +6)
    >>> BreakpointGraph.from_genome('(+1 -3 -6 -5)(+2 -4)')
    (+1 +2 +3 +4 +5 +6)
    >>> BreakpointGraph.from_genome('(+1 +2 +3 +4 +5 +6)').add_genomes('(+1 -3 -6 -5)(+2 -4)').calculate_cycles().cycles_num()
    3
    >>> BreakpointGraph.from_genome('(+1 +2 +3 +4 +5 +6)').add_genomes('(+1 -3 -6 -5)(+2 -4)').calculate_cycles().db_distance()
    3
    """

    def __init__(self):
        self.cycles = []
        self.nodes = {}
        self.blocks = set()

    @classmethod
    def from_genome(cls, genomes):
        graph = cls()
        if isinstance(genomes, str):
            genomes = Permutation.parse_repr(genomes)
        graph.add_genomes(genomes)
        return graph

    def add_genomes(self, genomes):
        if isinstance(genomes, str):
            genomes = Permutation.parse_repr(genomes)
        for genome in genomes:
            self.add_genome(genome)
        return self

    def add_genome(self, genome):
        if not isinstance(genome, CircularGenome):
            genome = CircularGenome(genome)
        for pair in genome.pairs():
            self.add_pair(*pair)

    def calculate_cycles(self):
        tracker = list(self.nodes.keys())
        while len(tracker):
            self.cycles.append(BreakpointCycle(self.get_node(tracker.pop())).traverse(tracker))
        return self

    def add_pair(self, node_a, node_b):
        a = self.get_node(node_a)
        b = self.get_node(node_b * -1)
        a.add_edge(b)
        b.add_edge(a)

    def get_node(self, node):
        """
        @return BreakpointNode
        """
        if not self.nodes.has_key(node):
            self.nodes[node] = BreakpointNode(node)
            self.blocks.add(abs(node))
        return self.nodes[node]

    def block_num(self):
        """
        >>> BreakpointGraph.from_genome('(+1 +2 +3 +4 +5 +6)').block_num()
        6
        """
        return len(self.blocks)

    def cycles_num(self):
        return len(self.cycles)

    def db_distance(self):
        return self.block_num() - self.cycles_num()

    def __repr__(self):
        return Permutation.out_perm(self.blocks)

class BreakpointNode:

    def __init__(self, node):
        self.node = node
        self.edges = []

    def __repr__(self):
        return "%i (%s)" % (self.node, ', '.join(map(lambda n: str(n.node), self.edges)))

    def add_edge(self, dst):
        self.edges.append(dst)

    def has_edges(self):
        return len(self.edges)

    def pop_edge(self):
        node = self.edges.pop()
        node.edges.remove(self)
        return node

class BreakpointCycle:

    def __init__(self, node):
        self.cycle = [node]

    def __repr__(self):
        return self.cycle.__repr__()

    def traverse(self, tracker):
        while self.cycle[-1].has_edges():
            node = self.cycle[-1].pop_edge()
            self.cycle.append(node)
            if node.node in tracker:
                tracker.remove(node.node)
        return self

class SyntenyConstructor:

    """
    >>> SyntenyConstructor('ACTG', 'ATGTA')
    SyntenyConstructor for dna_A of 4 bp and dna_B of 5 bp
    """

    def __init__(self, dna_a, dna_b):
        if not isinstance(dna_a, dna.Dna):
            dna_a = dna.Dna(dna_a)
        if not isinstance(dna_b, dna.Dna):
            dna_b = dna.Dna(dna_b)
        self.dna_a = dna_a
        self.dna_b = dna_b

    def __repr__(self):
        return "SyntenyConstructor for dna_A of %i bp and dna_B of %i bp" % (len(self.dna_a), len(self.dna_b))

    def shared_kmer_indices(self, k):
        """
        >>> list(SyntenyConstructor('AAACTCATC', 'TTTCAAATC').shared_kmer_indices(3))
        [(0, 4), (0, 0), (4, 2), (6, 6)]
        """
        b_kmer_tree = {}
        for index, kmer in enumerate(self.dna_b.kmer_generator(k)):
            if b_kmer_tree.has_key(kmer):
                b_kmer_tree[kmer].append(index)
            else:
                b_kmer_tree[kmer] = [index]

        for index, kmer in enumerate(self.dna_a.kmer_generator(k)):
            if b_kmer_tree.has_key(kmer):
                for jindex in b_kmer_tree[kmer]:
                    yield index, jindex
            # if b_kmer_tree.has_key(kmer[::-1]):
            #     yield index, b_kmer_tree[kmer[::-1]]
            if b_kmer_tree.has_key(dna.Dna.complementary_fragment(kmer)):
                for jindex in b_kmer_tree[dna.Dna.complementary_fragment(kmer)]:
                    yield index, jindex
            # if b_kmer_tree.has_key(dna.Dna.complementary_fragment(kmer)[::-1]):
            #     yield index, b_kmer_tree[dna.Dna.complementary_fragment(kmer)[::-1]]
