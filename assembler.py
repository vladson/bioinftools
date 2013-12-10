import dna

class Graph:

    def __init__(self):
        self.graph = {}

    def __repr__(self):
        return '\n'.join(self.output())

    @classmethod
    def from_dna(cls, seq, k):
        if not isinstance(seq, dna.Dna):
            seq = dna.Dna(seq)
        return cls.from_kmers(seq.kmer_generator(k))

    @classmethod
    def from_kmers(cls, kmeriterator):
        graph = cls()
        for kmer in iter(kmeriterator):
            graph.add_kmer(kmer)
        return graph

    def add_kmer(self, kmer):
        raise 'NotImplemented'

    def output(self):
        raise 'NotImplemented'


class DeBruijn(Graph):
    """
    >>> DeBruijn.from_dna('AAGATTCTCTAC', 4)
    AAG -> AGA
    TCT -> CTC,CTA
    GAT -> ATT
    AGA -> GAT
    ATT -> TTC
    CTA -> TAC
    CTC -> TCT
    TTC -> TCT
    >>> DeBruijn.from_kmers(["GAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"])
    GAG -> AGG
    AGG -> GGG
    GGG -> GGG,GGA
    CAG -> AGG
    GGA -> GAG
    """

    def add_kmer(self, kmer):
        start, end = DeBruijn.split(kmer)
        if start in self.graph:
            self.graph[start].append(end)
        else:
            self.graph[start] = [end]

    def output(self):
        for node, neighbours in self.graph.iteritems():
            yield "%s -> %s" % (node, ','.join(neighbours))

    @staticmethod
    def split(kmer):
        return kmer[:-1], kmer[1:]

    @staticmethod
    def edge(start, end):
        """
        >>> DeBruijn.edge('AAG', 'AGA')
        'AAGA'
        """
        assert len(start) == len(end)
        return start + end[-1]


class OverlapGraph(Graph):

    """
    >>> OverlapGraph.from_kmers(["ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT"])
    GGCAT -> GCATG
    AGGCA -> GGCAT
    CATGC -> ATGCG
    GCATG -> CATGC
    """

    def add_kmer(self, kmer):
        if not self.graph.has_key(kmer):
            checker, edges  = OverlapGraph.get_kmer_checker(kmer), []
            for node, (node_checker, node_edges) in self.graph.iteritems():
                if node_checker(kmer):
                    node_edges.append(kmer)
                if checker(node):
                    edges.append(node)
            self.graph[kmer] = (checker, edges)

    def output(self):

        for node, (_, edges) in self.graph.iteritems():
            for edge in edges:
                yield "%s -> %s" % (node, edge)

    @staticmethod
    def get_kmer_checker(pattern):
        """
        >>> c = OverlapGraph.get_kmer_checker('AGGCA')
        >>> c('GGCAT')
        True
        >>> c('GGCAG')
        True
        >>> c('GTCAG')
        False
        """
        def kmer_checker(kmer):
            assert len(pattern) == len(kmer)
            return pattern[1:] == kmer[:-1]
        return kmer_checker

