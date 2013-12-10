class OverlapGraph:

    def __init__(self, graph={}):
        self.graph = graph

    def __repr__(self):
        return '\n'.join(self.output())

    @classmethod
    def from_kmers(cls, kmeriterator):
        """
        >>> OverlapGraph.from_kmers(["ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT"])
        GGCAT -> GCATG
        AGGCA -> GGCAT
        CATGC -> ATGCG
        GCATG -> CATGC
        """
        graph = cls()
        for kmer in iter(kmeriterator):
            graph.add_kmer(kmer)
        return graph

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

