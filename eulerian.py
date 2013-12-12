import random, dna, assembler

class Graph:

    def __init__(self):
        self.nodes = {}
        self.edges = set()

    @classmethod
    def from_repr_lines(cls, lines):
        """
        #>>> Graph.from_repr_lines(["0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6"])
        1 -> 0
        0 -> 3
        3 -> 2
        2 -> 1,6
        5 -> 4
        4 -> 2
        7 -> 9
        6 -> 5,8
        9 -> 6
        8 -> 7
        """
        graph = cls()
        for line in lines:
            for start, end in Graph.parse_line(line):
                graph.add_edge(start, end)
        return graph

    @classmethod
    def from_kmers(cls, iterable):
        """
        #>>> Graph.from_kmers(['00', '01', '10', '11'])
        1 -> 1,0
        0 -> 1,0
        """
        graph = cls()
        for kmer in iter(iterable):
            graph.add_edge(*assembler.DeBruijn.split(kmer))
        return graph


    def eulerian_cycle(self):
        """
        Not representative output
        #>>> g = Graph.from_repr_lines(["0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6"])
        #>>> cycle = g.eulerian_cycle()
        #>>> cycle.reformat_start(g.get_node('6'))
        #>>> print cycle
        6->8->7->9->6->5->4->2->1->0->3->2->6
        >>> c = Graph.from_repr_lines(["CTT -> TTA", "ACC -> CCA", "TAC -> ACC", "GGC -> GCT", "GCT -> CTT", "TTA -> TAC"]).eulerian_cycle()
        >>> c
        GGC->GCT->CTT->TTA->TAC->ACC->CCA
        >>> print c.assemble(dna.Dna.nodes_overlap_assembler)
        GGCTTACCA
        """
        if not self.balanced():
            self.balance()
        cycle = Cycle(self)
        cycle.perform()
        while len(cycle.outways):
            cycle.reformat_start()
            cycle.perform()
        if hasattr(self, 'fake'):
            cycle.reformat_start(self.fake.end)
            cycle.traverse.pop()
        return cycle

    def add_edge(self, start, end):
        if isinstance(start, Node) and isinstance(end, Node):
            start_node, end_node = start, end
        else:
            start_node = self.get_node(start)
            end_node = self.get_node(end)
        edge = Edge(start_node, end_node)
        self.edges.add(edge)
        return edge

    def output(self):
        for _, node in self.nodes.iteritems():
            yield "%s -> %s" % (node.__repr__(), ','.join(map(lambda o: o.end.__repr__(), node.outbounds)))

    def get_node(self, key):
        if self.nodes.has_key(key):
            return self.nodes[key]
        else:
            node = Node(key)
            self.nodes[key] = node
            return node

    def get_random_node(self):
        return random.choice(self.nodes.values())

    def __repr__(self):
        return '\n'.join(self.output())

    def balanced(self):
        """
        >>> Graph.from_repr_lines(["0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6"]).balanced()
        True
        >>> Graph.from_repr_lines(["0 -> 2", "1 -> 3", "2 -> 1", "3 -> 0,4", "6 -> 3,7", "7 -> 8", "8 -> 9", "9 -> 6"]).balanced()
        False
        """
        return reduce(lambda was, node: was and node.balanced(), self.nodes.values(), True)

    def balance(self):
        """
        >>> graph = Graph.from_repr_lines(["0 -> 2", "1 -> 3", "2 -> 1", "3 -> 0,4", "6 -> 3,7", "7 -> 8", "8 -> 9", "9 -> 6"])
        >>> graph.balance()
        >>> graph.balanced()
        True
        """
        unbalanced = filter(lambda node: not node.balanced(), self.nodes.values())
        if len(unbalanced) > 2:
            raise StandardError("Heavily unbalanced")
        elif len(unbalanced) == 0:
            print "Graph balanced"
        else:
            self.fake = self.balance_pair(*unbalanced)


    @staticmethod
    def parse_line(line):
        """
        >>> list(Graph.parse_line('2 -> 1,6'))
        [('2', '1'), ('2', '6')]
        """
        start, ends = line.strip().split(' -> ')
        for end in ends.split(','):
            yield start, end

    def balance_pair(self, first, second):
        if first.disbalance() + second.disbalance() == 0:
            if first.disbalance() < 0:
                start, end = second.key, first.key
            else:
                start, end = first.key, second.key
        else:
            print first.disbalance() - second.disbalance()
            raise StandardError("Nodes are not equaly unbalanced")
        return self.add_edge(start, end)

class Node:

    def __init__(self, key):
        self.key = key
        self.inbounds = {}
        self.outbounds = {}

    def __repr__(self):
        """
        >>> str(Node(7))
        '7'
        """
        return str(self.key)

    def add_outbound(self, node):
        if not self.outbounds.has_key(Node):
            self.outbounds[node] = node
        else:
            raise StandardError("Node %s already in outbounds" % str(node))

    def add_inbound(self, node):
        if not self.inbounds.has_key(Node):
            self.inbounds[node] = node
        else:
            raise StandardError("Node %s already in inbounds" % str(node))

    def disbalance(self):
        return len(self.inbounds) - len(self.outbounds)

    def balanced(self):
        return self.disbalance() == 0

class Edge:

    def __init__(self, start, end, representation=''):
        self.start = start
        self.end = end
        if representation:
            self.representation = representation
        else:
            self.representation = "%s-%s" % (str(start), str(end))
        start.add_outbound(self)
        end.add_inbound(self)

    def __repr__(self):
        """
        >>> str(Edge(Node(3), Node(5)))
        '3-5'
        """
        return self.representation

class Cycle:

    def __init__(self, graph, start=False):
        self.graph = graph
        if not start:
            if hasattr(graph, 'fake'):
                start = graph.fake.end
            else:
                start = graph.get_random_node()
        self.traverse = [start]
        self.visited = {str(start): []}
        self.current = start
        self.outways = set()

    def assemble(self, funct):
        return funct(self.traverse)

    def perform(self):
        try:
            while self.has_ways():
                self.move()
        except KeyboardInterrupt:
            print self.traverse
            print self.visited
            print self.outways

    def move(self):
        new = self.get_way()
        self.visited[str(self.current)].append(new)
        self.traverse.append(new)
        if not self.visited.has_key(str(new)):
            self.visited[str(new)] = []
        if self.has_ways():
            self.outways.add(self.current)
        else:
            if self.current in self.outways:
                self.outways.remove(self.current)
        self.current = new

    def has_ways(self):
        return len(self.visited[str(self.current)]) < len(self.current.outbounds)

    def get_way(self):
        way = random.choice(self.current.outbounds.values()).end
        while way in self.visited[str(self.current)]:
            way = random.choice(self.current.outbounds.values()).end
        return way

    def reformat_start(self, outway=False):
        if outway:
            if not isinstance(outway, Node):
                outway = self.graph.get_node(outway)
        else:
            outway = self.outways.pop()
            self.outways.add(outway)
        pos = self.traverse.index(outway)
        self.traverse = self.traverse[pos:] + self.traverse[1:pos + 1]
        self.current = outway

    def __repr__(self):
        return '->'.join(map(lambda n: n.__repr__(), self.traverse))
