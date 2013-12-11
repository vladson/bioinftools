import random

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
        graph = Graph()
        for line in lines:
            for start, end in Graph.parse_line(line):
                graph.add_edge(start, end)
        return graph

    def eulerian_cycle(self):
        """
        >>> g = Graph.from_repr_lines(["0 -> 3", "1 -> 0", "2 -> 1,6", "3 -> 2", "4 -> 2", "5 -> 4", "6 -> 5,8", "7 -> 9", "8 -> 7", "9 -> 6"])
        >>> cycle = g.eulerian_cycle()
        >>> cycle.reformat_start(g.get_node('6'))
        >>> print cycle
        6->8->7->9->6->5->4->2->1->0->3->2->6
        """
        cycle = Cycle(self)
        cycle.perform()
        while len(cycle.outways):
            cycle.reformat_start()
            cycle.perform()
        return cycle

    def add_edge(self, start, end):
        start_node = self.get_node(start)
        end_node = self.get_node(end)
        self.edges.add(Edge(start_node, end_node))

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
        """
        return reduce(lambda was, node: was and node.balanced(), self.nodes.values(), True)

    @staticmethod
    def parse_line(line):
        """
        >>> list(Graph.parse_line('2 -> 1,6'))
        [('2', '1'), ('2', '6')]
        """
        start, ends = line.strip().split(' -> ')
        for end in ends.split(','):
            yield start, end


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
            raise "Node %s already in outbounds" % str(node)

    def add_inbound(self, node):
        if not self.inbounds.has_key(Node):
            self.inbounds[node] = node
        else:
            raise "Node %s already in inbounds" % str(node)

    def balanced(self):
        return len(self.inbounds) == len(self.outbounds)


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
            start = graph.get_random_node()
        self.traverse = [start]
        self.visited = {str(start): []}
        self.current = start
        self.outways = []


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
        if self.has_ways() and not self.current in self.outways:
            self.outways.append(self.current)
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
        if not outway:
            outway = random.choice(self.outways)
        pos = self.traverse.index(outway)
        self.traverse = self.traverse[pos:] + self.traverse[1:pos + 1]
        self.current = outway

    def __repr__(self):
        return '->'.join(map(lambda n: n.__repr__(), self.traverse))
