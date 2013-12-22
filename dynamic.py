import numpy
import re


class Dag:

    def __init__(self, vrtx_real_dtrm=False):
        if not self.n:
            self.n = 0
            self.m = 0
        self.vertices = [[DagVertex(0, coords=[i, j], real_dtrm=vrtx_real_dtrm) for i in range(self.m + 1)] for j in range(self.n + 1)]

    def __repr__(self):
        return numpy.matrix(self.vertices).__repr__()

    def longest_path(self):
        """
        #>>> strs = ["3 0 0 3 1 0 1 2 1 0 2 0", "3 3 4 1 0 1 2 1 3 0 2 1", "0 0 3 0 0 2 3 1 4 2 4 0", "0 3 4 0 2 4 0 2 4 4 2 1", "2 4 0 2 1 4 1 1 0 0 0 2", "-", "1 0 0 3 1 2 4 3 0 1 4", "1 2 3 2 4 1 4 2 3 4 0", "0 3 2 3 1 3 4 4 1 0 3", "0 4 1 1 4 3 0 0 2 4 4", "2 3 4 2 2 4 3 3 3 3 4", "1 1 4 0 4 4 3 4 1 3 3"]
        #>>> Manhattan.from_strings(strs).longest_path()
        47
        >>> m = Manhattan.from_strings(["1 0 2 4 3", "4 6 5 2 1", "4 4 5 2 1", "5 6 8 5 3", "-", "3 2 4 0", "3 2 4 2", "0 7 3 3", "3 3 0 2", "1 3 2 2"])
        >>> m.longest_path()
        34
        >>> Manhattan.from_strings(["1 0 2 4 3", "4 6 5 2 1", "4 4 5 2 1", "5 6 8 5 3", "-", "3 2 4 0", "3 2 4 2", "0 7 3 3", "3 3 0 2", "1 3 2 2", "-", "5 0 2 1", "8 4 3 0", "10 8 9 5", "5 6 4 7"]).longest_path()
        35
        """
        self.prepare_sides()

        for i in range(self.n):
            for j in range(self.m):
                self.vertices[i + 1][j + 1].update(max(self.inbounds(i, j)))
        return self.vertices[self.n][self.m]

    def backtrack(self, vertex, reslv_fnc=False):
        """
        Backtrack from
        """
        backtr = []
        if not reslv_fnc:
            reslv_fnc = self.resolve_vertex
        while vertex.ptr:
            reslv_fnc(self, vertex, backtr)
            vertex = vertex.ptr
        return ''.join(backtr[::-1])

    @staticmethod
    def resolve_vertex(dag, vrtx, backtr):
        raise StandardError('Not Implemented')

    def prepare_sides(self):
        raise StandardError('Not Implemented')

    def inbounds(self, i, j):
        raise StandardError('Not Implemented')


class DagVertex:

    real_dtrm = lambda slf, x: x > 0

    def __init__(self, weight, coords=[], ptr=0, real=False, direction=0, name = '', real_dtrm=False):
        """
        @type direction: int # if normal 1, if horizontal -1, if vertical -2
        @type coords: list
        @type weight: int
        """
        self.weight = weight
        self.inbounds = {}
        self.outbounds = {}
        self.coords = coords
        self.ptr = ptr
        self.real = real
        self.direction = direction
        self.name = name
        if real_dtrm:
            self.real_dtrm = real_dtrm

    def update(self, vrtx):
        self.weight = vrtx.weight
        self.ptr = vrtx.ptr
        self.real = vrtx.real
        self.direction = vrtx.direction

    def calculate_weight(self):
        if self.indegree():
            longest_in = max(self.inbounds.values(), key=lambda edg: edg.path_weight())
            self.weight += longest_in.path_weight()
            self.ptr = longest_in.src

    def indegree(self):
        return len(self.inbounds)

    def __repr__(self):
        if self.name:
            return self.name
        else:
            return self.weight.__repr__()

    def __add__(self, other, direction=0):
        if isinstance(other, int):
            return self.__class__(self.weight + other, ptr=self, direction=direction, real=(self.real_dtrm(other)))
        elif isinstance(other, self.__class__):
            return self.__class__(self.weight + other.weight, ptr=self, direction=direction)
        else:
            raise AttributeError('Not supported type')

    def __cmp__(self, other):
        return cmp(self.weight, other.weight)

class DagEdge:

    def __init__(self, src, dst, weight=1):
        self.src = src
        self.dst=dst
        self.weight = int(weight)
        src.outbounds[dst.name] = self
        dst.inbounds[src.name] = self

    def __repr__(self):
        return "%s-(%i)->%s" % (self.src.name, self.weight, self.dst.name)

    def path_weight(self):
        return self.src.weight + self.weight

class ArbitryDag(Dag):

    def __init__(self):
        self.vertices = {}
        self.ordering = []
        self.edges = []

    def __repr__(self):
        return ';'.join(map(lambda x: x.name, self.vertices.values()))

    @classmethod
    def from_lines(cls, lines=[]):
        """
        >>> ArbitryDag.from_lines(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3']).edges
        [0-(7)->1, 0-(4)->2, 2-(2)->3, 1-(1)->4, 3-(3)->4]
        """
        parser = re.compile('(\d*)->(\d*):(\d*)')
        graph = cls()
        for line in lines:
            graph.add_edge(*parser.search(line.strip()).groups())
        return graph

    def add_edge(self, src, dst, weight):
        src = self.get_vrtx(src)
        dst = self.get_vrtx(dst)
        edge = DagEdge(src, dst, weight)
        self.edges.append(edge)

    def get_vrtx(self, name):
        if self.vertices.has_key(name):
            return self.vertices[name]
        else:
            vrtx = DagVertex(0, name=name)
            self.vertices[name] = vrtx
            return vrtx

    def starting_nodes(self):
        return [vtx for vtx in self.vertices.values() if not vtx.indegree()]

    def topological_ordering(self):
        """
        >>> ag = ArbitryDag.from_lines(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3'])
        >>> ag.topological_ordering()
        >>> ag.ordering
        ['0', '2', '3', '1', '4']
        """
        self.ordering = []
        ordering = self.starting_nodes()
        indegrees = {}
        while len(ordering):
            src = ordering.pop()
            self.ordering.append(src.name)
            for edge in src.outbounds.values():
                if indegrees.has_key(edge.dst.name):
                    indegrees[edge.dst.name] += 1
                else:
                    indegrees[edge.dst.name] = 1
                if edge.dst.indegree() == indegrees[edge.dst.name]:
                    ordering.append(edge.dst)

    def longest_path(self, dst=''):
        """
        >>> ArbitryDag.from_lines(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3']).longest_path().weight
        Graph not odered. Ordering
        9
        """
        if not self.ordering:
            print 'Graph not odered. Ordering'
            self.topological_ordering()
        if not dst:
            dst = self.ordering[-1]
        for node in self.ordering:
            node = self.get_vrtx(node)
            node.calculate_weight()
        return self.get_vrtx(dst)

    @staticmethod
    def backtrack(vertex, src=''):
        """
        >>> d = ArbitryDag.from_lines(['0->1:7', '0->2:4', '2->3:2', '1->4:1', '3->4:3']).longest_path()
        Graph not odered. Ordering
        >>> print ArbitryDag.backtrack(d)
        0->2->3->4
        """
        back_track = vertex.name
        while vertex.ptr:
            vertex = vertex.ptr
            back_track = ("%s->" % vertex.name) + back_track
            if vertex.name == src:
                break
        return back_track


class Manhattan(Dag):
    def __init__(self, *directions):
        self.down = list(directions[0])
        self.right = list(directions[1])
        self.directions = directions[2:]
        self.n, self.m = len(self.down), len(self.right[0])
        Dag.__init__(self)

    @classmethod
    def from_strings(cls, iterable):
        """
        From lines of down and right separated with '-'
        >>> m = Manhattan.from_strings(["1 0 2 4 3", "4 6 5 2 1", "4 4 5 2 1", "5 6 8 5 3", "-", "3 2 4 0", "3 2 4 2", "0 7 3 3", "3 3 0 2", "1 3 2 2"])
        >>> m.down, m.right, m.n, m.m
        ([[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1], [5, 6, 8, 5, 3]], [[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2], [1, 3, 2, 2]], 4, 4)
        """
        directions = [[]]
        current = 0
        for line in iterable:
            if line == "-":
                current += 1
                if len(directions) <= current:
                    directions.append([])
                continue
            directions[current].append(map(lambda i: int(i), line.split(' ')))
        return cls(*directions)


    def prepare_sides(self):
        for i in range(self.n):
            self.vertices[i + 1][0].update(self.vertices[i][0] + self.down[i][0])
        for i in range(self.m):
            self.vertices[0][i + 1].update(self.vertices[0][i] + self.right[0][i])



    def inbounds(self, i, j):
        yield self.vertices[i][j + 1] + self.down[i][j + 1]
        yield self.vertices[i + 1][j] + self.right[i + 1][j]
        for table in self.directions:
            yield self.vertices[i][j] + table[i][j]


class Coins:
    def __init__(self, coins=[50, 25, 20, 10, 5, 1]):
        self.coins = [c for c in coins]
        self.min_coins = {0: 0}

    def change(self, money):
        """
        >>> Coins().change(40)
        2
        >>> Coins([20,15,9,8,5,3,1]).change(17015)
        851
        """
        for m in range(1, money + 1):
            if self.min_coins.has_key(m):
                continue
            local_min = m * 1000
            for coin in self.coins:
                if m >= coin:
                    if self.min_coins[m - coin] + 1 < local_min:
                        local_min = self.min_coins[m - coin] + 1
            self.min_coins[m] = local_min
        return self.min_coins[money]

    def min_values_ordered(self):
        """
        >>> c = Coins([6,5,1])
        >>> c.change(12)
        2
        >>> print ' '.join(map(lambda i: str(i), c.min_values_ordered()))
        0 1 2 3 4 1 1 2 3 4 2 2 2
        >>> c.change(23)
        4
        >>> print ' '.join(map(lambda i: str(i), c.min_values_ordered()))
        0 1 2 3 4 1 1 2 3 4 2 2 2 3 4 3 3 3 3 4 4 4 4 4
        """
        for i in range(len(self.min_coins)):
            yield self.min_coins[i]
