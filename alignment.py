import dynamic

class SequenceAlignment(dynamic.Dag):

    """
    >>> dag = SequenceAlignment('AACCTTGG', 'ACACTGTGA')
    >>> dag.backtrack(dag.longest_path())
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
            backtr.append(dag.v[vrtx.coords[0]])



