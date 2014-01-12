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
            perm = map(lambda s: int(s), representation.strip()[1:-1].split(' '))
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

