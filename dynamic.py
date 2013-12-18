import numpy


class Manhattan:

    def __init__(self, down=[[]], right=[[]], *direction_costs):
        self.down = list(down)
        self.right = list(right)
        self.directions = direction_costs
        self.n, self.m = len(down), len(right[0])
        self.vertices = [ [ 0 for _ in range(self.m+1)] for _ in range(self.n+1) ]

    def __repr__(self):
        return numpy.matrix(self.vertices).__repr__()

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
        for i in range(self.n):
            self.vertices[i+1][0] = self.vertices[i][0] + self.down[i][0]
        for i in range(self.m):
            self.vertices[0][i+1] = self.vertices[0][i] + self.right[0][i]
        for i in range(self.n):
            for j in range(self.m):
                self.vertices[i+1][j+1] = max(self.inbounds(i, j))
        return self.vertices[self.n][self.m]

    def inbounds(self, i, j):
        yield self.vertices[i][j+1] + self.down[i][j+1]
        yield self.vertices[i+1][j] + self.right[i+1][j]
        for table in self.directions:
            yield self.vertices[i][j] + table[i][j]

class Coins:

    def __init__(self, coins = [50,25,20,10,5,1]):
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
            local_min = m*1000
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
