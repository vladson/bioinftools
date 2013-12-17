import numpy


class Manhattan:

    def __init__(self, down=numpy.matrix([]), right=numpy.matrix([]), n=0, m=0):
        self.down = numpy.matrix(down)
        self.right = numpy.matrix(right)

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
