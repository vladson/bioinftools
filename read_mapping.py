import dna

class Trie:

    """
    >>> Trie('GGTA', 'CG', 'GGC').h_out()
    1 6 C
    1 2 G
    2 3 G
    3 8 C
    3 4 T
    4 5 A
    6 7 G
    """

    def __init__(self, *strings):
        self.root = [{}, {}]
        for string in strings:
            self.add_string(string)

    def __repr__(self):
        return self.root.__repr__()

    def add_string(self, string=""):
        c_index = 1
        for letter in string.strip():
            if self.root[c_index].has_key(letter):
                c_index = self.root[c_index][letter]
            else:
                root_len = len(self.root)
                self.root[c_index][letter] = root_len
                self.root.append({})
                c_index = root_len

    def long_string_matches_indices(self, long_string):
        """
        Fixme: use intertools.tee to improve computation speed
        @param long_string:
        @return: generator of int
        >>> list(Trie('ATCG', 'GGGT').long_string_matches_indices('AATCGGGTTCAATCGGGGT'))
        [1, 4, 11, 15]
        """
        for index in xrange(len(long_string) - 1):
            if self.match(long_string[index:]):
                yield index

    def match(self, iterable):
        """
        @param iterable: str
        @return: bool
        >>> trie = Trie('ATCG', 'GGGT')
        >>> trie.match('AATCGGGTTCAATCGGGGT')
        False
        >>> trie.match('ATCGGGTTCAATCGGGGT')
        True
        >>> trie.match('GGGTTCAATCGGGGT')
        True
        """
        index = 1
        for letter in iterable:
            if self.root[index].has_key(letter):
                index = self.root[index][letter]
                if not len(self.root[index]):
                    return True
            else:
                return False
        return False

    def triplets(self):
        for index, node in enumerate(self.root):
            if len(node):
                for key, value in node.iteritems():
                    yield index, value, key

    def h_out(self):
        for (i, j, l) in self.triplets():
            print i, j, l