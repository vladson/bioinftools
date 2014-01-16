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



class SuffixTrie:

    def __init__(self, string="$"):
        self.root = [{}]
        for index in xrange(len(string.strip()) - 1):
            try:
                self.add_string(string[index:], index)
            except KeyboardInterrupt:
                print "Currently we are at %i" % index
                raise

    def __repr__(self):
        return self.root.__repr__()

    def add_string(self, string, index):
        c_index = 0
        for letter in string:
            if self.root[c_index].has_key(letter):
                c_index = self.root[c_index][letter]
            else:
                root_len = len(self.root)
                self.root[c_index][letter] = root_len
                self.root.append({})
                c_index = root_len


    def match(self, iterable):
        """
        @param iterable: str
        @return: bool
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

class SuffixEdge:

    def __init__(self, text, children=[], starts=[]):
        self.label = text
        self.children = children
        self.starts = starts

    def __repr__(self):
        return "%s (%s)" % (self.label, ', '.join(map(lambda x: str(x), self.starts)))

    def append(self, text, start):
        if text:
            edges = [edge for edge in self.children if edge.label.startswith(text[0])]
        if text == '' or text == self.label:
            print 'none'
            self.starts.append(start)
        elif len(edges) > 0:
            child = edges[0]
            if text.startswith(child.label):
                print 'append'
                child.append(text[len(child.label):], start)
            else:
                print 'split'
                child.split(text, start)
        else:
            print 'append'
            self.children.append(SuffixEdge(text[len(self.label):], starts=[start]))

    def split(self, text, start):
        split_offset = 0
        for i in range(min(len(self.label), len(text))):
            split_offset = i
            if not self.label[i] == text[i]:
                break
        print split_offset, self.label, text
        edge_1 = SuffixEdge(self.label[split_offset:], self.children, self.starts)
        edge_2 = SuffixEdge(text[split_offset:], starts=[start])
        self.label = self.label[:split_offset]
        self.children = [edge_1, edge_2]

    def longest_repeat(self, parent=''):
        if self.children:
            child_seq = [child.longest_repeat(parent + self.label) for child in self.children]
            return max(child_seq, key=lambda x:len(x))
        else:
            return parent

    def traverse(self):
        print self.children.count()
        for edge in self.children:
            edge.traverse()


class SuffixTree:

    """
    >>> SuffixTree('ATAAATG$')
    """

    def __init__(self, text):
        self.root = SuffixEdge('')
        for i in xrange(len(text)):
            self.root.append(text[i:], start=i)

    def __repr__(self):
        return self.root.children.__repr__()

    def longest_repeat(self):
        return self.root.longest_repeat()


