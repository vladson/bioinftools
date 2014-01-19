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

    def __init__(self, text, children=None, starts=None, seq_ids=None):
        if children == None:
            children = []
        if starts == None:
            starts = []
        if seq_ids == None:
            self.seq_ids = set()
        elif isinstance(seq_ids, int) or isinstance(seq_ids, str):
            self.seq_ids = set(seq_ids)
        else:
            self.seq_ids = seq_ids
        self.label = text
        self.children = children
        self.starts = starts

    def __repr__(self):
        return "%s (%s)" % (self.label, ', '.join(map(lambda x: str(x), self.starts)))

    def append(self, text, start, seq_id=None):
        if text == '':
            # print 'none'
            self.starts.append(start)
            if seq_id:
                self.seq_ids.add(seq_id)
        else:
            try:
                child = next(edge for edge in self.children if edge.label.startswith(text[0]))
                if text.startswith(child.label):
                    # print "child append %s" % text
                    child.append(text[len(child.label):], start, seq_id)
                else:
                    # print "split %s" % text
                    child.split(text, start, seq_id)
            except StopIteration:
                # print "append %s" % text
                self.children.append(SuffixEdge(text, starts=[start], seq_ids=seq_id))

    def split(self, text, start, seq_id=None):
        split_offset = 0
        for i in range(min(len(self.label), len(text))):
            split_offset = i
            if not self.label[i] == text[i]:
                break
        # print split_offset, self.label, text, self.children
        edge_1 = SuffixEdge(self.label[split_offset:], self.children[:], self.starts, self.seq_ids)
        edge_2 = SuffixEdge(text[split_offset:], starts=[start], seq_ids=seq_id)
        self.label = self.label[:split_offset]
        if seq_id:
            self.seq_ids.add(seq_id)
        self.children = [edge_1, edge_2]

    def longest_repeat(self, parent):
        if self.children:
            child_seq = [child.longest_repeat(parent + self.label) for child in self.children]
            return max(child_seq, key=lambda x:len(x))
        else:
            return parent

    def traverse(self, buffer=None):
        if self.label:
            if isinstance(buffer, list):
                buffer.append(self.label)
            else:
                print self.label
        for edge in self.children:
            edge.traverse(buffer)

    def bfs_paths(self, condition=lambda n, c, p: c.label[-1] == '$', out=lambda n, c, p: p + n.label + c.label):
        """
        @param comparer: function, lambda taking node child and path
        @return: generator of path, can
        >>> list(SuffixTree('TCGGTAGATTGCGCCCACTC$', 'A').root.bfs_paths())
        ['$', 'AGATTGCGCCCACTC$', 'ATTGCGCCCACTC$', 'ACTC$', 'GGTAGATTGCGCCCACTC$', 'GTAGATTGCGCCCACTC$', 'GATTGCGCCCACTC$', 'GCGCCCACTC$', 'GCCCACTC$', 'CACTC$', 'CTC$', 'C$', 'CCCACTC$', 'CCACTC$', 'CGGTAGATTGCGCCCACTC$', 'CGCCCACTC$', 'TAGATTGCGCCCACTC$', 'TTGCGCCCACTC$', 'TGCGCCCACTC$', 'TCGGTAGATTGCGCCCACTC$', 'TC$']
        """
        queue = [(self, '')]
        while len(queue):
            (node, path) = queue.pop()
            for child in node.children:
                if condition(node, child, path):
                    yield out(node, child, path)
                else:
                    queue.append((child, path + node.label))


    def shared_suffix(self, sequence):
        """
        >>> SuffixTree('TCGGTAGATTGCGCCCACTC$').root.shared_suffix('AGAA')
        'AGA'
        """
        shared_part = ''
        if sequence:
            if len(self.label) == 0 or sequence.startswith(self.label):
                try:
                    child = next(child for child in self.children if child.label[0] == sequence[len(self.label)])
                    shared_part = self.label + child.shared_suffix(sequence[len(self.label):])
                except StopIteration:
                    shared_part = self.label
                except IndexError:
                    shared_part = self.label[:len(sequence)]
            else:
                for i in xrange(min(len(sequence), len(self.label))):
                    if not sequence[i] == self.label[i]:
                        shared_part = sequence[:i]
                        break
        # print 'returning %s' % shared_part
        return shared_part


class SuffixTree:

    """
    >>> SuffixTree('ATAAATG$')
    [A (0), T (1), G$ (6), $ (7)]
    >>> SuffixTree('ATAAATG$').traverse()
    A
    T
    AAATG$
    G$
    A
    ATG$
    TG$
    T
    AAATG$
    G$
    G$
    $
    >>> SuffixTree('ATAAATG$').traverse([])
    ['A', 'T', 'AAATG$', 'G$', 'A', 'ATG$', 'TG$', 'T', 'AAATG$', 'G$', 'G$', '$']
    >>> SuffixTree('ATATCGTTTTATCGTT').longest_repeat()
    'TATCGT'
    """

    def __init__(self, text, id=None):
        self.root = SuffixEdge('')
        self.append_named_sequence(text, id)

    def __repr__(self):
        return self.root.children.__repr__()

    def append_named_sequence(self, sequence, id=None):
        for i in xrange(len(sequence)):
            self.root.append(sequence[i:], start=i, seq_id=id)
        return self

    def to_suffix_array(self):
        """
        @return: suffix array
        >>> SuffixTree('panamabananas$').to_suffix_array()
        [13, 5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12]
        """
        out = lambda n, c, p: (p + n.label + c.label, c.starts[0])
        return map(lambda (s, i): i, sorted(self.root.bfs_paths(out=out), key=lambda x: x[0]))

    def longest_repeat(self):
        return self.root.longest_repeat('')

    def traverse(self, buffer=None):
        self.root.traverse(buffer)
        return buffer

    @classmethod
    def longest_shared_between(cls, seq_1, seq_2):
        """
        >>> SuffixTree.longest_shared_between('TCGGTAGATTGCGCCCACTC', 'AGGGGCTCGCAGTGTAAGAA')
        'AGA'
        """
        tree = cls(seq_2 + '$', 'B').append_named_sequence(seq_1 + '$', 'A')
        condition = lambda n, c, p: 'B' not in c.seq_ids
        out = lambda n, c, p: p + n.label
        unshared = tree.root.bfs_paths(condition, out)
        shortest = next(unshared)
        common = [shortest]
        for possible in unshared:
            if len(possible) > len(shortest) and not possible[-1] == '$':
                shortest = possible
                common = [possible]
            elif len(possible) == len(shortest) and not possible[-1] == '$':
                common.append(possible)
        return min(common)

    @classmethod
    def shortest_unshared_between(cls, seq_1, seq_2):
        """
        >>> SuffixTree.shortest_unshared_between('CCAAGCTGCTAGAGG', 'CATGCTGGGCTGGCT')
        'AA'
        """
        tree = cls(seq_2 + '$', 'B').append_named_sequence(seq_1 + '$', 'A')
        condition = lambda n, c, p: 'B' not in c.seq_ids
        out = lambda n, c, p: p + n.label + c.label[0]
        unshared = tree.root.bfs_paths(condition, out)
        shortest = next(unshared)
        common = [shortest]
        for possible in unshared:
            if len(possible) < len(shortest) and not possible[-1] == '$':
                shortest = possible
                common = [possible]
            elif len(possible) == len(shortest) and not possible[-1] == '$':
                common.append(possible)
        return min(common)