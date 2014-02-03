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

class SuffixArray:

    def __init__(self, repr):
        if isinstance(repr, list):
            self.indices = repr
        else:
            self.indices = []


    def __repr__(self):
        return self.indices.__repr__()

    @classmethod
    def from_sequence(cls, sequence):
        """
        >>> SuffixArray.from_sequence('AACGATAGCGGTAGA$')
        [15, 14, 0, 1, 12, 6, 4, 2, 8, 13, 3, 7, 9, 10, 11, 5]
        """
        proto_indices = list(range(len(sequence)))
        def seq_sort(i, j):
            return cmp(sequence[i:], sequence[j:])
        return cls(sorted(proto_indices, cmp=seq_sort ))

    def h_out(self):
        return ', '.join(map(lambda x: str(x), self.indices))

class BWT:

    """
    >>> BWT('GCGTGCCTGGTCA$')
    'ACTGGCT$TGCGGC'
    >>> BWT(last_col='enwbpeoseu$llt')
    'enwbpeoseu$llt'
    >>> BWT('panamabananas$').partial_suffix_array
    {1: 5, 11: 10, 12: 0}
    """

    def __init__(self, sequence=None, last_col=None, sa_slice=5):
        def bwt_cmp(i, j):
            return cmp(sequence[i:], sequence[j:])
        if sequence:
            self.tl = len(sequence)
            indices = sorted(list(range(self.tl)), cmp=bwt_cmp)
            self.last_col = ''.join([sequence[indx-1] for indx in indices])
            if sa_slice:
                self.partial_suffix_array = {ndx: vl for (ndx, vl) in enumerate(indices) if not vl % sa_slice}
        elif last_col:
            self.tl = len(last_col)
            self.last_col = last_col
            self.partial_suffix_array = None
        else:
            raise StandardError('No valid sequence provided')
        self.first_col = None
        self.first_alfa_index = None
        self.last_alfa_index = None

    def first_column(self):
        """
        @return: str
        >>> BWT(last_col='ard$rcaaaabb').first_column()
        '$aaaaabbcdrr'
        """
        if not self.first_col:
            self.first_col = ''.join(sorted(self.last_col))
        return self.first_col

    def __repr__(self):
        return self.last_col.__repr__()

    def reconstruct(self):
        """
        @return: str
        >>> BWT('panamabananas$').reconstruct()
        'panamabananas$'
        >>> BWT(last_col='TTCCTAACG$A').reconstruct()
        'TACATCACGT$'
        """
        self.first_column()
        first_alfa_index, last_alfa_index = self.setup_mapping()

        previos_index = last_alfa_index['$'][0]
        reconstruction = self.first_col[previos_index]
        for i in xrange(1, self.tl):
            last = reconstruction[-1]
            next_index = first_alfa_index[last].index(previos_index)
            previos_index = last_alfa_index[last][next_index]
            reconstruction += self.first_col[previos_index]
        return reconstruction

    def match_num(self, pattern):
        """
        @param pattern: str
        @return: int
        >>> b = BWT(last_col='TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC')
        >>> b.match_num('CCT')
        2
        >>> b.match_num('CAG')
        0
        >>> b.match_num('ATC')
        1
        """
        self.setup_mapping()
        pattern = list(pattern)
        top, bottom = 0, self.tl
        while top <= bottom:
            if pattern:
                current = pattern.pop()
                occurencies = filter(lambda x: top <= x <= bottom, self.last_alfa_index[current])
                if occurencies:
                    top = self.first_alfa_index[current][self.last_alfa_index[current].index(occurencies[0])]
                    bottom = self.first_alfa_index[current][self.last_alfa_index[current].index(occurencies[-1])]
                else:
                    return 0
            else:
                return bottom - top + 1

    def better_matching(self, pattern, offsets=False, slice_len=5):
        """
        >>> b = BWT(last_col='TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC')
        >>> b.better_matching('CCT')
        2
        >>> b = BWT(last_col='GGCGCCGC$TAGTCACACACGCCGTA')
        >>> b.better_matching('ACC')
        1
        >>> b.better_matching('CCG')
        2
        >>> b.better_matching('CAG')
        0
        >>> b.better_matching('GAT')
        0
        >>> b.better_matching('ACC', slice_len=5, offsets=False)
        1
        >>> b = BWT('AATCGGGTTCAATCGGGGT')
        >>> list(b.better_matching('ATCG', True))
        [11, 1]
        """
        first_occurencies, counts_at = self.setup_better_match(slice_len)
        pattern = list(pattern)
        top, bottom = 0, self.tl
        while top <= bottom:
            if pattern:
                current = pattern.pop()
                if counts_at(bottom)[current] - counts_at(top)[current] > 0:
                    top = first_occurencies[current] + counts_at(top)[current]
                    bottom = first_occurencies[current] + counts_at(bottom + 1)[current] - 1
                else:
                    if offsets:
                        return []
                    else:
                        return 0
            else:
                if offsets:
                    return self.offsets(top, bottom)
                else:
                    return bottom - top + 1


    def offsets(self, top, bottom):
        """
        >>> b = BWT('panamabananas$')
        >>> s = b.setup_better_match(5)
        >>> list(b.offsets(5,5))
        [9]
        >>> list(b.offsets(3,5))
        [1, 7, 9]
        """
        if not self.partial_suffix_array:
            return
        for i in xrange(top, bottom + 1):
            addition = 0
            curr = i
            while not self.partial_suffix_array.has_key(curr):
                curr = self.step_back(curr)
                addition += 1
            yield self.partial_suffix_array[curr] + addition

    def setup_mapping(self):
        self.first_column()
        if not self.first_alfa_index:
            if not hasattr(self, 'alfa_list'):
                self.alfa_list = list(set(self.last_col))
            self.first_alfa_index = dict(map(lambda x: (x, []), self.alfa_list))
            self.last_alfa_index = dict(map(lambda x: (x, []), self.alfa_list))
            for index, alfa in enumerate(self.first_col):
                self.first_alfa_index[alfa].append(index)
            for index, alfa in enumerate(self.last_col):
                self.last_alfa_index[alfa].append(index)
        return self.first_alfa_index, self.last_alfa_index

    def setup_better_match(self, slice_len=1):
        """
        >>> BWT('panamabananas$').setup_better_match()[0]
        {'a': 1, 'b': 7, '$': 0, 'm': 8, 'n': 9, 'p': 12, 's': 13}
        >>> b = BWT('panamabananas$')
        >>> first_o,counts = b.setup_better_match(5)
        >>> counts(0)
        {'a': 0, 'p': 0, 's': 0, 'b': 0, '$': 0, 'm': 0, 'n': 0}
        >>> counts(13)
        {'a': 5, 'p': 1, 's': 1, 'b': 1, '$': 1, 'm': 1, 'n': 3}
        >>> counts(14)
        {'a': 6, 'p': 1, 's': 1, 'b': 1, '$': 1, 'm': 1, 'n': 3}
        >>> first_o['$']
        0
        >>> first_o['s']
        13
        """
        if not hasattr(self, 'first_occurencies'):
            self.slice_len = slice_len
            if not hasattr(self, 'alfa_list'):
                self.alfa_list = sorted(list(set(self.last_col)))
            counts = []
            current = {i : 0 for i in self.alfa_list}
            for index, char in enumerate(self.last_col):
                if not index % slice_len:
                    counts.append(current.copy())
                    current[char] += 1
                else:
                    current[char] += 1
            counts.append(current.copy())
            self.counts = counts
            first_occurencies = {char: 0 for char in self.alfa_list}
            current_offset = 0
            for char in sorted(self.alfa_list):
                first_occurencies[char] = current_offset
                current_offset += counts[-1][char]
            self.first_occurencies = first_occurencies

            def counts_at(num):
                start_slice_pos = num / slice_len
                start = slice_len * start_slice_pos
                end = min(self.tl, start + num % slice_len)
                current_counts = counts[start_slice_pos].copy()
                for i in xrange(start, end):
                    current_counts[self.last_col[i]] += 1
                return current_counts

            self.counts_at = counts_at

        if slice_len != self.slice_len:
            print 'WARNING: better match set up with slice_len = %i' % self.slice_len
        return self.first_occurencies, self.counts_at

    def step_back(self, curr):
        """
        >>> b = BWT('panamabananas$')
        >>> s = b.setup_better_match(5)
        >>> b.step_back(4)
        7
        """
        char = self.last_col[curr]
        return self.first_occurencies[char] + self.counts_at(curr)[char]






