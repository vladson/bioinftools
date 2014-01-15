import dna, rna, protein, assembler, eulerian, universal_string, dynamic, alignment, rearrangements
import read_mapping


class Stepic:

    test = True

    #
    # Read Mapping
    #

    @staticmethod
    def trie_matching(path):
        data = open(path)
        reference_genome = data.readline().strip()
        print 'Constructing'
        trie = read_mapping.Trie(*data.readlines())
        print 'Matching'
        print ' '.join(map(lambda x: str(x), trie.long_string_matches_indices(reference_genome)))

    @staticmethod
    def trie_construction(path):
        data = open(path)
        print 'Constructing'
        trie = read_mapping.Trie(*data.readlines())
        data.close()
        output = open('output/trie.txt', 'w')
        for (i, j, l) in trie.triplets():
            output.write("%i %i %s\n" % (i,j,l))
        output.close()

    #
    # Genetic rearrangements
    #

    @staticmethod
    def shared_kmers(path, test=False):
        data = open(path)
        if test:
            data.readline()
        k = int(data.readline().strip())
        genome1 = data.readline().strip()
        genome2 = data.readline().strip()
        constructor = rearrangements.SyntenyConstructor(genome1, genome2)
        if test:
            data.readline()
            test_data = map(lambda s: s.strip(), data.readlines())
        data.close()
        print constructor
        otp = open('output/shared_kmers.txt', 'w')
        print "Beginning shared kmers search with k = %i" % k
        for item in constructor.shared_kmer_indices(k):
            if test and str(item) in test_data:
                test_data.remove(str(item))
            otp.write(str(item) + '\n')
        otp.close()
        print "Done. Written."
        if test:
            print 'Test data mismatch:'
            print test_data

    @staticmethod
    def double_break_distance(path, test=False):
        data = open(path)
        if test:
            data.readline()
        genome1 = data.readline().strip()
        genome2 = data.readline().strip()
        graph = rearrangements.BreakpointGraph.from_genome(genome1).add_genomes(genome2)
        print "Graph of %i nodes. Beginning cycles calculation" % graph.block_num()
        graph.calculate_cycles()
        print "Calculated! Cycles: %i, blocks %i, Double Break distance %i" % (graph.cycles_num(), graph.block_num(), graph.db_distance())

    @staticmethod
    def greedy_sorting(path, test= False):
        data = open(path)
        if test:
            data.readline()
        permutation = rearrangements.Permutation(data.readline())
        print "Permutation is %i long" % len(permutation.perm)
        otp = open('output/greedy_sorting.txt', 'w')
        print "Beginning"
        for item in permutation.greedy_sorting():
            otp.write(rearrangements.Permutation.out_perm(item) + '\n')
        otp.close()
        print "Done. Written."

    @staticmethod
    def breakpoints_count(path, test=False):
        data = open(path)
        if test:
            data.readline()
        permutation = rearrangements.Permutation(data.readline())
        print "Permutation is %i long" % len(permutation.perm)
        print "Beginning breakpoints"
        breakpoints = permutation.breakpoints()
        print "Breakpoints"
        print breakpoints
        print "Breakpoints count is %i" % len(breakpoints)

    #
    #   Sequence alignment
    #

    @staticmethod
    def edit_distance(path, test=False):
        data = open(path)
        if test:
            data.readline()
        src = data.readline().strip()
        dst = data.readline().strip()
        if test:
            data.readline()
            check = data.readline()
        data.close()
        graph = alignment.ScoredAlign(src, dst, -5)
        print 'Beginning calculations of edit distance'
        results, node = graph.edit_distance()
        print results
        if test:
            otp = open('output/edit_distance.txt', 'w')
            print "should be %s" % check
            otp.write(str(results)+"\n")
            otp.write(check)
            otp.write(src + "\n")
            otp.write(dst + "\n")
            otp.write(graph.backtrack(node, graph.__class__.resolve_vertex_levenstein))
            otp.close()
        #return graph, node

    @staticmethod
    def local_align(path, test=False):
        data = open(path)
        src = data.readline().strip()
        dst = data.readline().strip()
        graph = alignment.ScoredAlign(src, dst, -5, scorer=alignment.ScoredAlign.pam250)
        print 'Beginning calculations'
        node, align_1st, align_2nd = graph.local_alignment()
        print 'Score'
        print node
        print 'Backtrack'
        print align_1st
        print align_2nd
        output = open('./output/local_align.txt', 'w')
        output.write(str(node.weight) + "\n")
        if test:
            output.write(src + "\n")
        output.write(align_1st + "\n")
        if test:
            output.write("\n")
            output.write(dst + "\n")
        output.write(align_2nd + "\n")
        print 'results written to ./output/local_align.txt'
        #return graph, path_node
        exit()

    @staticmethod
    def global_align(path):
        data = open(path)
        src = data.readline().strip()
        dst = data.readline().strip()
        graph = alignment.ScoredAlign(src, dst, -5)
        print 'Beginning calculations'
        node, align_1st, align_2nd = graph.global_alignment()
        print 'Score'
        print node
        print 'Backtrack'
        print align_1st
        print align_2nd
        output = open('./output/global_align.txt', 'w')
        output.write(str(node.weight) + "\n")
        output.write(align_1st + "\n")
        output.write(align_2nd + "\n")
        print 'results written to ./output/global_align.txt'
        #return graph, path_node
        exit()

    @staticmethod
    def dag_longest(path):
        data = open(path)
        src = data.readline().strip()
        dst = data.readline().strip()
        graph = dynamic.ArbitryDag.from_lines(data.readlines())
        print graph
        print "Starting search in the graph with source: %s, sink: %s" % (src, dst)
        node = graph.longest_path(dst)
        src_node = graph.get_vrtx(src)
        print "Node weight: %i, src.weight: %i, diff_path_weight: %i" % (node.weight, src_node.weight, node.weight - src_node.weight)
        print "Length and backtrack"
        print node.weight - src_node.weight
        print graph.backtrack(node, src)
        return graph


    #
    #   Assemble genome
    #

    @staticmethod
    def contigs(path, test=False):
        data = open(path)
        if test:
            data.readline()
        if test:
            lines = []
            while True:
                line = data.readline()
                if line.strip() == 'Output:':
                    break
                lines.append(line)
            output = data.readlines()
        else:
            lines = data.readlines()
        data.close()
        graph = eulerian.Graph.from_kmers(lines)
        if not test:
            results = open('./output/contigs.txt', 'w')
            for contig in graph.contigs():
                results.write(contig+'\n')
            results.close()
            print "Results saved"
        else:
            contigs = list(graph.contigs())
            print "Graph edges %i, results len %i, answer len %i" % (len(graph.edges), len(contigs), len(output))
            print 'Results'
            print '\n'.join(contigs)
            print 'Output'
            print output


    @staticmethod
    def read_pairs_reconstructor(path, test=False, split_func=dna.Dna.kdmer_splitter):
        data = open(path)
        if test:
            data.readline()
        d = int(data.readline().strip())
        if test:
            lines = []
            while True:
                line = data.readline()
                if line.strip() == 'Output:':
                    break
                lines.append(line)
            output = data.readline()
        else:
            lines = data.readlines()
        data.close()
        graph = eulerian.Graph.from_kmers(lines, split_func)
        cycle = graph.eulerian_cycle()
        print "Graph edges %i, results len %i" % (len(graph.edges), len(cycle.traverse))
        while len(graph.edges) > len(cycle.traverse):
            cycle = graph.eulerian_cycle()
            print "Graph edges %i, results len %i" % (len(graph.edges), len(cycle.traverse))
        if not test:
            results = open('./output/readpairs_reconstruction.txt', 'w')
            results.write(cycle.assemble(dna.Dna.get_kdmer_nodes_overlap_assembler(d)))
            results.close()
            print "Results saved"
        else:
            print "Graph edges %i, results len %i, answer len %i" % (len(graph.edges), len(cycle.traverse), len(output.split('->')))
            print 'Results'
            print cycle.assemble(dna.Dna.get_kdmer_nodes_overlap_assembler(d))
            print 'Output'
            print output

    @staticmethod
    def universal_k_binary_string(k, write=False):
        """
        #>>> Stepic.universal_k_binary_string(3)
        00011101
        #>>> Stepic.universal_k_binary_string(4)
        0000110010111101
        """
        graph = eulerian.Graph.from_kmers(universal_string.UniversalBinaryString.all_kmer_generator(k))
        cycle = graph.eulerian_cycle()
        string = cycle.assemble(universal_string.UniversalBinaryString.nodes_cyrcular_assembler)
        if write:
            results = open('./output/universal_k_binary_string.txt', 'w')
            results.write(string)
            results.close()
        print string

    @staticmethod
    def eulerian_cycle(path, assembler_func=False, test = False):
        data = open(path)
        if test:
            data.readline()
            lines = []
            while True:
                line = data.readline()
                if line.strip() == 'Output':
                    break
                lines.append(line)
            graph = eulerian.Graph.from_repr_lines(lines)
            output = data.readline()
        else:
            graph = eulerian.Graph.from_repr_lines(data.readlines())
        data.close()
        cycle = graph.eulerian_cycle()
        print "Graph edges %i, results len %i" % (len(graph.edges), len(cycle.traverse))
        while len(graph.edges) > len(cycle.traverse):
            cycle = graph.eulerian_cycle()
            print "Graph edges %i, results len %i" % (len(graph.edges), len(cycle.traverse))
        if not test:
            results = open('./output/eulerian_cycle.txt', 'w')
            if assembler_func:
                results.write(cycle.assemble(assembler_func))
            else:
                results.write(cycle.__repr__())
            results.close()
        else:
            print "Graph edges %i, results len %i, answer len %i" % (len(graph.edges), len(cycle.traverse), len(output.split('->')))
            print 'Results'
            if assembler_func:
                print cycle.assemble(assembler_func)
            else:
                print cycle
            print 'Output'
            print output


    @staticmethod
    def debrujin_string(path):
        data = open(path)
        k = int(data.readline().strip())
        seq = data.readline().strip()
        print "Assembling DeBrujin from %s with k %i" % (seq, k)
        graph = assembler.DeBruijn.from_dna(seq, k)
        results = open('./output/debrujin_string.txt', 'w')
        results.write(graph.__repr__())
        results.close()
        print graph

    @staticmethod
    def debrujin_kmers(path):
        data = open(path)
        kmers = []
        for line in data.readlines():
            kmers.append(line.strip())
        data.close()
        print "Assembling from %i kmers" % len(kmers)
        graph = assembler.DeBruijn.from_kmers(kmers)
        results = open('./output/debrujin_kmers.txt', 'w')
        results.write(graph.__repr__())
        results.close()
        print graph

    @staticmethod
    def overlap(path):
        data = open(path)
        kmers = []
        for line in data.readlines():
            kmers.append(line.strip())
        data.close()
        print "Assembling from %i kmers" % len(kmers)
        graph = assembler.OverlapGraph.from_kmers(kmers)
        results = open('./output/overlap.txt', 'w')
        results.write(graph.__repr__())
        results.close()
        print graph

    @staticmethod
    def composition(path, test=False):
        data = open(path)
        if test:
            data.readline()
        k = int(data.readline().strip())
        genome = dna.Dna(data.readline().strip())
        if test:
            data.readline()
            print "Answers"
            for line in data.readlines():
                print line
        data.close()
        print "Results"
        results = open('./output/composition.txt', 'w')
        for kmer in genome.kmer_composition(k):
            print kmer
            results.write(kmer+'\n')
        results.close()

    #   End of Assemble Genome
    #   Motiff finding
    #

    @staticmethod
    def profile_most_probable_kmer(path):
        data = open(path)
        seq = dna.Dna(data.readline().strip())
        k = int(data.readline().strip())
        profile = dna.Profile.from_lines_legend(data.readlines())
        print "For dna %s and profile %s" % (seq, profile)
        print seq.most_probable_motiff(profile)

    @staticmethod
    def greedy_motiff_search(path, pseudocounts=True, test=False):
        data = open(path)
        if test:
            data.readline()
        k, t = map(lambda x: int(x), data.readline().strip().split())
        motiffs = []
        for i in range(t):
            motiffs.append(data.readline().strip())
        if test:
            data.readline()
            answers = []
            for line in data.readlines():
                answers.append(line.strip())
        data.close()
        print "For k %i and t %i with %i motiffs" % (k, t, len(motiffs))
        if pseudocounts:
            print 'Pseudocounts ON'
        results = dna.Dna.greedy_motiff_search(motiffs, k, t, pseudocounts)
        print '\n'.join(results)
        if test:
            print 'ANSWERS'
            print '\n'.join(answers)
            print 'DIFF'
            for i in range(len(answers)):
                if results[i] != answers[i]:
                    print answers[i]

    @staticmethod
    def randomised_motiff_search(path, treshold=10, test=False):
        data = open(path)
        if test:
            data.readline()
        k, t = map(lambda x: int(x), data.readline().strip().split())
        motiffs = []
        for i in range(t):
            motiffs.append(data.readline().strip())
        if test:
            data.readline()
            answers = []
            for line in data.readlines():
                answers.append(line.strip())
            treshold = dna.Dna.motiff_scorer(answers)
            print "Treshold set to answers score and is %i" % treshold
        data.close()
        print "For k %i and treshold %i with %i motiffs" % (k, treshold, len(motiffs))
        results = dna.Dna.iterative_randomized_motiff_search(motiffs, k, t, treshold, True)
        print '\n'.join(results)
        if test:
            print 'ANSWERS'
            print '\n'.join(answers)
            print 'DIFF'
            for i in range(len(answers)):
                res = sorted(results)[i]
                ans = sorted(answers)[i]
                if res != ans:
                    print ans

    @staticmethod
    def gibbs_sampler(path, treshold=10, test=False):
        data = open(path)
        if test:
            data.readline()
        k, t, n = map(lambda x: int(x), data.readline().strip().split())
        motiffs = []
        for i in range(t):
            motiffs.append(data.readline().strip())
        if test:
            data.readline()
            answers = []
            for line in data.readlines():
                answers.append(line.strip())
            treshold = dna.Dna.motiff_scorer(answers)
            print "Treshold set to answers score and is %i" % treshold
        data.close()
        print "For k %i and treshold %i with %i motiffs and %i iterations per call" % (k, treshold, len(motiffs), n)
        results = dna.Dna.iterative_gibbs_sampler(motiffs, k, t, n, treshold, True)
        print '\n'.join(results)
        if test:
            print 'ANSWERS'
            print '\n'.join(answers)
            print 'DIFF'
            for i in range(len(answers)):
                res = sorted(results)[i]
                ans = sorted(answers)[i]
                if res != ans:
                    print ans