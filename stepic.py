import dna, rna, protein, assembler, eulerian, universal_string

class Stepic:

    test = True

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