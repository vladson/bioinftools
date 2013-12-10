import assembler
import dna, rna, protein
class Stepic:

    test = True

    #
    #   Assemble genome
    #

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