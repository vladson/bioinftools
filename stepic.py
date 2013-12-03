import dna, rna, protein
class Stepic:

    test = True

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

