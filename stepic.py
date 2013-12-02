import dna, rna, protein
class Stepic:

    test = True

    @staticmethod
    def profile_most_probable_kmer(path):
        data = open(path)
        seq = dna.Dna(data.readline().strip().split())
        k = int(data.readline().strip().split())
        profile = dna.Profile.from_lines_legend(data.readlines())
        print "For dna %s and profile %s" % (seq, profile)
        print seq.most_probable_motiff(profile)

