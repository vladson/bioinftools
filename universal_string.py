import dna

class UniversalBinaryString(dna.Dna):
    """
    >>> list(UniversalBinaryString.all_kmer_generator(2))
    ['00', '01', '10', '11']
    """

    alfabet = '0', '1'