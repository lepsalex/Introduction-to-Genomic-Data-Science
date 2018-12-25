from functools import reduce
from math import log2

# Implement GreedyMotifSearch() below along with any subroutines that you need.
# Your function should return a list of strings.
def GreedyMotifSearch(Dna, k, t):
    return GreedyMotifSearchWithScore(Dna, k, t)[0]


def GreedyMotifSearchWithScore(Dna, k, t):
    bestMotifs = [x[:k] for x in Dna]
    bestMotifsScore = computeEntropyMotif(getFrequencyMatrix(bestMotifs))
    firstLineKmers = [Dna[0][x:x + k] for x in range(0, len(Dna[0]) - k + 1)]
    linesToAnalyze = Dna[1:t]
    return reduce(lambda x, y: GmsReducer(x, y, linesToAnalyze, k), firstLineKmers, (bestMotifs, bestMotifsScore))


def GmsReducer(acc, curr, linesToAnalyze, k):
    # Pass to the next reducer the lines[1:t] and the starting motif
    currMotifs = reduce(lambda x, y: MotifReducer(
        x, y, k), linesToAnalyze, [curr])
    currScore = computeEntropyMotif(getFrequencyMatrix(currMotifs))
    return (currMotifs, currScore) if currScore < acc[1] else acc


def MotifReducer(acc, curr, k):
    # Create a profile from all motifs in our list
    # first pass will have the passed in starting motif
    profile = getFrequencyMatrix(acc)
    return acc + [ProfileMostProbableKmer(curr, k, profile)]


# Reduce to kmer with highest score tuple, return the kmer string
def ProfileMostProbableKmer(text, k, profile):
    return reduce(lambda acc, curr: ProfileReduce(acc, curr, profile), [
        text[x:x + k] for x in range(0, len(text) - k + 1)], (text[0:k], 0))[0]

# Compare pm score and return higher as tuple (kmer, score)
def ProfileReduce(acc, curr, profile):
    pm = ComputeProfileMatrix(curr, profile)
    return (curr, pm) if pm > acc[1] else acc

# ComputeProfileMatrix(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753
def ComputeProfileMatrix(pattern, profile):
    return reduce(lambda x, y: x * y, (profile[nuc][idx] for idx, nuc in enumerate(pattern)))


def getFrequencyMatrix(input):
    counts = generateEmptyCountMatrix(len(input[0]))

    for line in input:
        for idx, nuc in enumerate(line):
            counts[nuc][idx] += 1 / len(input)

    return counts


def generateEmptyCountMatrix(length):
    return {
        'A': [0] * length,
        'C': [0] * length,
        'G': [0] * length,
        'T': [0] * length
    }


def computeEntropyMotif(freqMatrix):
    entropy = [0] * len(list(freqMatrix.values())[0])

    for frequencies in freqMatrix.values():
        for idx, freq in enumerate(frequencies):
            entropy[idx] += (freq * log2(freq) * -1) if freq > 0 else 0

    return sum(entropy)


# Sample input
dna_sample = [
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG'
]
# sample = GreedyMotifSearchWithScore(dna_sample, 3, 5)
# print("Sample: ", sample)

# Test One
dna_1 = [
    'GCCCAA',
    'GGCCTG',
    'AACCTA',
    'TTCCTT',
]
# test_one = GreedyMotifSearch(dna_1, 3, 4)
# print("Test One: ", test_one)
# assert test_one == ['GCC', 'GCC', 'AAC', 'TTC']

# Test Two
dna_2 = [
    'GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC',
    'TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC',
    'TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT',
    'GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
    'GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT',
    'TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT',
    'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
    'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA',
]
# test_two = GreedyMotifSearch(dna_2, 5, 8)
# print("Test Two: ", test_two)
# assert test_two == [
#     'GAGGC',
#     'TCATC',
#     'TCGGC',
#     'GAGTC',
#     'GCAGC',
#     'GCGGC',
#     'GCGGC',
#     'GCATC',
# ]

# Test Four
dna_4 = [
    'GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA',
    'TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA',
    'GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA',
    'GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC',
    'AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA',
    'AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA',
    'GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA',
    'TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT',
]
test_four = GreedyMotifSearch(dna_4, 5, 8)
print("Test Four: ", test_four)
assert test_four == [
    'GCAGC',
    'TCATT',
    'GGAGT',
    'TCATC',
    'GCATC',
    'GCATC',
    'GGTAT',
    'GCAAC',
]