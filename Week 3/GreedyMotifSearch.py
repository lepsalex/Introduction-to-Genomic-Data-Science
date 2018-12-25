from functools import reduce

# Implement GreedyMotifSearch() below along with any subroutines that you need.
# Your function should return a list of strings.


def GreedyMotifSearch(Dna, k, t):
    return GreedyMotifSearchWithScore(Dna, k, t)[0]


def GreedyMotifSearchWithScore(Dna, k, t):
    bestMotifs = [x[:k] for x in Dna]
    bestMotifsScore = GetScore(GetConcensusKmer(bestMotifs), bestMotifs)
    firstLineKmers = [Dna[0][x:x + k] for x in range(0, len(Dna[0]) - k + 1)]
    linesToAnalyze = Dna[1:t]
    return reduce(lambda x, y: GmsReducer(x, y, linesToAnalyze, k), firstLineKmers, (bestMotifs, bestMotifsScore))


def GmsReducer(acc, curr, linesToAnalyze, k):
    # Pass to the next reducer the lines[1:t] and the starting motif
    currMotifs = reduce(lambda x, y: MotifReducer(
        x, y, k), linesToAnalyze, [curr])
    currScore = GetScore(GetConcensusKmer(currMotifs), currMotifs)
    return (currMotifs, currScore) if currScore < acc[1] else acc


def MotifReducer(acc, curr, k):
    # Create a profile from all motifs in our list
    # first pass will have the passed in starting motif
    profile = GetFrequencyMatrix(acc)
    return acc + [ProfileMostProbableKmer(curr, k, profile)]


# Reduce to kmer with highest score tuple, return the kmer string
def ProfileMostProbableKmer(text, k, profile):
    return reduce(lambda acc, curr: ProfileReduce(acc, curr, profile), [text[x:x + k] for x in range(0, len(text) - k + 1)], (text[0:k], 0))[0]


# Compare pm score and return higher as tuple (kmer, score)
def ProfileReduce(acc, curr, profile):
    pm = ComputeProfileMatrix(curr, profile)
    return (curr, pm) if pm > acc[1] else acc


# ComputeProfileMatrix(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753
def ComputeProfileMatrix(pattern, profile):
    return reduce(lambda x, y: x * y, (profile[nuc][idx] for idx, nuc in enumerate(pattern)))


def GetFrequencyMatrix(input):
    counts = GenerateEmptyCountMatrix(len(input[0]))

    for line in input:
        for idx, nuc in enumerate(line):
            counts[nuc][idx] += 1 / len(input)

    return counts


def GenerateEmptyCountMatrix(length):
    return {
        'A': [0] * length,
        'C': [0] * length,
        'G': [0] * length,
        'T': [0] * length
    }


def GetConcensusKmer(kmers):
    matrix = GetFrequencyMatrix(kmers)
    kmerLen = len(kmers[0])
    concensus = ''
    for i in range(kmerLen):
        char = ('', 0)
        for k, v in matrix.items():
            if v[i] > char[1]:
                char = (k, v[i])
        concensus += char[0]
    
    return concensus


def GetScore(conKmer, motifs):
    return reduce(lambda acc, curr: acc + HammingDistance(conKmer, curr), motifs, 0)


def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
        if char_1 != char_2:
            mismatches += 1
    return mismatches


# # Sample input
# dna_sample = [
#     'GGCGTTCAGGCA',
#     'AAGAATCAGTCA',
#     'CAAGGAGTTCGC',
#     'CACGTCAATCAC',
#     'CAATAATATTCG'
# ]
# sample = GreedyMotifSearchWithScore(dna_sample, 3, 5)
# print("Sample: ", sample)

# Test One
dna_1 = [
    'GCCCAA',
    'GGCCTG',
    'AACCTA',
    'TTCCTT',
]
test_one = GreedyMotifSearch(dna_1, 3, 4)
print("Test One: ", test_one)
assert test_one == ['GCC', 'GCC', 'AAC', 'TTC']

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
test_two = GreedyMotifSearch(dna_2, 5, 8)
print("Test Two: ", test_two)
assert test_two == [
    'GAGGC',
    'TCATC',
    'TCGGC',
    'GAGTC',
    'GCAGC',
    'GCGGC',
    'GCGGC',
    'GCATC',
]

# Test Three
dna_3 = [
    'GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT',
    'CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT',
    'GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT',
    'AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC',
    'ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG',
]
test_three = GreedyMotifSearch(dna_3, 6, 5)
print("Test Three: ", test_three)
assert test_three == [
    'GTGCGT',
    'GTGCGT',
    'GCGCCA',
    'GTGCCA',
    'GCGCCA',
]

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
