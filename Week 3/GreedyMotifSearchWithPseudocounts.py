from functools import reduce

# Implement GreedyMotifSearchWithPseudocounts() below along with any subroutines that you need.
# Your function should include a pseudocount parameter that is added to each element of the count matrix.
# (When your function is called, the pseudocount parameter will be equal to 1.)


def GreedyMotifSearchWithPseudocounts(Dna, k, t, pseudocount = 1):
    return GMSPseudocountsWithScore(Dna, k, t, pseudocount)[0]


def GMSPseudocountsWithScore(Dna, k, t, pseudocount):
    bestMotifs = [x[:k] for x in Dna]
    bestMotifsScore = GetScore(GetConcensusKmer(
        bestMotifs, pseudocount), bestMotifs)
    firstLineKmers = [Dna[0][x:x + k] for x in range(0, len(Dna[0]) - k + 1)]
    linesToAnalyze = Dna[1:t]
    return reduce(lambda x, y: GmsReducer(x, y, linesToAnalyze, k, pseudocount), firstLineKmers, (bestMotifs, bestMotifsScore))


def GmsReducer(acc, curr, linesToAnalyze, k, pseudocount):
    # Pass to the next reducer the lines[1:t] and the starting motif
    currMotifs = reduce(lambda x, y: MotifReducer(
        x, y, k, pseudocount), linesToAnalyze, [curr])
    currScore = GetScore(GetConcensusKmer(currMotifs, pseudocount), currMotifs)
    return (currMotifs, currScore) if currScore < acc[1] else acc


def MotifReducer(acc, curr, k, pseudocount):
    # Create a profile from all motifs in our list
    # first pass will have the passed in starting motif
    profile = GetFrequencyMatrix(acc, pseudocount)
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


def GetFrequencyMatrix(input, pseudocount):
    counts = GenerateEmptyCountMatrix(len(input[0]), pseudocount)

    for line in input:
        for idx, nuc in enumerate(line):
            counts[nuc][idx] += 1 / len(input)

    return counts


def GenerateEmptyCountMatrix(length, pseudocount=1):
    return {
        'A': [pseudocount] * length,
        'C': [pseudocount] * length,
        'G': [pseudocount] * length,
        'T': [pseudocount] * length
    }


def GetConcensusKmer(kmers, pseudocount):
    matrix = GetFrequencyMatrix(kmers, pseudocount)
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


# Sample input
dna_sample = [
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG',
]
sample = GreedyMotifSearchWithPseudocounts(dna_sample, 3, 5, 1)
print("Sample: ", sample)
