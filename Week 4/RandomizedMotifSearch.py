from random import randint
from functools import reduce

# Import the random package, then write a function RandomizedMotifSearch() here along with any subroutines you need.
# RandomizedMotifSearch() should return a list of strings.
def RandomizedMotifSearch(Dna, k, t):
    pseudocount = 1
    bestMotifs = GenerateRandomMotifs(Dna, k, t)
    while True:
        profile = GetFrequencyMatrix(bestMotifs, pseudocount)
        newMotifs = Motifs(profile, Dna, k)
        if Score(newMotifs) < Score(bestMotifs):
            bestMotifs = newMotifs
        else:
            return bestMotifs

# Then, write a function here called RepeatedRandomizedMotifSearch() that takes Dna, k, t, and a parameter
# N that returns the best collection of motifs resulting from running RandomizedMotifSearch()
# N times.  It should return a list of strings.
def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    bestMotifs = GenerateRandomMotifs(Dna, k, t)
    for _ in range(N):
        newMotifs = RandomizedMotifSearch(Dna, k, t)
        if Score(newMotifs) < Score(bestMotifs):
            bestMotifs = newMotifs

    return bestMotifs


def GenerateRandomMotifs(Dna, k, t):
    return [x[y:y+k] for x, y in zip(Dna, [randint(0, len(Dna[0]) - k) for z in range(t)])]


def Motifs(Profile, Dna, k):
    # insert your code here
    return [ProfileMostProbableKmer(Dna[i], k, Profile) for i in range(len(Dna))]


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


def Score(motifs):
    pseudocount = 1
    conKmer = GetConcensusKmer(motifs, pseudocount)
    return reduce(lambda acc, curr: acc + HammingDistance(conKmer, curr), motifs, 0)


def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
        if char_1 != char_2:
            mismatches += 1
    return mismatches


# Sample
sample_dna = [
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
]

sample = RepeatedRandomizedMotifSearch(sample_dna, 8, 5, 1000)
print(sample)