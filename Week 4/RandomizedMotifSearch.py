from random import randint
from functools import reduce

# Import the random package, then write a function RandomizedMotifSearch() here along with any subroutines you need.
# RandomizedMotifSearch() should return a list of strings.
def RandomizedMotifSearch(Dna, k, t):
    return [x[y:y+k] for x, y in zip(Dna, [randint(0, len(Dna[0]) - k) for z in range(t)])]


# Then, write a function here called RepeatedRandomizedMotifSearch() that takes Dna, k, t, and a parameter
# N that returns the best collection of motifs resulting from running RandomizedMotifSearch()
# N times.  It should return a list of strings.
def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    pseudocount = 1
    bestMotifs = RandomizedMotifSearch(Dna, k, t)
    bestMotifsScore = GetScore(GetConcensusKmer(bestMotifs, pseudocount), bestMotifs)
    for _ in range(N):
        newMotifs = RandomizedMotifSearch(Dna, k, t)
        newMotifsScore = GetScore(GetConcensusKmer(newMotifs, pseudocount), newMotifs)
        if newMotifsScore < bestMotifsScore:
            bestMotifs = newMotifs
            bestMotifsScore = newMotifsScore

    return bestMotifs


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


def GetScore(conKmer, motifs):
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