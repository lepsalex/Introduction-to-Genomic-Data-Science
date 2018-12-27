from random import uniform, randint
from functools import reduce

# Write a GibbsSampler() function here along with any subroutines that you need.
# GibbsSampler() should return a list of strings
def GibbsSampler(Dna, k, t, N):
    pseudocount = 1
    bestMotifs = GenerateRandomMotifs(Dna, k, t)
    for _ in range(N):
        i = randint(0, t - 1)
        reducedDna = Dna[:]
        extractedText = reducedDna.pop(i)
        profile = GetFrequencyMatrix(bestMotifs, pseudocount)
        newMotifs = Motifs(profile, reducedDna, k)
        newMotifs.insert(i, ProfileGeneratedString(extractedText, profile, k))
        if Score(newMotifs) < Score(bestMotifs):
            bestMotifs = newMotifs

    return bestMotifs

# Fill in your RepeatedGibbsSampler() function here. You should simply call
# GibbsSampler() n times and return the best Motifs. This function should 
# return a list of strings.
def RepeatedGibbsSampler(Dna, k, t, N, n):
    bestMotifs = GenerateRandomMotifs(Dna, k, t)
    for _ in range(n):
        newMotifs = GibbsSampler(Dna, k, t, N)
        if Score(newMotifs) < Score(bestMotifs):
            bestMotifs = newMotifs

    return bestMotifs

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    roll = uniform(0, 1)
    last = 0
    for k, v in Probabilities.items():
        if roll <= (last + v):
            return k
        else:
            last = last + v

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    total = sum(Probabilities.values())
    for k, v in Probabilities.items():
        Probabilities[k] = v / total
    return Probabilities


def GenerateRandomMotifs(Dna, k, t):
    return [x[y:y+k] for x, y in zip(Dna, [randint(0, len(Dna[0]) - k) for z in range(t)])]


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

def Motifs(Profile, Dna, k):
    return [ProfileMostProbableKmer(Dna[i], k, Profile) for i in range(len(Dna))]


# Reduce to kmer with highest score tuple, return the kmer string
def ProfileMostProbableKmer(text, k, profile):
    return reduce(lambda acc, curr: ProfileReduce(acc, curr, profile), [text[x:x + k] for x in range(0, len(text) - k + 1)], (text[0:k], 0))[0]


# Compare pm score and return higher as tuple (kmer, score)
def ProfileReduce(acc, curr, profile):
    pm = Pr(curr, profile)
    return (curr, pm) if pm > acc[1] else acc


# ComputeProfileMatrix(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753
def Pr(pattern, profile):
    return reduce(lambda x, y: x * y, (profile[nuc][idx] for idx, nuc in enumerate(pattern)))


def Score(motifs):
    pseudocount = 1
    conKmer = GetConcensusKmer(motifs, pseudocount)
    return reduce(lambda acc, curr: acc + HammingDistance(conKmer, curr), motifs, 0)


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

def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
        if char_1 != char_2:
            mismatches += 1
    return mismatches


# Sample
sample_dna = [
    'CAGAATCGGAGCCGCTTACTCTTTCCAAGTACCAGACAGACGACGCAGAA',
    'TCGGAGCCGCTTACTCTTTCCAAGTACCAGACAGACGATCCTGCGCAGAA',
    'TCGCACGCTTTCTTGCACACGCCTACAAACTACTCAGAAGGGCCCCAGAT',
    'AGGCTACGACGCTGCAGTGTCACTAATGCCTGTGAAAACTTCCGGTTTGC',
    'CGGGGTGTTGGTGAACGGCGCGCGGGCTGTCTTCCAGCCCCCACCCACTC',
    'CATTATGCTACCAATAGATGCCTGGCAACGCGTTGCAGACATTCACCTAT',
    'GGGCACTCCCTTCTGTCCCAACCTTCAGCAGATCCCCAACACTCAACGGG',
    'AAAAGCAATGCCCCCCACACTCCTGACAGTTATCAACGAGCTACCCGGCC',
    'GGCATGGCGTTCGCTTCTTGTCAGCTCGGAACTTGCTCATCGGGAATCAT',
    'GTGGACCTCTCGTTAAGCCTGCTACGATTGGCCATATATGCTGGTCCTGC',
]

sample = GibbsSampler(sample_dna, 5, 10, 1500)
print(sample)