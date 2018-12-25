from itertools import product
from sys import maxsize

def MedianString(dna, k):
    distance = maxsize
    median = ''

    for kmer in GenerateAllKmersOfLength(k):
        dp = DistanceBetweenPatternAndStrings(kmer, dna)
        if distance > dp:
            distance = dp
            median = kmer
    
    return median


def DistanceBetweenPatternAndStrings(pattern, dna):
    k = len(pattern)
    distance = 0

    for line in dna:
        hd = maxsize
        for kmer in [line[x:x + k] for x in range(0, len(line) - k + 1)]:
            patternHd = HammingDistance(pattern, kmer)
            if hd > patternHd:
                hd = patternHd
        distance += hd

    return distance


def GenerateAllKmersOfLength(k):
    dna = ['A', 'C', 'G', 'T']
    return [''.join(i) for i in product(dna, repeat=k)]


def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
        if char_1 != char_2:
            mismatches += 1
    return mismatches
