# Week 3 Questions

"""
MotifEnumeration(Dna, k, d)
        Patterns ← an empty set
        for each k-mer Pattern in the first string in Dna
            for each k-mer Pattern’ differing from Pattern by at most d mismatches
                if Pattern' appears in each string from Dna with at most d mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns
"""

from functools import reduce

# Write your MotifEnumeration() function here along with any subroutines you need.
# This function should return a list of strings.
def MotifEnumeration(dna, k, d):
    Patterns = set()
    initLine = dna[0]
    for initialKmer in [initLine[x:x + k] for x in range(0, len(initLine) - k + 1)]:
        for neighborK in Neighbors(initialKmer, d):
            isValid = reduce(lambda acc, curr: False if acc == False else True if ApproximatePatternCount(curr, neighborK, d) > 0 else False, dna[1:], True)
            if (isValid):
                Patterns.add(neighborK)
    return sorted(Patterns)

def Neighbors(Pattern, d):
    nucleotides = {"A", "C", "G", "T"}

    if d == 0:
        return set([Pattern])

    if len(Pattern) == 1:
        return nucleotides

    Neighborhood = set()

    SuffixNeighbors = Neighbors(Pattern[1:], d)

    for text in SuffixNeighbors:
        if HammingDistance(Pattern[1:], text) < d:
            for nucleotide in nucleotides:
                Neighborhood.add(nucleotide + text)
        else:
            Neighborhood.add(Pattern[:1] + text)

    return Neighborhood

def ApproximatePatternCount(Text, Pattern, d):
    count = 0
    text_len = len(Text)
    pattern_len = len(Pattern)
    for i in range(text_len - pattern_len + 1):
        sample = Text[i:i+pattern_len]
        if HammingDistance(Pattern, sample) <= d:
            count += 1
    return count

def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
        if char_1 != char_2:
            mismatches += 1
    return mismatches


### test ###

dna1 = ['ACGT', 'ACGT', 'ACGT']
assert MotifEnumeration(dna1, 3, 0) == ['ACG', 'CGT']

dna2 = ['AAAAA', 'AAAAA', 'AAAAA']
assert MotifEnumeration(dna2, 3, 1) == ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'AGA', 'ATA', 'CAA', 'GAA', 'TAA']

dna3 = ['AACAA', 'AAAAA', 'AAAAA']
assert MotifEnumeration(dna3, 3, 0) == []

dna4 = ['AAAAA', 'AAAAA', 'AACAA']
assert MotifEnumeration(dna4, 3, 0) == []

dna5 = ['AAAAA', 'AAAAA', 'AAAAA']
assert MotifEnumeration(dna5, 3, 3) == ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]