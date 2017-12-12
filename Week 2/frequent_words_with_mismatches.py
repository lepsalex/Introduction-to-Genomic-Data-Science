# Write your FrequentWordsWithMismatches() function here, along with any subroutines you need.
# Your function should return a list.
def FrequentWordsWithMismatches(Text, k, d):
    FrequentPatterns = set()
    Neighborhoods = []
    NeighborhoodArray = []
    
    Index = []
    Count = [0] * ((4**k))

    for i in range(len(Text) - k + 1):
        Neighborhoods.append(Neighbors(Text[i:(i + k)], d))

    for Neighborhood in Neighborhoods:
        for Pattern in Neighborhood:
            NeighborhoodArray.append(Pattern)

    for i in range(len(NeighborhoodArray)):
        Pattern = NeighborhoodArray[i]
        PatternNumber = PatternToNumber(Pattern)
        Count[PatternNumber] = 1
        Index.append(PatternNumber)

    SortedIndex = sorted(Index)
    for i in range(len(SortedIndex) - 1):
        if SortedIndex[i] == SortedIndex[i + 1]:
            patternNum = SortedIndex[i]
            Count[patternNum] = Count[patternNum] + 1

    maxCount = max(Count)

    for i in range(len(Count)):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.add(Pattern)

    return FrequentPatterns


def Neighbors(Pattern, d):
    nucleotides = {"A", "C", "G", "T"}

    if d == 0:
        return {Pattern}

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


def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
            if char_1 != char_2:
                    mismatches += 1
    return mismatches


def PatternToNumber(pattern):
    if len(pattern) is 0:
        return 0
    symbol = LastSymbol(pattern)
    prefix = Prefix(pattern)
    return 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)


def LastSymbol(pattern):
    return pattern[-1:]


def Prefix(pattern):
    return pattern[:-1]


def SymbolToNumber(symbol):
    mapping = ['A', 'C', 'G', 'T']
    return mapping.index(symbol)


def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    prefixIndex = index // 4
    r = index % 4
    symbol = NumberToSymbol(r)
    return NumberToPattern(prefixIndex, k - 1) + symbol


def NumberToSymbol(num):
    mapping = ['A', 'C', 'G', 'T']
    return mapping[num]


# TEST 0
text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k = 4
d = 1
assert FrequentWordsWithMismatches(text, k, d) == {'ATGC', 'ATGT', 'GATG'}

# TEST 1
text = "AAAAAAAAAA"
k = 2
d = 1
assert FrequentWordsWithMismatches(text, k, d) == {'AA', 'AC', 'AG', 'CA', 'AT', 'GA', 'TA'}

# TEST 2
text = "AGTCAGTC"
k = 4
d = 2
assert FrequentWordsWithMismatches(text, k, d) == {'TCTC', 'CGGC', 'AAGC', 'TGTG', 'GGCC', 'AGGT', 'ATCC', 'ACTG', 'ACAC', 'AGAG', 'ATTA', 'TGAC', 'AATT', 'CGTT', 'GTTC', 'GGTA', 'AGCA', 'CATC'}

# TEST 3
text = "AATTAATTGGTAGGTAGGTA"
k = 4
d = 0
assert FrequentWordsWithMismatches(text, k, d) == {'GGTA'}

# TEST 4
text = "ATA"
k = 3
d = 1
# print(FrequentWordsWithMismatches(text, k, d))
assert FrequentWordsWithMismatches(text, k, d) == {'GTA', 'ACA', 'AAA', 'ATC', 'ATA', 'AGA', 'ATT', 'CTA', 'TTA', 'ATG'}

# TEST 5
text = "AAT"
k = 3
d = 0
# print(sorted(FrequentWordsWithMismatches(text, k, d)))
assert FrequentWordsWithMismatches(text, k, d) == {'AAT'}

# TEST 5
text = "TAGCG"
k = 2
d = 1
assert FrequentWordsWithMismatches(text, k, d) == {'GG', 'TG'}
