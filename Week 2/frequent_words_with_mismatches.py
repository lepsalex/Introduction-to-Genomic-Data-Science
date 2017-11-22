# Write your FrequentWordsWithMismatches() function here, along with any subroutines you need.
# Your function should return a list.
def FrequentWordsWithMismatches(Text, k, d):
    FrequentPatterns = set()
    Neighborhoods = []
    NeighborhoodArray = []

    Index = []
    Count = []

    for i in range(len(Text) - k):
        Neighborhoods.append(Neighbors(Text[i:i + k], d))

    for Neighborhood in Neighborhoods:
        NeighborhoodArray = NeighborhoodArray + list(Neighborhood)

    for i in range(len(Neighborhoods)):
        Pattern = NeighborhoodArray[i]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)

    SortedIndex = sorted(Index)

    for i in range(len(Neighborhoods) - 1):
        if SortedIndex[i] == SortedIndex[i + 1]:
            Count[i + 1] = Count[i] + 1

    maxCount = max(Count)
    print(maxCount)

    for i in range(len(Neighborhoods)):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(SortedIndex[i], k)
            FrequentPatterns.add(Pattern)

    return FrequentPatterns


def Neighbors(Pattern, d):
    nucleotides = {"A", "C", "G", "T"}

    if d == 0:
        return set(Pattern)

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


# TEST
text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"

k = 4
d = 1

result = FrequentWordsWithMismatches(text, k, d)
print(result)
