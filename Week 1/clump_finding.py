# place your ClumpFinding() function here along with any subroutines you need.
def ClumpFinding(genome, k, L, t):
    frequent_patterns = set()
    clump = [0] * ((4**k))
    text = genome[0:L]
    frequency_array = ComputingFrequencies(text, k)
    
    for i in range(4**k):
        if frequency_array[i] >= t:
            clump[i] = 1

    for i in range(1, len(genome) - L):

        first_pattern = genome[(i-1):(i+k-1)]
        index = PatternToNumber(first_pattern)
        frequency_array[index] = frequency_array[index] - 1

        if frequency_array[index] >= t:
            clump[index] = 1
        
        last_pattern = genome[(i+L-k): i+L]
        index = PatternToNumber(last_pattern)
        frequency_array[index] = frequency_array[index] + 1
        
        if frequency_array[index] >= t:
            clump[index] = 1
    
    for i in range(4**k):
        if clump[i] == 1:
            pattern = NumberToPattern(i, k)
            frequent_patterns.add(pattern)
    
    return frequent_patterns


def ComputingFrequencies(text, k):
    frequencies = [0] * ((4**k))
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        j = PatternToNumber(pattern)
        frequencies[j] = frequencies[j] + 1
    return frequencies


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


# # Test
# genome = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
# output = ClumpFinding(genome, 5, 50, 4)
# print(output)

# CGACA GAAGA


# Test 2
# genome = "AAAACGTCGAAAAA"
# output = ClumpFinding(genome, 2, 4, 2)
# print(output)

# # # Test 3 (Large Data)
# with open("E_coli.txt", "r") as genomeFile:
#     genome = genomeFile.read().replace("\n", "")
#     output = ClumpFinding(genome, 9, 500, 3)
#     print(len(output))
