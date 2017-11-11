# Place your ComputingFrequencies() function here along with any subroutines you need.
# ComputingFrequencies() should return a list.
def ComputingFrequencies(text, k):
    frequencies = [0] * ( (4**k) )
    for i in range(len(text) - k + 1):
        pattern = text[i:i+k]
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


print( ComputingFrequencies('ACGCGGCTCTGAAA', 2) )
