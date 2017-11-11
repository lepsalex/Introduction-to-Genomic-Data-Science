# fill in your PatternToNumber() function here along with any subroutines you need.
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


print(PatternToNumber('ATGCAA'))
