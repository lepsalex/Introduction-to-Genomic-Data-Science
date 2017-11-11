# fill in your NumberToPattern() function here along with any subroutines you need.
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
