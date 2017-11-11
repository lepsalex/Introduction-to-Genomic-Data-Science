# fill in your ReverseComplement() function here along with any subroutines that you need.
def ReverseComplement(Pattern):
    mapping = {
        'A' : 'T',
        'C' : 'G',
        'G' : 'C',
        'T' : 'A'
    }

    complement = ""
    for char in Pattern:
        complement += mapping[char]

    return complement[::-1]


print(ReverseComplement('TTGTGTC'))
