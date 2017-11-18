# Place your Neighbors() function here, along with any subroutines you need.
# Neighbors() should return a list.
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


# TEST
pattern = "ACG"
output = Neighbors(pattern, 1)
print(output)
