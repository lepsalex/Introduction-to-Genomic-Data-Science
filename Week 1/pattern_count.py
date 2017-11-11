haystack = 'CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC'
needle = 'CGCG'

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count += 1

    return count

# print(PatternCount(haystack, needle))

text = 'CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA'
kmer = 3

# fill in your FrequentWords() function here along with any subroutines you need.
def FrequentWords(Text, k):
    counts = {}
    max_count=  0
    max_kmers = set()
    for i in range(len(Text) - k):
        pattern = Text[i:i + k]
        counts[pattern] = PatternCount(Text, pattern)
        if counts[pattern] > max_count:
            max_count = counts[pattern]
        
    for kmer, count in counts.items():
        if count == max_count:
            max_kmers.add(kmer)

    return max_kmers

print(FrequentWords(text, kmer))
