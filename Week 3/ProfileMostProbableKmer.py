from functools import reduce

# Write your ProfileMostProbableKmer() function here along with any subroutines you need.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats

# Reduce to kmer with highest score tuple, return the kmer string
def ProfileMostProbableKmer(text, k, profile):
    return reduce(lambda acc, curr: ProfileReduce(acc, curr, profile), [
        text[x:x + k] for x in range(0, len(text) - k + 1)], (text[0:k], 0))[0]

# Compare pm score and return higher as tuple (kmer, score)
def ProfileReduce(acc, curr, profile):
    pm = ComputeProfileMatrix(curr, profile)
    return (curr, pm) if pm > acc[1] else acc

# ComputeProfileMatrix(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753
def ComputeProfileMatrix(pattern, profile):
    return reduce(lambda x, y: x * y, (profile[nuc][idx] for idx, nuc in enumerate(pattern)))


matrix_1 = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}

sample = ProfileMostProbableKmer(
    'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, matrix_1)
print(sample)