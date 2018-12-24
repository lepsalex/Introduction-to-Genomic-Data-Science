from itertools import product

def GenerateAllKmersOfLength(k):
    dna = ['A', 'C', 'G', 'T']
    return [''.join(i) for i in product(dna, repeat = k)]

print(GenerateAllKmersOfLength(3))