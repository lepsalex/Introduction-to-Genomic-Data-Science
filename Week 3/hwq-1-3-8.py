# Entropy Matric Calc
from math import log2

inputList = [
    ['T', 'C', 'G', 'G', 'G', 'G', 'G','T','T','T','T','T'],
    ['C', 'C', 'G', 'G', 'T', 'G', 'A','C','T','T','A','C'],

    ['A', 'C', 'G', 'G', 'G', 'G', 'A','T','T','T','T','C'],
    ['T', 'T', 'G', 'G', 'G', 'G', 'A','C','T','T','T','T'],

    ['A', 'A', 'G', 'G', 'G', 'G', 'A','C','T','T','C','C'],
    ['T', 'T', 'G', 'G', 'G', 'G', 'A','C','T','T','C','C'],
    
    ['T', 'C', 'G', 'G', 'G', 'G', 'A','T','T','C','A','T'],
    ['T', 'C', 'G', 'G', 'G', 'G', 'A','T','T','C','C','T'],

    ['T', 'A', 'G', 'G', 'G', 'G', 'A','A','C','T','A','C'],
    ['T', 'C', 'G', 'G', 'G', 'T', 'A','T','A','A','C','C']
]

def generateEmptyCountMatrix(length):
    return {
    'A': [0] * length,
    'C': [0] * length,
    'G': [0] * length,
    'T': [0] * length
}

def getFrequencyMatrix(input):
    counts = generateEmptyCountMatrix(len(input[0]))

    for line in input:
        for idx, nuc in enumerate(line):
            counts[nuc][idx] += 1 / len(input)

    return counts

def computeEntropyMotif(freqMatrix):
    entropy = [0] * len(list(freqMatrix.values())[0])
    
    for frequencies in freqMatrix.values():
        for idx, freq in enumerate(frequencies):
            entropy[idx] += (freq * log2(freq) * -1) if freq > 0 else 0

    return sum(entropy)

freqMat = getFrequencyMatrix(inputList)
ent = computeEntropyMotif(freqMat)

print(ent)