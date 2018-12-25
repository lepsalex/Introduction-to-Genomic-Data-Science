from functools import reduce

def ComputeProfileMatrix(pattern, profile):
    return reduce(lambda x, y: x * y, (profile[nuc][idx] for idx, nuc in enumerate(pattern)), 1)

matrix = {
    'A': [.2, .2, 0, 0, 0, 0, .9, .1, .1, .1, .3, 0],
    'C': [.1, .6, 0, 0, 0, 0, 0, .4, .1, .2, .4, .6],
    'G': [0, 0, 1, 1, .9, .9, .1, 0, 0, 0, 0, 0],
    'T': [.7, .2, 0, 0, .1, .1, 0, .5, .8, .7, .3, .4]
}

test = ComputeProfileMatrix('ACGGGGATTACC', matrix)

print(test)
