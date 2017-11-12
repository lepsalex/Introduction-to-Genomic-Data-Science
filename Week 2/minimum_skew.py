# Place your MinimumSkew() function here along with any subroutines you need.
# MinimumSkew() should return a list.
def MinimumSkew(Genome):
    skew_arr = Skew(Genome)
    sorted_vals = sorted(skew_arr)
    lowest_value = sorted_vals[0]
    min_points = []
    for index, skew in enumerate(skew_arr):
        if skew == lowest_value:
            min_points.append(index)
    return min_points

def Skew(genome):
    genome_length = len(genome)
    skew_arr = [0] * (genome_length + 1)
    for i in range(genome_length):
        if genome[i] == "G":
            skew_arr[i + 1] = skew_arr[i] + 1
        elif genome[i] == "C":
            skew_arr[i + 1] = skew_arr[i] - 1
        else:
            skew_arr[i + 1] = skew_arr[i]
    return skew_arr

genome = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
min_skew = MinimumSkew(genome)
print(min_skew)
