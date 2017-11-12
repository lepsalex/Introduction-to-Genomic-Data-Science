def Skew(genome):
    genome_length = len(genome)
    skew_arr = [0] * (genome_length + 1)
    for i in range(genome_length):
        if genome[i] == "G":
            skew_arr[i + 1] = skew_arr[i] + 1
        elif genome[i] == "C":
            skew_arr[i + 1] = skew_arr[i] - 1
        else :
            skew_arr[i + 1] = skew_arr[i]
    return skew_arr
            

genome = "GAGCCACCGCGATA"
print(' '.join(str(e) for e in Skew(genome)))
