# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []  # initializing list of positions
    # your code here
    text_len = len(Text)
    pattern_len= len(Pattern)
    for i in range(text_len - pattern_len + 1):
        sample = Text[i:i+pattern_len]
        if HammingDistance(Pattern, sample) <= d:
            positions.append(i)
    return positions


# Insert your Hamming distance function on the following line.
def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
            if char_1 != char_2:
                    mismatches += 1
    return mismatches


# # TEST
# text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
# pattern = "ATTCTGGA"

# TEST 3
text = "CCAAATCCCCTCATGGCATGCATTCCCGCAGTATTTAATCCTTTCATTCTGCATATAAGTAGTGAAGGTATAGAAACCCGTTCAAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGCGCCATAATCCAAACA"
pattern = "AATCCTTTCA"

# # TEST 4
# text = "CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC"
# pattern = "CCGTCATCC"

positions = ApproximatePatternMatching(text, pattern, 3)
print(positions)

