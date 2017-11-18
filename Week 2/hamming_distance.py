# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
            if char_1 != char_2:
                    mismatches += 1
    return mismatches


# TEST
string_1 = "GGGCCGTTGGT"
string_2 = "GGACCGTTGAC"

mismatches = HammingDistance(string_1, string_2)

print(mismatches)
