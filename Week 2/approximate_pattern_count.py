def ApproximatePatternCount(Text, Pattern, d):
    count = 0
    text_len = len(Text)
    pattern_len= len(Pattern)
    for i in range(text_len - pattern_len + 1):
        sample = Text[i:i+pattern_len]
        if HammingDistance(Pattern, sample) <= d:
            count += 1
    return count


# Insert your Hamming distance function on the following line.
def HammingDistance(p, q):
    mismatches = 0
    for char_1, char_2 in zip(p, q):
            if char_1 != char_2:
                    mismatches += 1
    return mismatches


# TEST
text = "TTTAGAGCCTTCAGAGG"
pattern = "GAGG"
count = ApproximatePatternCount(text, pattern, 2)
print(count)
