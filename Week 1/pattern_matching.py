# fill in your PatternMatching() function along with any subroutines that you need.
def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome)):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(str(i))
    
    return positions

with open("Vibrio_cholerae.txt", "r") as genomeFile:
    genome = genomeFile.read().replace("\n", "")
    pattern = "CTTGATCAT"
    positions = PatternMatching(pattern, genome)
    print(" ".join(positions))
