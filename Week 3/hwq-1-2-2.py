prob = 0.25 ** 9
occ = (1000 - 9 + 1) * prob
final = 500 * occ

print(final)

x = 1000
end = 1000 - 9 # after the "k-th" time we are repeating strings
y = 0
while x > end:
    y = y + (x // 9)
    x = x - 1

print(y * prob * 500)