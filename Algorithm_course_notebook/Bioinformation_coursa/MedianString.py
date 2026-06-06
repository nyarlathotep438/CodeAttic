from itertools import product
import sys
from pathlib import Path

def hamming_distance(p, q):
    return sum(1 for a, b in zip(p, q) if a != b)

def d(pattern, text):
    k = len(pattern)
    return min(hamming_distance(pattern, text[i:i+k]) for i in range(len(text) - k + 1))

def d_pattern_dna(pattern, dna):
    return sum(d(pattern, seq) for seq in dna)

def median_string(dna, k):
    best_pattern = None
    best_distance = float("inf")

    for pattern_tuple in product("ACGT", repeat=k):
        pattern = "".join(pattern_tuple)
        distance = d_pattern_dna(pattern, dna)

        if distance < best_distance:
            best_distance = distance
            best_pattern = pattern

    return best_pattern

if __name__ == "__main__":
    file_path = Path(sys.argv[1])
    with file_path.open() as f:
        k = int(f.readline().strip())
        dna = f.readline().strip().split()

    print(median_string(dna, k))
