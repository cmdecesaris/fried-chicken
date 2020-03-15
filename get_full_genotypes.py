import sys
import numpy as np
from collections import defaultdict


def generate_compatible(genotype_samples, starred_positions):
    compatible = defaultdict(list)
    frequency = defaultdict(float)
    for idx, genotype in enumerate(genotype_samples):
        starred_list = starred_positions[genotype]
        temp = [genotype]


if __name__ == "__main__":
    if not len(sys.argv) == 2:
        sys.exit('Usage: python3 get_full_genotypes.py <masked_file.txt>')

    masked_file = sys.argv[1]

    f = open(masked_file,  "r")
    lines = f.readlines()
    num_samples = len(lines[0].strip("\n").split(" "))
    genotype_samples = np.full(num_samples, "", dtype=object)
    starred_positions = defaultdict(list)
    for line_idx, line in enumerate(lines):
        # print(line_idx)
        line = line.strip("\n").split(" ")
        #print("Currently considering line", line)
        for idx, value in enumerate(line):
            genotype_samples[idx] = genotype_samples[idx] + value
            if value == "*":
                starred_positions[idx].append(line_idx)
    max_count = 0
    for l in starred_positions.items():
        print(len(l[1]), type(l))
        if len(l[1]) > max_count:
            max_count = len(l[1])
    print("Max number of stars is", max_count)

    # print(genotype_samples)
    # print(starred_positions)
    # print(len(genotype_samples), genotype_samples[0])

    generate_compatible(genotype_samples, starred_positions)
