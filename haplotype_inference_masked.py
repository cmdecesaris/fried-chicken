import numpy as np
import time
import sys
import os


def load_data(filename):
    data = np.loadtxt(filename, dtype=str, delimiter=" ")
    return data


def possible(h_list, genotype, haplotype, index):
    if len(genotype) == index:
        h_list[haplotype] = True
        return h_list

    if genotype[index] == "1":
        possible(h_list, genotype, haplotype + "0", index+1)
        possible(h_list, genotype, haplotype + "1", index+1)

    if genotype[index] == "0":
        possible(h_list, genotype, haplotype + "0", index+1)

    if genotype[index] == "2":
        possible(h_list, genotype, haplotype + "1", index+1)

    if genotype[index] == "*":
        possible(h_list, genotype, haplotype + "*", index+1)


def compute_pairs(h_list, genotype):
    pairs = []
    seen = {}
    for haplotype in h_list.keys():
        # Get the other pair
        # print len(genotype) == len(haplotype)

        # print "Considering haplotype", haplotype
        complement = compute_complement(genotype, haplotype)

        if h_list.get(complement) and not seen.get(complement) and not seen.get(haplotype):
            # Ex: If 101*1 has complement 100*1
            # We'll generate all the possible haplotype pairs here,
            # then remove the complement 100*1 because it won't match to anything else,
            # so we don't deal with it again.
            seen[complement] = True
            # print "seen list", seen
            pairs.extend(compute_starred_pairs(haplotype, complement))
        else:
            pass
            # print "Uh oh, haplotype", haplotype, "doesn't have a complement"

        seen[haplotype] = True

    # print "pairs found:", pairs
    # print "number of pairs:", len(pairs)
    return pairs


def expand_star(haplotype, combinations, index, partial_haplotype):
    if index == len(haplotype):
        combinations.append(partial_haplotype)
        return combinations

    if haplotype[index] == "*":
        expand_star(haplotype, combinations, index+1,
                    partial_haplotype + "0")
        expand_star(haplotype, combinations, index+1,
                    partial_haplotype + "1")
    else:
        expand_star(haplotype, combinations, index+1,
                    partial_haplotype + haplotype[index])


def compute_starred_pairs(haplotype, complement):
    partners1 = []
    partners2 = []
    pairs = []
    expand_star(haplotype, partners1, 0, "")
    expand_star(complement, partners2, 0, "")

    # print "Partners 1 and 2", partners1, partners2
    for husband in partners1:
        for wife in partners2:
            pairs.append((husband, wife))
    # print "I found these pairs", pairs
    return pairs


def compute_complement(genotype, haplotype):
    # print(haplotype, genotype, type(haplotype), type(genotype))
    complement = ""
    for i in range(len(haplotype)):
        SNP1 = haplotype[i]
        SNP2 = genotype[i]
        if SNP1 != "*" and SNP2 != "*":
            complement += str(int(SNP2) - int(SNP1))
        else:
            complement += "*"
    # print "Complement is ", complement
    return complement


def EM_algorithm(genotype_list, iterations):  # 2-D array
    # genotype_list = genotype_list.astype(int)
    genotype_to_pair = {}  # genotypes map to compatible pairs of haplotypes
    unique_haplotypes = {}  # haplotypes map to a probability
    hap_pair_to_expect = {}  # haplotype pair maps to expected value

    ret_dict = {}  # genotype to haplotype pair
    ret_list = []

    n = genotype_list.shape[0]  # gets the number of genotypes

    for genotype in genotype_list:
        # print "Genotype list is:\n", genotype_list, type(genotype_list), type(genotype_list[0])
        # created unique_haplotypes
        h_list = {}
        possible(h_list, genotype, "", 0)

        # print "Haplotype list is:", h_list
        # Convert genotype to a string
        genotype = ''.join(str(e) for e in genotype)

        # pairs is list of tuples of strings
        pairs = compute_pairs(h_list, genotype)
        genotype_to_pair[genotype] = pairs
        for pair in pairs:
            unique_haplotypes[pair[0]] = 0
            unique_haplotypes[pair[1]] = 0

    # intialize initial probabilites
    num = len(unique_haplotypes)
    unique_haplotypes = dict.fromkeys(unique_haplotypes, 1/float(num))

    # print(unique_haplotypes)

    # Want to get a list of tuples
    for genotypes in genotype_to_pair.keys():
        for i, element in enumerate(genotype_to_pair[genotypes]):
            genotype_to_pair[genotypes][i] = tuple(element)

    # Now done with initial setup

    # print(genotype_to_pair)
    # print(unique_haplotypes)
    # print(hap_pair_to_expect)

    for _ in range(iterations):
        # now map the pair to a excpectation
        # get sum of compatible pairs C(g)
        # This is the E step
        for genotype in genotype_to_pair.keys():
            denominator_sum = 0
            for pair in genotype_to_pair[genotype]:
                denominator_sum += (unique_haplotypes[pair[0]]) * \
                    (unique_haplotypes[pair[1]])
            # now we have the denominator sum

            # calculate expected for every pair
            for pair in genotype_to_pair[genotype]:
                hap_pair_to_expect[pair] = (
                    (unique_haplotypes[pair[0]]) * (unique_haplotypes[pair[1]])) / denominator_sum

        # This is the M step
        # now update ph's
        # go through all unique haps and search for where it appears in hap_pair_to_expect
        for hap in unique_haplotypes.keys():
            # serach for haplotype in hap_paor
            expt = 0
            for genotype in genotype_to_pair.keys():
                for pair in genotype_to_pair[genotype]:
                    if hap in pair:
                        expt += hap_pair_to_expect[pair]
            unique_haplotypes[hap] = expt/(2*float(n))

        # Check if algorithm has converged and return the correct haplotypes
        # There is a bug here
        count = 0
        for genotype in genotype_to_pair.keys():
            for pair in genotype_to_pair[genotype]:
                if hap_pair_to_expect[pair] == 1:
                    ret_dict[genotype] = pair
                    count += 1
        if (count == len(genotype_to_pair)):
            break

    if (len(genotype_to_pair) != len(ret_dict)):
        # some haplotypes did not converge to expectation 1
        # get the max expected pair
        for genotype in genotype_to_pair.keys():
            if (genotype not in ret_dict.keys()):
                max_val = 0
                for pair in genotype_to_pair[genotype]:
                    if (hap_pair_to_expect[pair] >= max_val):
                        max_val = hap_pair_to_expect[pair]
                        max_pair = pair
                ret_dict[genotype] = max_pair

    # print(hap_pair_to_expect)
    # print(genotype_to_pair)
    # print(len(ret_dict))
    # print(len(genotype_to_pair))
    # Now get in correct formate
    for i in range(genotype_list.shape[0]):  # for every genotype
        genotype = ''.join(str(e)
                           for e in genotype_list[i])  # now genotype is string
        ret_list.append(ret_dict[genotype][0])
        ret_list.append(ret_dict[genotype][1])

    return ret_list


def EM(X_1, bucket_size, iterations):
    num_snps = X_1.shape[1]
    num_ppl = X_1.shape[0]
    print "Input data shape:", num_snps, num_ppl

    # Dimensions is (num_ppl * 2) by num_snps. Will transpose later
    final_phase_list = [""]*num_ppl*2

    start = time.clock()
    window_count = 0

    i = 0
    while i < num_snps:
        window = X_1[:num_ppl, i:i+bucket_size]  # this may go out of range?
        time_passed = (time.clock()-start)/60
        window_count += 1

        # EM helper does the whole EM algorithm to find most likely phasing
        # Returns a phase list of num_ppl * 2 by bucket size
        # Append buckets together to create final phase list of
        # num_ppl * 2 by nums_SNPs

        phase_list = EM_algorithm(window, iterations)
        # print "Phase list dimensions:", len(phase_list), len(phase_list[0])
        for j, phase in enumerate(phase_list):
            final_phase_list[j] = final_phase_list[j] + phase_list[j]

        i = i + bucket_size

        # Output stats every 100 windows (1000 SNPs if window is 10)
        if window_count % 100 == 0 and i != 0:
            print '{} SNPs looked at'.format(i), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/(i)*num_snps - time_passed
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)

    output = np.array(final_phase_list).transpose()
    return output


def write_to_file(filename, output, bucket_size, num_EM_iterations):
    fout = open("output/sol_" + str(bucket_size) +
                "_" + str(num_EM_iterations) + "_" + os.path.basename(filename), 'w+')

    for i in range(0, len(output[0])):
        for person in output:
            fout.write("%s " % person[i])
        fout.write("\n")

    fout.close()


if __name__ == "__main__":
    # This is the Python main function.
    if len(sys.argv) != 4:
        print('Usage: python ' +
              sys.argv[0] + ' <data_to_load_path> <bucket_size> <num_EM_iterations>')
        quit()

    filename = sys.argv[1]
    bucket_size = int(sys.argv[2])
    num_EM_iterations = int(sys.argv[3])

    # 50 individuals, 39496 or 75503 SNPS
    X_1 = load_data(filename).transpose()

    output = EM(X_1, bucket_size, num_EM_iterations)

    # One column is one haplotype, with length num of SNPs.
    # Two columns next to each other adds up to one person's genotype.
    # Result is an array of num_snps by (num_ppl*2)
    write_to_file(filename, output, bucket_size, num_EM_iterations)
