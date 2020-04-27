# -*- coding: utf-8 -*-
"""
utils for differential cryptanalysis of the basic SPN cipher
"""


def difference_distribution_table(sbox):
    """
    Return the difference distribution table for a given sbox
    """
    nrows = len(sbox)
    ncols = 1 << (max(sbox).bit_length())
    if nrows & (nrows - 1) != 0:
        raise TypeError("sbox length is not a power of 2")

    DDTable = [[0 for j in range(ncols)] for i in range(nrows)]
    for i in range(nrows):
        S_i = sbox[i]
        for input_diff in range(nrows):
            output_diff = S_i ^ sbox[i^input_diff]
            DDTable[input_diff][output_diff] += 1

    return DDTable

def print_table(table):
    """
    Print the difference distribution table
    """
    nrows = len(table)
    ncols = len(table[0])
    print("    |", end=' ')
    for output_mask in range(ncols):
        print("{:3x}".format(output_mask), end=' ')
    print()
    print(' ' + '-' * (4*ncols + 4))
    for input_mask in range(nrows):
        print("{:3x} |".format(input_mask), end=' ')
        for output_mask in range(ncols):
            print("{:3d}".format(table[input_mask][output_mask]), end=' ')
        print()


if __name__ == '__main__':
    from spn import SBOX

    DDTable = difference_distribution_table(SBOX)
    print_table(DDTable)