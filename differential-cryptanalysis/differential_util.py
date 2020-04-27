# -*- coding: utf-8 -*-
"""
utils for differential cryptanalysis of the basic SPN cipher
"""


def byte_xor(a: bytes, b: bytes) -> bytes:
    """
    return a XOR b
    """
    if len(a) != len(b):
        raise TypeError("Both a and b must have the same length.")
    return bytes([x^y for x, y in zip(a, b)])


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
    for output_diff in range(ncols):
        print("{:3x}".format(output_diff), end=' ')
    print()
    print(' ' + '-' * (4*ncols + 4))
    for input_diff in range(nrows):
        print("{:3x} |".format(input_diff), end=' ')
        for output_diff in range(ncols):
            print("{:3d}".format(table[input_diff][output_diff]), end=' ')
        print()


if __name__ == '__main__':
    from spn import SBOX

    DDTable = difference_distribution_table(SBOX)
    print_table(DDTable)
