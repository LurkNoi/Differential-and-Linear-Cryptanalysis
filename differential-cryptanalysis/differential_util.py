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

def difference_distribution_table(sbox, truncated=False):
    """
    Return the difference distribution table for a given sbox
    """
    nrows = len(sbox)
    ncols = 1 << (max(sbox).bit_length())
    if nrows & (nrows - 1) != 0:
        raise TypeError("sbox length is not a power of 2")

    diff_dist_table = [[0 for j in range(ncols)] for i in range(nrows)]
    for i in range(nrows):
        s_i = sbox[i]
        for input_diff in range(nrows):
            output_diff = s_i ^ sbox[i^input_diff]
            diff_dist_table[input_diff][output_diff] += 1

    if truncated:
        truncated_DDT = [
            [1 if diff_dist_table[i][j] else 0 for j in range(ncols)]
            for i in range(nrows)
        ]
        return truncated_DDT

    return diff_dist_table

def print_table(table, nonzero=False):
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
            v = table[input_diff][output_diff]
            if nonzero and (v == 0):
                print(' '*3, end=' ')
            else:
                print("{:3d}".format(v), end=' ')
        print()


if __name__ == '__main__':
    SBOX = [
        0x0E, 0x04, 0x0D, 0x01,
        0x02, 0x0F, 0x0B, 0x08,
        0x03, 0x0A, 0x06, 0x0C,
        0x05, 0x09, 0x00, 0x07,
    ]

    DDT = difference_distribution_table(SBOX)
    print_table(DDT)
    DDT = difference_distribution_table(SBOX, truncated=True)
    print_table(DDT, nonzero=True)
