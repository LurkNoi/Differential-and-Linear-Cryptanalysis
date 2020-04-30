# -*- coding: utf-8 -*-
"""
utils for linear cryptanalysis of the basic SPN cipher
"""
from itertools import product


def linear_approximation_table(sbox, absolute=False):
    """
    Return the linear approximation table for a given sbox
    """
    nrows = len(sbox)
    ncols = 1 << (max(sbox).bit_length())
    if nrows & (nrows - 1) != 0:
        raise TypeError("sbox length is not a power of 2")

    absolute_bias_table = [[0 for j in range(ncols)] for i in range(nrows)]
    for input_mask, output_mask in product(range(nrows), range(ncols)):
        total = 0
        for idx in range(nrows):
            input_masked = idx & input_mask
            output_masked = sbox[idx] & output_mask
            if (bin(input_masked).count("1") + bin(output_masked).count("1"))&1 == 0:
                total += 1
        count = total - (nrows//2)
        if absolute:
            count = abs(count)
        absolute_bias_table[input_mask][output_mask] = count

    return absolute_bias_table

def print_table(table, nonzero=False):
    """
    Print the linear approximation table
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
            v = table[input_mask][output_mask]
            if nonzero and v == 0:
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

    LAT = linear_approximation_table(SBOX)
    print_table(LAT)
