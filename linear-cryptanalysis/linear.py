# -*- coding: utf-8 -*-
"""
implementation of linear cryptanalysis of the basic SPN cipher

use example mentioned in section 3.4
    U_{4,6}  ^ U_{4,8}  ^ U_{4,14}  ^ U_{4,16}  ^ P_{5}  ^ P_{7}  ^ P_{8} = 0
    with a bias of 1/32
"""
import os

from spn import byte_to_int, basicSPN


KEY = os.urandom(10)
SPN_CIPHER = basicSPN(KEY)
SBOX_INV = SPN_CIPHER.sbox_dct['sbox_inv']
print(f'key: {KEY.hex()}')

# a known plaintext attack
BS = 2 # block size (in bytes)
NPAIRS = 10000 # number of plaintext-ciphertext pairs
plain = os.urandom(BS * NPAIRS)
cipher = SPN_CIPHER.encrypt(plain)
data_pairs = [(plain[i*BS : (i+1)*BS], cipher[i*BS : (i+1)*BS])
              for i in range(NPAIRS)]

partial_subkey_nbits = 8
count_table = [-NPAIRS//2 for _ in range(1<<partial_subkey_nbits)]
for pt, ct in data_pairs:
    pt_i, ct_i = byte_to_int(pt), byte_to_int(ct)
    # p_partial_sum = p_5 ^ p_7 ^ p_8
    p_partial_sum = ((pt_i >> 11) ^ (pt_i >> 9) ^ (pt_i >> 8))&0b1
    ct_5_to_8 = (ct_i >> 8)&0b1111
    ct_13_to_16 = ct_i & 0b1111
    for partial_subkey in range(1<<partial_subkey_nbits):
        k5_5_to_8 = (partial_subkey >> 4)&0b1111
        k5_13_to_16 = partial_subkey & 0b1111

        # reverse last round key mixing
        v4_5_to_8 = ct_5_to_8 ^ k5_5_to_8
        v4_13_to_16 = ct_13_to_16 ^ k5_13_to_16

        # reverse last round substitution
        u4_5_to_8 = SBOX_INV[v4_5_to_8]
        u4_13_to_16 = SBOX_INV[v4_13_to_16]

        # count how many times the eqution holds
        u4_6 = (u4_5_to_8 >> 2)&0b1
        u4_8 = u4_5_to_8 & 0b1
        u4_14 = (u4_13_to_16 >> 2)&0b1
        u4_16 = u4_13_to_16 & 0b1
        if u4_6 ^ u4_8 ^ u4_14 ^ u4_16 ^ p_partial_sum == 0:
            count_table[partial_subkey] += 1

bias_table = [{'partial_subkey': partial_subkey,
               'bias': abs(count_table[partial_subkey])}
              for partial_subkey in range(1<<partial_subkey_nbits)]


KEY_i = byte_to_int(KEY)
target_partial_subkey = ((KEY_i >> 4)&0b11110000) + (KEY_i & 0b1111)
print(f'Target partial subkey: {target_partial_subkey:02x} ' \
      f'with bias: {1/32.0:.4f}')
for prob in sorted(bias_table, key=lambda x: x['bias'], reverse=True)[:10]:
    print(f"  guess {prob['partial_subkey']:02x} " \
          f"with bias {prob['bias']/float(NPAIRS):.4f}")
