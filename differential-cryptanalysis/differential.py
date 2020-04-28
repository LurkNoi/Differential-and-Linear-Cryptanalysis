# -*- coding: utf-8 -*-
r"""
implementation of differential cryptanalysis of the basic SPN cipher

use example mentioned in section 4.3
    delta P = delta U1 = [0000 1011 0000 0000] = b"\x0b\x00",
    delta U4 = [0000 0110 0000 0110] = b"\x06\x06"
    with a bias of 27/1024
"""
import os

from spn import byte_to_int, basicSPN
from differential_util import byte_xor

KEY = os.urandom(10)
SPN_CIPHER = basicSPN(KEY)
SBOX_INV = SPN_CIPHER.sbox_dct['sbox_inv']
print(SPN_CIPHER)

# a choosen plaintext attack
BS = 2 # block_size
NPAIRS = 5000 # number of pairs of ciphertexts with related plaintext
plain_1 = os.urandom(BS * NPAIRS)
cipher_1 = SPN_CIPHER.encrypt(plain_1)
plain_2 = byte_xor(plain_1, b"\x0b\x00" * NPAIRS)
cipher_2 = SPN_CIPHER.encrypt(plain_2)
data_pairs = [(cipher_1[i*BS : (i+1)*BS], cipher_2[i*BS : (i+1)*BS])
              for i in range(NPAIRS)]

partial_subkey_nbits = 8
count_table = [0 for _ in range(1<<partial_subkey_nbits)]
for data_pair_i in data_pairs:
    ct = list(map(byte_to_int, data_pair_i))
    ct_5_to_8 = [(ct[0] >> 8)&0b1111,
                 (ct[1] >> 8)&0b1111]
    ct_13_to_16 = [ct[0] & 0b1111,
                   ct[1] & 0b1111]
    for partial_subkey in range(1<<partial_subkey_nbits):
        k5_5_to_8 = (partial_subkey >> 4)&0b1111
        k5_13_to_16 = partial_subkey & 0b1111

        # reverse last round key mixing
        v4_5_to_8 = [ct_5_to_8[0] ^ k5_5_to_8,
                     ct_5_to_8[1] ^ k5_5_to_8]
        v4_13_to_16 = [ct_13_to_16[0] ^ k5_13_to_16,
                       ct_13_to_16[1] ^ k5_13_to_16]

        # reverse last round substitution
        u4_5_to_8 = [SBOX_INV[v4_5_to_8[0]],
                     SBOX_INV[v4_5_to_8[1]]]
        u4_13_to_16 = [SBOX_INV[v4_13_to_16[0]],
                       SBOX_INV[v4_13_to_16[1]]]

        # count how many times the expected difference occurs
        delta_u4_5_to_8 = u4_5_to_8[0] ^ u4_5_to_8[1]
        delta_u4_13_to_16 = u4_13_to_16[0] ^ u4_13_to_16[1]
        if 0x6 == delta_u4_5_to_8 == delta_u4_13_to_16:
            count_table[partial_subkey] += 1

bias_table = [{'partial_subkey': partial_subkey,
               'bias': count_table[partial_subkey]}
              for partial_subkey in range(1<<partial_subkey_nbits)]


KEY_i = byte_to_int(KEY)
target_partial_subkey = ((KEY_i >> 4)&0b11110000) + (KEY_i & 0b1111)
print(f'Target partial subkey: {target_partial_subkey:02x} ' \
      f'with bias: {27/1024:.4f}')
for prob in sorted(bias_table, key=lambda x: x['bias'], reverse=True)[:3]:
    print(f"  guess {prob['partial_subkey']:02x} " \
          f"with bias {prob['bias']/float(NPAIRS):.4f}")
