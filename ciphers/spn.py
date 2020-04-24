# -*- coding: utf-8 -*-
"""
A python implement of basic SPN cipher as proposed in

- Heys, Howard M. "A tutorial on linear and differential cryptanalysis."
  Cryptologia 26.3 (2002): 189-221.
"""

SBOX = [
    0x0E, 0x04, 0x0D, 0x01,
    0x02, 0x0F, 0x0B, 0x08,
    0x03, 0x0A, 0x06, 0x0C,
    0x05, 0x09, 0x00, 0x07,
]

SBOX_INV = [
    0x0E, 0x03, 0x04, 0x08,
    0x01, 0x0C, 0x0A, 0x0F,
    0x07, 0x0D, 0x09, 0x06,
    0x0B, 0x02, 0x00, 0x05,
]

PBOX = [
    0x00, 0x04, 0x08, 0x0C,
    0x01, 0x05, 0x09, 0x0D,
    0x02, 0x06, 0x0A, 0x0E,
    0x03, 0x07, 0x0B, 0x0F,
]

# since P-box is symmetric
PBOX_INV = PBOX

def byte_to_int(b):
    return int.from_bytes(b, 'big')

def int_to_byte(n, blocksize=0):
    b = int.to_bytes(n, (n.bit_length() + 7) // 8, 'big')

    if blocksize > 0 and len(b) % blocksize:
        b = (blocksize - len(b) % blocksize) * b'\x00' + b
    return b


class SPN:
    """Basic Substitution-Permutation Network (SPN) Cipher"""

    def __init__(self, masterkey: bytes):
        if len(masterkey) != 10:
            raise TypeError('key length must be 10')
        self.subkey_ = [byte_to_int(masterkey[i:i+2])
                        for i in range(0, 10, 2)]

    def substitution(self, block: int, SBox: list):
        box_ = [(block >> (4*i))&0x0F for i in range(4)]
        block = sum([(SBox[box_[i]] << (4*i))
                     for i in range(4)])
        return block

    def permutation(self, block: int, PBox: list):
        bits = [(block >> i)&1 for i in range(16)]
        bits.reverse()
        block = sum([(bits[PBox[i]] << (15-i))
                     for i in range(16)])
        return block

    def key_mixing(self, block: int, r: int):
        block ^= self.subkey_[r]
        return block & 0xFFFF

    def encrypt_block(self, block: bytes):
        block = byte_to_int(block) & 0xFFFF
        for r in range(3):
            block = self.key_mixing(block, r)
            block = self.substitution(block, SBOX)
            block = self.permutation(block, PBOX)
        block = self.key_mixing(block, 3)
        block = self.substitution(block, SBOX)
        block = self.key_mixing(block, 4)
        return int_to_byte(block, 2)

    def decrypt_block(self, block: bytes):
        block = byte_to_int(block) & 0xFFFF
        block = self.key_mixing(block, 4)
        block = self.substitution(block, SBOX_INV)
        block = self.key_mixing(block, 3)
        for r in range(2, -1, -1):
            block = self.permutation(block, PBOX_INV)
            block = self.substitution(block, SBOX_INV)
            block = self.key_mixing(block, r)
        return int_to_byte(block, 2)

    def encrypt(self, data: bytes):
        if len(data) % 2 != 0:
            raise TypeError('data length must be even')
        cipher = b''
        for i in range(0, len(data), 2):
            cipher += self.encrypt_block(data[i:i+2])
        return cipher

    def decrypt(self, cipher: bytes):
        if len(cipher) % 2 != 0:
            raise TypeError('cipher length must be even')
        data = b''
        for i in range(0, len(cipher), 2):
            data += self.decrypt_block(cipher[i:i+2])
        return data


if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('key')
    parser.add_argument('input_file')
    parser.add_argument('--dec', action='store_true')
    parser.add_argument('-o', '--out', default=None)

    args = parser.parse_args()
    key = bytes.fromhex(args.key)
    is_dec = args.dec
    data = open(args.input_file, 'rb').read()
    if args.out is None:
        write_output = False
    else:
        write_output = True
        output_file = args.out

    spn = SPN(key)
    if is_dec:
        output = spn.decrypt(data)
    else:
        output = spn.encrypt(data)
    if not write_output:
        print('output: {}'.format(output.hex()))
    else:
        open(output_file, 'wb').write(output)
