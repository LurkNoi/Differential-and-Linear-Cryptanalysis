# -*- coding: utf-8 -*-
"""
python implement of SPN ciphers
"""


def byte_to_int(b):
    return int.from_bytes(b, 'big')

def int_to_byte(n, blocksize=0):
    b = int.to_bytes(n, (n.bit_length() + 7) // 8, 'big')

    if blocksize > 0 and len(b) % blocksize:
        b = (blocksize - len(b) % blocksize) * b'\x00' + b
    return b


MODE_ENC = 0
MODE_DEC = 1

class SPN:
    """Substitution-Permutation Network (SPN) Cipher"""

    def __init__(self, sbox: list, pbox: list,
                 nrounds: int, block_size: int):
        sbox_length = len(sbox)
        pbox_length = len(pbox)
        sbox_size = (sbox_length - 1).bit_length()
        if sbox_length & (sbox_length - 1) != 0:
            raise TypeError("sbox length is not a power of 2")
        if (8 * block_size) % sbox_size != 0:
            raise TypeError("sbox size must divide block size")
        if set(sbox) != set(range(sbox_length)):
            raise TypeError("sbox is not invertible")
        if pbox_length != block_size * 8:
            raise TypeError("pbox length disagrees with block size")
        if set(pbox) != set(range(pbox_length)):
            raise TypeError("pbox is not invertible")

        self.Nr = nrounds
        self.block_size = block_size
        self.sbox_dct = {
            'sbox': sbox,
            'sbox_inv': [sbox.index(i) for i in range(sbox_length)],
            'sbox_size': sbox_size,
            'nsboxes': 8 * block_size // sbox_size,
            'mask': (1 << sbox_size) - 1
        }
        self.pbox = pbox
        self.pbox_inv = [pbox.index(i) for i in range(len(pbox))]
        self.block_mask = (1 << (8 * block_size)) - 1

    def substitution(self, block_i, mode):
        sbox_dct = self.sbox_dct
        if mode == MODE_ENC:
            sbox = sbox_dct['sbox']
        elif mode == MODE_DEC:
            sbox = sbox_dct['sbox_inv']
        sbox_size = sbox_dct['sbox_size']
        nsboxes = sbox_dct['nsboxes']
        mask = sbox_dct['mask']
        box_ = [(block_i >> (sbox_size*i))&mask for i in range(nsboxes)]
        block_i = sum([(sbox[box_[i]] << (sbox_size*i))
                       for i in range(nsboxes)])
        return block_i

    def permutation(self, block_i, mode):
        if mode == MODE_ENC:
            pbox = self.pbox
        elif mode == MODE_DEC:
            pbox = self.pbox_inv
        block_nbits = 8 * self.block_size
        bits = [(block_i >> i)&1 for i in range(block_nbits)]
        bits.reverse()
        block_i = sum([(bits[pbox[i]] << (block_nbits-1-i))
                       for i in range(block_nbits)])
        return block_i

    def key_mixing(self, block_i, subkey_i):
        block_i ^= subkey_i
        return block_i & self.block_mask

    def encrypt_block(self, block):
        raise NotImplementedError

    def decrypt_block(self, block):
        raise NotImplementedError

    def encrypt(self, data):
        block_size = self.block_size
        if len(data) % block_size != 0:
            raise TypeError('data length must be divided by block size')
        cipher = b''
        for i in range(0, len(data), block_size):
            cipher += self.encrypt_block(data[i:i+block_size])
        return cipher

    def decrypt(self, cipher):
        block_size = self.block_size
        if len(cipher) % block_size != 0:
            raise TypeError('cipher length must be divided by block size')
        data = b''
        for i in range(0, len(cipher), block_size):
            data += self.decrypt_block(cipher[i:i+block_size])
        return data


class basicSPN(SPN):
    """
    Basic Substitution-Permutation Network (SPN) Cipher as proposed in

    - Heys, Howard M. "A tutorial on linear and differential cryptanalysis."
      Cryptologia 26.3 (2002): 189-221.
    """
    SBOX = [
        0x0E, 0x04, 0x0D, 0x01,
        0x02, 0x0F, 0x0B, 0x08,
        0x03, 0x0A, 0x06, 0x0C,
        0x05, 0x09, 0x00, 0x07,
    ]
    PBOX = [
        0x00, 0x04, 0x08, 0x0C,
        0x01, 0x05, 0x09, 0x0D,
        0x02, 0x06, 0x0A, 0x0E,
        0x03, 0x07, 0x0B, 0x0F,
    ]

    def __init__(self, key, nrounds=4):
        block_size = 2
        if len(key) != block_size * (nrounds + 1):
            raise TypeError('key length must be {}'.format(
                block_size * (nrounds + 1)
            ))
        self.KEY = key
        self.subkey_ = [byte_to_int(key[i:i+block_size])
                        for i in range(0, len(key), block_size)]
        SPN.__init__(self, sbox=self.SBOX, pbox=self.PBOX,
                     nrounds=nrounds, block_size=block_size)

    def __repr__(self):
        return "basicSPN ({} rounds) with key = {}".format(self.Nr, self.KEY.hex())

    def encrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        subkey_ = self.subkey_
        for r in range(Nr - 1):
            block_i = self.key_mixing(block_i, subkey_[r])
            block_i = self.substitution(block_i, MODE_ENC)
            block_i = self.permutation(block_i, MODE_ENC)
        block_i = self.key_mixing(block_i, subkey_[Nr-1])
        block_i = self.substitution(block_i, MODE_ENC)
        block_i = self.key_mixing(block_i, subkey_[Nr])
        return int_to_byte(block_i, self.block_size)

    def decrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        subkey_ = self.subkey_
        block_i = self.key_mixing(block_i, subkey_[Nr])
        block_i = self.substitution(block_i, MODE_DEC)
        block_i = self.key_mixing(block_i, subkey_[Nr-1])
        for r in range(Nr-2, -1, -1):
            block_i = self.permutation(block_i, MODE_DEC)
            block_i = self.substitution(block_i, MODE_DEC)
            block_i = self.key_mixing(block_i, subkey_[r])
        return int_to_byte(block_i, self.block_size)


class PRESENT_80(SPN):
    """
    PRESENT-80 CIPHER
    """

    SBOX = [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD,
            0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2]
    PBOX = [
        0x00, 0x10, 0x20, 0x30, 0x01, 0x11, 0x21, 0x31,
        0x02, 0x12, 0x22, 0x32, 0x03, 0x13, 0x23, 0x33,
        0x04, 0x14, 0x24, 0x34, 0x05, 0x15, 0x25, 0x35,
        0x06, 0x16, 0x26, 0x36, 0x07, 0x17, 0x27, 0x37,
        0x08, 0x18, 0x28, 0x38, 0x09, 0x19, 0x29, 0x39,
        0x0a, 0x1a, 0x2a, 0x3a, 0x0b, 0x1b, 0x2b, 0x3b,
        0x0c, 0x1c, 0x2c, 0x3c, 0x0d, 0x1d, 0x2d, 0x3d,
        0x0e, 0x1e, 0x2e, 0x3e, 0x0f, 0x1f, 0x2f, 0x3f,
    ]

    def __init__(self, key, nrounds=31):
        block_size = 8
        if len(key) != 10:
            raise TypeError('key length must be 10')
        self.KEY = key
        self.subkey_ = self.key_schedule(key, nrounds)
        SPN.__init__(self, sbox=self.SBOX, pbox=self.PBOX,
                     nrounds=nrounds, block_size=block_size)

    def __repr__(self):
        return "PRESENT-80 ({} rounds) with key = {}".format(self.Nr, self.KEY.hex())

    def key_schedule(self, key, nrounds):
        Sbox = self.SBOX
        key_i = byte_to_int(key)
        roundkeys = []
        for i in range(1, nrounds+1):
            roundkeys.append(key_i >> 16)
            key_i = ((key_i & (2**19-1)) << 61) + (key_i >> 19)
            key_i = (Sbox[key_i >> 76] << 76) + (key_i & (2**76-1))
            key_i ^= i << 15
        return roundkeys

    def encrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        subkey_ = self.subkey_
        for r in range(Nr - 1):
            block_i = self.key_mixing(block_i, subkey_[r])
            block_i = self.substitution(block_i, MODE_ENC)
            block_i = self.permutation(block_i, MODE_ENC)
        block_i = self.key_mixing(block_i, subkey_[Nr-1])
        return int_to_byte(block_i, self.block_size)

    def decrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        subkey_ = self.subkey_
        block_i = self.key_mixing(block_i, subkey_[Nr-1])
        for r in range(Nr-2, -1, -1):
            block_i = self.permutation(block_i, MODE_DEC)
            block_i = self.substitution(block_i, MODE_DEC)
            block_i = self.key_mixing(block_i, subkey_[r])
        return int_to_byte(block_i, self.block_size)


if __name__ == '__main__':
    import os

    nrounds = 5
    key = os.urandom(2*nrounds + 2)
    cipher = basicSPN(key, nrounds=nrounds)
    print(cipher)
    plaintext = os.urandom(10)
    ciphertext = cipher.encrypt(plaintext)
    if cipher.decrypt(ciphertext) != plaintext:
        raise ValueError("sth wrong")
    print("  PASS")

    nrounds = 32
    key = os.urandom(10)
    cipher = PRESENT_80(key, nrounds=nrounds)
    print(cipher)
    plaintext = os.urandom(16)
    ciphertext = cipher.encrypt(plaintext)
    if cipher.decrypt(ciphertext) != plaintext:
        raise ValueError("sth wrong")
    print("  PASS")