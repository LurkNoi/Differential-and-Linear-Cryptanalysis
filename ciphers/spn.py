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


MODE_ENC = 0
MODE_DEC = 1

class SPN:
    """Basic Substitution-Permutation Network (SPN) Cipher"""

    def __init__(self, masterkey: bytes,
                 sbox: list,
                 pbox: list,
                 nrounds: int = 4,
                 block_size: int = 2):
        if len(masterkey) != block_size * (nrounds + 1):
            raise TypeError('key length must be {}'.format(
                block_size * (nrounds + 1)
            ))
        sbox_length = len(sbox)
        pbox_length = len(pbox)
        sbox_size = (sbox_length - 1).bit_length()
        if sbox_length & (sbox_length - 1) != 0:
            raise TypeError("sbox length is not a power of 2")
        if (8 * block_size) % sbox_size != 0:
            raise TypeError("sbox size must divide block size")
        if set(sbox) != set(range(sbox_length)):
            raise TypeError("sbox is not invertable")
        if pbox_length != block_size * 8:
            raise TypeError("pbox length disagrees with block size")
        if set(pbox) != set(range(pbox_length)):
            raise TypeError("pbox is not invertable")

        self.Nr = nrounds
        self.block_size = block_size
        self.subkey_ = [byte_to_int(masterkey[i:i+block_size])
                        for i in range(0, len(masterkey), block_size)]
        self.sbox = {
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
        sbox = self.sbox
        if mode == MODE_ENC:
            S = sbox['sbox']
        elif mode == MODE_DEC:
            S = sbox['sbox_inv']
        sbox_size = sbox['sbox_size']
        nsboxes = sbox['nsboxes']
        mask = sbox['mask']
        box_ = [(block_i >> (sbox_size*i))&mask for i in range(nsboxes)]
        block_i = sum([(S[box_[i]] << (sbox_size*i))
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

    def key_mixing(self, block_i, round_r):
        block_i ^= self.subkey_[round_r]
        return block_i & self.block_mask

    def encrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        for r in range(Nr - 1):
            block_i = self.key_mixing(block_i, r)
            block_i = self.substitution(block_i, MODE_ENC)
            block_i = self.permutation(block_i, MODE_ENC)
        block_i = self.key_mixing(block_i, Nr-1)
        block_i = self.substitution(block_i, MODE_ENC)
        block_i = self.key_mixing(block_i, Nr)
        return int_to_byte(block_i, self.block_size)

    def decrypt_block(self, block):
        block_i = byte_to_int(block) & self.block_mask
        Nr = self.Nr
        block_i = self.key_mixing(block_i, Nr)
        block_i = self.substitution(block_i, MODE_DEC)
        block_i = self.key_mixing(block_i, Nr-1)
        for r in range(Nr-2, -1, -1):
            block_i = self.permutation(block_i, MODE_DEC)
            block_i = self.substitution(block_i, MODE_DEC)
            block_i = self.key_mixing(block_i, r)
        return int_to_byte(block_i, self.block_size)

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


if __name__ == '__main__':
    import os

    nrounds = 5
    key = os.urandom(2*nrounds + 2)
    SPN_CIPHER = SPN(key, SBOX, PBOX, nrounds=nrounds)
    plaintext = os.urandom(10)
    ciphertext = SPN_CIPHER.encrypt(plaintext)
    if SPN_CIPHER.decrypt(ciphertext) != plaintext:
        raise ValueError("sth wrong")
    print("PASS")
