# -*- coding: utf-8 -*-
"""
Use milp_util.py to find the best linear trails
of round reduced PRESENT-80 cipher
"""
import os
from math import gcd, log2
from fractions import Fraction
from itertools import chain
from functools import reduce

import mip
from milp_util import *


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
PBOX_INV = [PBOX.index(i) for i in range(len(PBOX))]


def lcm(lst):
    """
    return lcm of a given list
    """
    if len(lst) == 1:
        return lst[0]
    return reduce(lambda x, y: x*y//gcd(x, y), lst)


def separate_lat(lin_appr_table):
    """
    separate LAT by linear bias
    """
    nrows = len(lin_appr_table)
    ncols = len(lin_appr_table[0])
    lin_appr_table[0][0] = 0
    biases = set(chain.from_iterable(lin_appr_table))
    biases.remove(0)
    biases = list(biases)
    result = []
    for bias in biases:
        truncated_table = [
            [1 if lin_appr_table[i][j] == bias else 0 for j in range(ncols)]
            for i in range(nrows)
        ]
        to_pla(truncated_table, f'./result/lat{bias}.pla')
        simplify(f'./result/lat{bias}.pla', f'./result/lat{bias}_sp.pla')
        DATA_SP = from_pla(f'./result/lat{bias}_sp.pla')
        print(f"bias({bias}) requires {len(DATA_SP)} constraints")
        CONSTR = parse_constr(DATA_SP)
        result.append({
            'b': bias,
            'bias': Fraction(bias, nrows),
            'weight': (lcm(biases)//bias)**2,
            'constr': CONSTR,
        })
    return result


def add_sbox(model, Q_t, xs, ys=None):
    # non-zero input `xs` active the corresponding sbox
    model += mip.xsum(xs) >= Q_t
    for xs_i in xs:
        model += Q_t >= xs_i
    if ys is not None:
        # Because impossible constraints of SBox is conditional,
        # we need additional constraint for non-zero output `ys`
        # resulting in non-zero input `xs`
        model += len(ys) * mip.xsum(xs) >= mip.xsum(ys)


def add_cond_impos(model, X, Y, Qb_t):
    """
    add conditional constraints
    """
    sz_i = sz_o = 4
    # M is a sufficiently big integer to build the conditional constraint
    M = sz_i + sz_o + 1
    model += mip.xsum(
        X[i] * (1-pattern.get(f"X{i}", 1))
        + (1-X[i]) * pattern.get(f"X{i}", 0)
        for i in range(sz_i)
    ) + mip.xsum(
        Y[j] * (1-pattern.get(f"Y{j}", 1))
        + (1-Y[j]) * pattern.get(f"Y{j}", 0)
        for j in range(sz_o)
    ) + M*(1 - Qb_t) >= 1


def print_path(Nr, Q, U):
    global PBOX_INV
    bs = 64
    nsbox = 16
    for r in range(Nr):
        for i in range(bs):
            if int(U[r][i].x) == 1:
                print("1", end='')
            else:
                print(".", end='')
        print()
        for i in range(nsbox):
            if int(Q[r][i].x) == 1:
                print("|SB|", end='')
            else:
                print('-'*(bs//nsbox), end='')
        print()
        if r != Nr - 1:
            for i in range(bs):
                if int(U[r+1][PBOX_INV[i]].x) == 1:
                    print("1", end='')
                else:
                    print(".", end='')
            print()
        print()



LAT = linear_approximation_table(SBOX, absolute=True)
LATS = separate_lat(LAT)

max_seconds = 3600
Nr = 5
sz = 4 # 4-bit S-Box
bs = 64 # block size
nsbox = bs // sz # number of sbox at each round

SLVR = MilpOptim(sense=mip.MINIMIZE, solver_name=mip.CBC)
SLVR.verbose = 0
SLVR.preprocess = 1
SLVR.threads = 6

U_ = [SLVR.new_var("U{}".format(r), bs) for r in range(Nr)]
Q_ = [SLVR.new_var("Q{}".format(r), nsbox) for r in range(Nr)]
for table in LATS:
    table['Qb'] = [SLVR.new_var("Q{}_{}".format(table['b'], r), nsbox)
                   for r in range(Nr)]
for r in range(Nr):
    for i in range(nsbox):
        SLVR += mip.xsum(
            table['Qb'][r][i] for table in LATS
        ) == Q_[r][i]

SLVR.objective = mip.xsum(
    table['weight'] * mip.xsum(chain.from_iterable(table['Qb'][:-1]))
    for table in LATS
)

SLVR += mip.xsum(Q_[-1]) >= 1 # avoid trivial solution
# SLVR += mip.xsum(Q_[-1]) <= 4 # brute force at most 2 bytes at a time


r = Nr - 1
Ur = U_[r]
for i in range(nsbox):
    add_sbox(model=SLVR, Q_t=Q_[r][i], xs=Ur[sz*i : sz*(i+1)])

for table in LATS:
    for r in range(Nr - 1):
        Ur = U_[r]
        Ur_next = U_[r+1] # U_{r+1}
        Vr = [Ur_next[PBOX_INV[i]] for i in range(bs)]
        for i in range(nsbox):
            X = Ur[sz*i : sz*(i+1)]
            Y = Vr[sz*i : sz*(i+1)]
            add_sbox(model=SLVR, Q_t=Q_[r][i], xs=X, ys=Y)
            Qt = table['Qb'][r][i]
            for pattern in table['constr']:
                add_cond_impos(model=SLVR, X=X, Y=Y, Qb_t=Qt)

# SLVR.write('present.lp')
print(SLVR)
SLVR.check(max_seconds=max_seconds)

print_path(Nr, Q_, U_)

# get bias
bias = Fraction(1, 1)
for table in LATS:
    for r in range(Nr - 1):
        for i in range(nsbox):
            if int(table['Qb'][r][i].x) == 1:
                # Piling-Up Lemma
                bias *= 2*table['bias']
bias /= 2
print(f"linear bias = {bias} ({log2(bias):.6f})")
# N_L \approx 1/(e^2)
print(f"number of paris N_L = {int(1/(bias**2))}\n")

"""
result:

Using Python-MIP package version 1.8.1
bias(4) requires 37 constraints
bias(2) requires 12 constraints
MIP (CBC) model with 4065 constraints
Coin3009W Conflict graph built in 0.001 seconds, density: 0.232%
optimal solution cost 6.0 found
............1..1............................................1111
------------|SB|--------------------------------------------|SB|
...............1...............................................1

............................................................1..1
------------------------------------------------------------|SB|
...............................................................1

...............................................................1
------------------------------------------------------------|SB|
.............................................................1.1

.......................................................1.......1
----------------------------------------------------|SB|----|SB|
.....................................................1.1.....111

.......................1.......1.......................1...1...1
--------------------|SB|----|SB|--------------------|SB||SB||SB|

linear bias = 1/128 (-7.000000)
number of paris N_L = 16384


real    3m13.098s
user    3m11.094s
sys     0m1.969s
"""
