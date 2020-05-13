# -*- coding: utf-8 -*-
"""
example using milp_util.py to solve zer0SPN
"""
from math import log2
from fractions import Fraction
from itertools import chain

import mip
from milp_util import *


def linear_approximation_table(sbox, absolute=False, lbound=0):
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
        absolute_bias_table[input_mask][output_mask] = count if count >= lbound else 0

    return absolute_bias_table


def add_sbox(model, Q_t, xs):
    # non-zero input `xs` active the corresponding sbox
    model += mip.xsum(xs) >= Q_t
    for xs_i in xs:
        model += Q_t >= xs_i


def get_prob(Nr, sol, LAT, ptable_inv):
    """
    caculate the prob.
    """
    nsbox = 8
    sz = 8
    prob = Fraction(1, 1)
    for r in range(Nr-1):
        for t in range(nsbox):
            if prob == 0:
                raise ValueError("prob. shouldn't be zero")
            if int(sol[f"Q{r}_{t}"]) == 1:
                in_bits = [int(sol[f"U{r}_{k}"]) for k in range(sz*t, sz*(t+1))]
                in_mask = sum(in_bits[k] << (sz-1-k) for k in range(sz))
                out_bits = [int(sol[f"U{r+1}_{ptable_inv[k]}"]) for k in range(sz*t, sz*(t+1))]
                out_mask = sum(out_bits[k] << (sz-1-k) for k in range(sz))
                # print(r, t, in_mask, out_mask, LAT[in_mask][out_mask])
                prob *= LAT[in_mask][out_mask]
                prob /= 1 << sz
                prob *= 2
    return prob/2


def print_path(Nr, sol):
    bs = 64
    nsbox = 8
    for r in range(Nr):
        for i in range(bs):
            if int(sol[f'U{r}_{i}']) == 1:
                print("1", end='')
            else:
                print(".", end='')
        print()
        for i in range(nsbox):
            if int(sol[f'Q{r}_{i}']) == 1:
                print("| SBox |", end='')
            else:
                print('-'*(bs//nsbox), end='')
        print()
        if r != Nr - 1:
            for i in range(bs):
                if int(sol[f'U{r+1}_{ptable[i]}']) == 1:
                    print("1", end='')
                else:
                    print(".", end='')
            print()
        print()



sbox = [62, 117, 195, 179, 20, 210, 41, 66, 116, 178, 152, 143, 75, 105, 254, 1, 158, 95, 101, 175, 191, 166, 36, 24, 50, 39, 190, 120, 52, 242, 182, 185, 61, 225, 140, 38, 150, 80, 19, 109, 246, 252, 40, 13, 65, 236, 124, 186, 214, 86, 235, 100, 97, 49, 197, 154, 176, 199, 253, 69, 88, 112, 139, 77, 184, 45, 133, 104, 15, 54, 177, 244, 160, 169, 82, 148, 73, 30, 229, 35, 79, 137, 157, 180, 248, 163, 241, 231, 81, 94, 165, 9, 162, 233, 18, 85, 217, 84, 7, 55, 63, 171, 56, 118, 237, 132, 136, 22, 90, 221, 103, 161, 205, 11, 255, 14, 122, 47, 71, 201, 99, 220, 83, 74, 173, 76, 144, 16, 155, 126, 60, 96, 44, 234, 17, 215, 107, 138, 159, 183, 251, 3, 198, 0, 89, 170, 131, 151, 219, 29, 230, 32, 187, 125, 134, 64, 12, 202, 164, 247, 25, 223, 222, 119, 174, 67, 147, 146, 206, 51, 243, 53, 121, 239, 68, 130, 70, 203, 211, 111, 108, 113, 8, 106, 57, 240, 21, 93, 142, 238, 167, 5, 128, 72, 189, 192, 193, 92, 10, 204, 87, 145, 188, 172, 224, 226, 207, 27, 218, 48, 33, 28, 123, 6, 37, 59, 4, 102, 114, 91, 23, 209, 34, 42, 2, 196, 141, 208, 181, 245, 43, 78, 213, 216, 232, 46, 98, 26, 212, 58, 115, 194, 200, 129, 227, 249, 127, 149, 135, 228, 31, 153, 250, 156, 168, 110]
ptable = [
    0, 8, 16, 24, 32, 40, 48, 56,
    1, 9, 17, 25, 33, 41, 49, 57,
    2, 10, 18, 26, 34, 42, 50, 58,
    3, 11, 19, 27, 35, 43, 51, 59,
    4, 12, 20, 28, 36, 44, 52, 60,
    5, 13, 21, 29, 37, 45, 53, 61,
    6, 14, 22, 30, 38, 46, 54, 62,
    7, 15, 23, 31, 39, 47, 55, 63
]
ptable_inv = [ptable.index(i) for i in range(len(ptable))]



# ignore small (< 40) bias
LAT = linear_approximation_table(sbox, absolute=True, lbound=40)

to_pla(LAT, './result/lat.pla')
simplify('./result/lat.pla', './result/lat_sp.pla')
DATA_SP = from_pla('./result/lat_sp.pla')
print(f"each sbox requires {len(DATA_SP)} constraints")
CONSTR = parse_constr(DATA_SP)
# exit()


max_solutions = 3
max_seconds = 3600
Nr = 4
sz = 8 # 8-bit S-Box
bs = 64 # block size
nsbox = bs // sz # number of sbox at each round


SLVR = MilpOptim(sense=mip.MINIMIZE, solver_name=mip.CBC)
SLVR.verbose = 0
SLVR.preprocess = 1
SLVR.threads = 6

x_ = [SLVR.new_var("U{}".format(r), bs) for r in range(Nr)]
Q_ = [SLVR.new_var("Q{}".format(r), nsbox) for r in range(Nr)]

SLVR += mip.xsum(Q_[-1]) <= 2 # brute force at most 2 bytes at a time
SLVR += mip.xsum(Q_[-1]) >= 1 # remove the trivial solution
# ignore the last round
SLVR.objective = mip.xsum(chain.from_iterable(Q_[:-1]))


for r in range(Nr):
    Ur = x_[r]
    if r == Nr - 1:
        for i in range(nsbox):
            add_sbox(model=SLVR, Q_t=Q_[r][i], xs=Ur[sz*i : sz*(i+1)])
        break
    Ur_next = x_[r+1] # U_{r+1}
    Vr = [Ur_next[ptable_inv[i]] for i in range(bs)]
    for i in range(nsbox):
        X = Ur[sz*i : sz*(i+1)]
        Y = Vr[sz*i : sz*(i+1)]
        add_sbox(model=SLVR, Q_t=Q_[r][i], xs=X)
        for pattern in CONSTR:
            SLVR.add_impos(pattern=pattern, X=X, Y=Y)


sols = SLVR.solve(num_solutions=max_solutions, max_seconds=max_seconds)

for sol in sols:
    print_path(Nr, sol)
    prob = get_prob(Nr, sol, LAT, ptable_inv)
    print(f"linear bias = {prob} ({log2(prob):.6f})")
    # N_L \approx 1/(e^2)
    print(f"number of paris N_L = {int(1/(prob**2))}\n")


"""
result:

...........11..1................................................
--------| SBox |------------------------------------------------
.........1......................................................

.........1......................................................
--------| SBox |------------------------------------------------
........1.....1.................................................

.1...............................................1..............
| SBox |----------------------------------------| SBox |--------
1.....1.........................................1.....1.........

1.....1.........................................1.....1.........
| SBox |----------------------------------------| SBox |--------

linear bias = 1500625/33554432 (-4.482868)
number of paris N_L = 499

..1111.1........................................................
| SBox |--------------------------------------------------------
1...............................................................

1...............................................................
| SBox |--------------------------------------------------------
...1...1........................................................

........................1...............................1.......
------------------------| SBox |------------------------| SBox |
...........................1...1...........................1...1

...........................1...1...........................1...1
------------------------| SBox |------------------------| SBox |

linear bias = 1029/65536 (-5.992973)
number of paris N_L = 4056

...11..1........................................................
| SBox |--------------------------------------------------------
.1..............................................................

........1.......................................................
--------| SBox |------------------------------------------------
...........1...1................................................

.........................1...............................1......
------------------------| SBox |------------------------| SBox |
........................1.....1.........................1.....1.

...1...1...........................................1...1........
| SBox |----------------------------------------| SBox |--------

linear bias = 300125/8388608 (-4.804796)
number of paris N_L = 781
"""
