# -*- coding: utf-8 -*-
"""
try to use MILP finding smallest number of active S-boxes

REFERENCES:

 - Sun S. et al. Automatic Security Evaluation of Block Ciphers with
   S-bP Structures Against Related-Key Differential Attacks.
   https://doi.org/10.1007/978-3-319-12087-4_3
"""
import os
from fractions import Fraction

import mip
from differential_util import difference_distribution_table
from spn import basicSPN


def add_xor(model, a, b, c):
    r"""
    add constraints for xor: ``a`` \oplus ``b`` = ``c``
    by removing four impossible case
    """
    model += a + b + (1-c) >= 1
    model += a + (1-b) + c >= 1
    model += (1-a) + b + c >= 1
    model += (1-a) + (1-b) + (1-c) >= 1

def hamming_weight(num):
    cnt = 0
    while num:
        cnt += 1
        num &= num - 1
    return cnt

def branch_number(DDT):
    """
    caculate the branch number from DDT of sbox
    """
    input_size = (len(DDT)-1).bit_length()
    output_size = (len(DDT[0])-1).bit_length()
    min_weight = input_size + output_size
    for i in range(1 << input_size):
        for j in range(1 << output_size):
            if min_weight == 2:
                return 2
            if (DDT[i][j] != 0) and (i != j):
                wt = hamming_weight((i << output_size) + j)
                if wt < min_weight:
                    min_weight = wt
    return min_weight

def add_sbox(model, xs, A_t, ys=None, B_S=None):
    """
    add constraints for sbox

    INPUT:

    - ``model`` - a mip Model instance
    - ``B_s`` - branch number of S-Box
    - ``xs`` - list of S-Box input bits, None for last round
    - ``ys`` - list of S-Box output bits, None for last round
    """
    input_size = len(xs)
    sum_xs = mip.xsum(xs[i] for i in range(input_size))
    model += sum_xs >= A_t
    for i in range(input_size):
        model += A_t >= xs[i]
    if (ys is None) or (B_S is None):
        return
    output_size = len(ys)
    sum_ys = mip.xsum(ys[j] for j in range(output_size))
    sum_xs_ys = sum_xs + sum_ys
    for i in range(input_size):
        model += sum_xs_ys >= B_S * xs[i]
    for j in range(output_size):
        model += sum_xs_ys >= B_S * ys[j]
    model += input_size * sum_ys >= sum_xs
    model += output_size * sum_xs >= sum_ys


def print_path(Nr, sol_x_, sol_y_, sol_A_):
    """
    print the path for basicSPN
    """
    for r in range(Nr):
        for i in range(16*r, 16*(r+1)):
            if i%4 == 0:
                print(' ', end='')
            if sol_x_[i] == 1:
                print('1', end=' ')
            else:
                print('.', end=' ')
        print()
        for i in range(4):
            t = 4*r + i
            if sol_A_[t] == 1:
                print('| S-BOX |', end='')
            else:
                print(' '*9, end='')
        print()
        if r == Nr - 1:
            break
        for i in range(16*r, 16*(r+1)):
            if i%4 == 0:
                print(' ', end='')
            if sol_y_[i] == 1:
                print('1', end=' ')
            else:
                print('.', end=' ')
        print()
        print('X'*36)

def get_prob(sol_x_, sol_A_, DDT):
    """
    caculate the prob.
    """
    prob = Fraction(1, 1)
    for t in range(len(sol_A_) - 4):
        if sol_A_[t] == 1:
            in_diff = int(8*sol_x_[4*t] + 4*sol_x_[4*t+1]
                          + 2*sol_x_[4*t+2] + sol_x_[4*t+3])
            out_diff = int(8*sol_y_[4*t] + 4*sol_y_[4*t+1]
                           + 2*sol_y_[4*t+2] + sol_y_[4*t+3])
            prob *= DDT[in_diff][out_diff]
            prob /= 16
    return prob

# Calculating the Minimum Number of Active S-boxes
NR = 4
key = os.urandom(2*NR + 2)
cipher = basicSPN(key, nrounds=NR)
print(cipher)
SBOX = cipher.sbox_dct['sbox']
PBOX_INV = cipher.pbox_inv
DDT = difference_distribution_table(SBOX, truncated=True)
B_S = branch_number(DDT)
number_of_sbox = 4 * NR
number_of_vars = 16 * NR
max_solutions = 1

m = mip.Model(sense=mip.MINIMIZE, solver_name=mip.CBC)
m.verbose = 0
m.preprocess = 1
m.threads = 4

x_ = [m.add_var(name=f'x{i}', var_type=mip.BINARY)
      for i in range(number_of_vars)]
A_ = [m.add_var(name=f'S{t}', var_type=mip.BINARY)
      for t in range(number_of_sbox)]
sum_At = mip.xsum(A_[t] for t in range(number_of_sbox))
m.objective = sum_At
m += sum_At >= 1

# add constraints
for r in range(NR):
    block_in = x_[16*r : 16*(r+1)]
    if r == NR - 1:
        for i in range(4):
            t = 4*r + i
            add_sbox(model=m, xs=block_in[4*i : 4*(i+1)], A_t=A_[t])
        break
    next_block_in = x_[16*(r+1) : 16*(r+2)]
    block_out = [next_block_in[PBOX_INV[i]] for i in range(16)]
    for i in range(4):
        t = 4*r + i
        add_sbox(model=m, xs=block_in[4*i : 4*(i+1)], A_t=A_[t],
                 ys=block_out[4*i : 4*(i+1)], B_S=B_S)

# m.write('spn.lp')
sol_x_ = None
num_sol = 0
while num_sol < max_solutions:
    if sol_x_ is not None:
        # remove previous solution
        m += mip.xsum(
            sol_x_[i]*(1-x_[i]) + (1-sol_x_[i])*x_[i]
            for i in range(number_of_vars)
        ) >= 1
    status = m.optimize(max_seconds=30)
    if status == mip.OptimizationStatus.OPTIMAL:
        print('optimal solution cost {} found'.format(
            m.objective_value))
    elif status == mip.OptimizationStatus.FEASIBLE:
        print('sol.cost {} found, best possible: {}'.format(
            m.objective_value, m.objective_bound))
    elif status == mip.OptimizationStatus.NO_SOLUTION_FOUND:
        raise ValueError('no feasible solution found, '\
                         'lower bound is: {}'.format(m.objective_bound))

    sol_x_ = [int(x_[i].x) for i in range(number_of_vars)]
    sol_A_ = [int(A_[t].x) for t in range(number_of_sbox)]
    sol_y_ = sol_x_[16:]
    for r in range(NR - 1):
        next_block_in = sol_y_[16*r : 16*(r+1)]
        sol_y_[16*r : 16*(r+1)] = [next_block_in[PBOX_INV[i]]
                                   for i in range(16)]

    prob = get_prob(sol_x_, sol_A_, DDT)
    if prob == 0:
        continue
    num_sol += 1
    print_path(NR, sol_x_, sol_y_, sol_A_)
    print('diff. prob. = {}\n'.format(prob))

"""result-1 -- not feasible
optimal solution cost 6.0 found
 . . . .  1 1 1 1  . . . .  . . . .
         | S-BOX |
 . . . .  1 . . 1  . . . .  . . . .
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 . 1 . .  . . . .  . . . .  . 1 . .
| S-BOX |                  | S-BOX |
 1 1 . .  . . . .  . . . .  1 1 . .
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 1 . . 1  1 . . 1  . . . .  . . . .
| S-BOX || S-BOX |
 . 1 . .  . 1 . .  . . . .  . . . .
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 . . . .  1 1 . .  . . . .  . . . .
         | S-BOX |
diff. prob. = 1/1048576
"""
