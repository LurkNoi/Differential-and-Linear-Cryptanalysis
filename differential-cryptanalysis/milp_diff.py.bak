# -*- coding: utf-8 -*-
"""
use MILP to find the best differential characteristic

REFERENCES:

 - Abdelkhalek, A. et al. MILP Modeling for (Large) S-boxes to Optimize
   Probability of Differential Characteristics.
   https://doi.org/10.13154/tosc.v2017.i4.99-129
"""
import os
import re
from math import gcd
from fractions import Fraction
from itertools import chain
from functools import reduce

import mip
from differential_util import difference_distribution_table, print_table
from spn import basicSPN


def lcm(lst):
    """
    return lcm of a given list
    """
    if len(lst) == 1:
        return lst[0]
    return reduce(lambda x, y: x*y//gcd(x, y), lst)

def separate_ddt(diff_dist_table):
    """
    separate DDT by diffreence
    """
    nrows = len(diff_dist_table)
    ncols = len(diff_dist_table[0])
    max_count = diff_dist_table[0][0]
    diff_dist_table[0][0] = 0
    biases = set(chain.from_iterable(diff_dist_table))
    biases.remove(0)
    biases = list(biases)
    result = []
    for bias in biases:
        truncated_table = [
            [1 if diff_dist_table[i][j] == bias else 0 for j in range(ncols)]
            for i in range(nrows)
        ]
        result.append({
            'b': bias,
            'probility': Fraction(bias, max_count),
            'weight': lcm(biases)//bias,
            'input_size': (nrows-1).bit_length(),
            'output_size': (ncols-1).bit_length(),
            'table': truncated_table,
        })
    return result

def to_product_of_sum(dct):
    """
    convert DDT zero entries to POS format
    """
    b = dct['b']
    sz_i = dct['input_size']
    sz_o = dct['output_size']
    table = dct['table']
    POS = "F{} = ".format(b)
    for x in range(1 << sz_i):
        for y in range(1 << sz_o):
            if table[x][y] == 0:
                # from lsb to msb
                X = [(x>>k)&1 for k in range(sz_i)]
                Y = [(y>>k)&1 for k in range(sz_o)]
                constraint = ""
                for k in reversed(range(sz_i)):
                    constraint += "+X{}".format(k)
                    if X[k] == 1:
                        constraint += "'"
                    constraint += "+Y{}".format(k)
                    if Y[k] == 1:
                        constraint += "'"
                POS += '(' + constraint[1:] + ')'
    return POS + ';'

def from_product_of_sum(minimized_pos):
    """
    build constraints for a given (minimized) POS
    """
    constraints = []
    for constr in re.findall(r'\([^\)]+\)', minimized_pos):
        pattern = {}
        for v in constr[1:-1].split('+'):
            if v.endswith("'"):
                pattern[v[:-1]] = 1
            else:
                pattern[v] = 0
        constraints.append(pattern.copy())
    return constraints


class MilpOptim:
    """
    MILP Optimizer for finding the best differential characteristic of SPN cipher
    """
    def __init__(self, Nr, DDTs, pbox_inv, **kwargs):
        self.Nr = Nr
        self.DDTs = DDTs
        self.pbox_inv = pbox_inv
        model = mip.Model(sense=mip.MINIMIZE, solver_name=mip.CBC)
        model.verbose = kwargs.get('verbose', 0)
        model.preprocess = kwargs.get('preprocess', 1)
        model.threads = kwargs.get('threads', 4)
        input_size = DDTs[0]['input_size']
        output_size = DDTs[0]['output_size']
        # M is a sufficiently big integer to build the conditional constraint
        self.M = kwargs.get('M', input_size+output_size+1)
        block_size = len(pbox_inv)
        nsbox = block_size // input_size
        x_ = [model.add_var(name=f'x_{i}', var_type=mip.BINARY)
              for i in range(block_size*Nr)]
        Q_ = [model.add_var(name=f'Q_{t}', var_type=mip.BINARY)
              for t in range(nsbox*Nr)]
        for dct in DDTs:
            dct['Qb'] = [model.add_var(name=f"Q{dct['b']}_{t}", var_type=mip.BINARY)
                         for t in range(nsbox*Nr)]
        for t in range(nsbox*Nr):
            model += mip.xsum(
                dct['Qb'][t] for dct in DDTs
            ) == Q_[t]
        objective = 0
        for dct in DDTs:
            objective += dct['weight'] * mip.xsum(dct['Qb'][:-nsbox])
        model.objective = objective
        # avoid trivial solution
        model += mip.xsum(Q_) >= 1
        self.model = model
        self.x = x_
        self.Q = Q_
        self.block_size = block_size
        self.nsbox = nsbox

    def add_constr(self, ddt_dct, X, Y, t):
        M = self.M
        m = self.model
        sz_i = ddt_dct['input_size']
        sz_o = ddt_dct['output_size']
        for pattern in ddt_dct['CONSTR']:
            m += mip.xsum(
                X[i]*(1-pattern.get(f"X{sz_i-1-i}", 1))
                + (1-X[i])*pattern.get(f"X{sz_i-1-i}", 0)
                for i in range(sz_i)
            ) + mip.xsum(
                Y[j]*(1-pattern.get(f"Y{sz_o-1-j}", 1))
                + (1-Y[j])*pattern.get(f"Y{sz_o-1-j}", 0)
                for j in range(sz_o)
            ) + M*(1 - ddt_dct['Qb'][t]) >= 1

    def add_sbox(self, xs, t, **kwargs):
        ys = kwargs.get('ys', None)
        m = self.model
        Q_ = self.Q
        sum_xs = mip.xsum(xs)
        m += sum_xs >= Q_[t]
        for xs_i in xs:
            m += Q_[t] >= xs_i
        if ys is None:
            return
        sum_ys = mip.xsum(ys)
        m += len(ys) * sum_xs >= sum_ys
        m += len(xs) * sum_ys >= sum_xs
        for dct in self.DDTs:
            self.add_constr(dct, xs, ys, t)

    def get_prob(self):
        prob = Fraction(1, 1)
        for dct in self.DDTs:
            for t in range((self.Nr-1)*self.nsbox):
                if dct['Qb'][t].x == 1:
                    prob *= dct['probility']
        return prob

    def print_path(self):
        bs = self.block_size
        x_ = self.x
        Nr = self.Nr
        for r in range(Nr):
            for i in range(bs):
                if x_[bs*r+i].x == 1:
                    print("1", end=' ')
                else:
                    print(".", end=' ')
            print()

    def solve(self, **kwargs):
        max_seconds = kwargs.get('max_seconds', 30)
        max_solutions = kwargs.get('max_solutions', 1)
        write_file = kwargs.get('write_file', None)
        sol_x_ = None
        num_sol = 0
        PBOX_INV = self.pbox_inv
        m = self.model
        Nr = self.Nr
        x_ = self.x
        bs = self.block_size
        nsbox = self.nsbox
        sz = bs // nsbox
        self.solutions = []
        for r in range(Nr):
            block_in = x_[bs*r : bs*(r+1)]
            if r == Nr-1:
                for i in range(nsbox):
                    t = nsbox*r + i
                    self.add_sbox(xs=block_in[sz*i : sz*(i+1)], t=t)
                break
            next_block_in = x_[bs*(r+1) : bs*(r+2)]
            block_out = [next_block_in[PBOX_INV[i]] for i in range(bs)]
            for i in range(4):
                t = nsbox*r + i
                self.add_sbox(xs=block_in[sz*i : sz*(i+1)], t=t,
                              ys=block_out[sz*i : sz*(i+1)])

        while num_sol < max_solutions:
            if sol_x_ is not None:
                # remove previous solution
                m += mip.xsum(
                    sol_x_[i]*(1-x_[i]) + (1-sol_x_[i])*x_[i]
                    for i in range(bs*Nr)
                ) >= 1
            if write_file is not None:
                m.write(write_file + '.lp')
            try:
                status = m.optimize(max_seconds=max_seconds)
                if status == mip.OptimizationStatus.OPTIMAL:
                    print('optimal solution cost {} found'.format(
                        m.objective_value))
                elif status == mip.OptimizationStatus.FEASIBLE:
                    print('sol.cost {} found, best possible: {}'.format(
                        m.objective_value, m.objective_bound))
                elif status == mip.OptimizationStatus.NO_SOLUTION_FOUND:
                    raise ValueError('no feasible solution found, '\
                                     'lower bound is: {}'.format(m.objective_bound))
            except Exception as e:
                print(e)
                break

            sol_x_ = [int(x_[i].x) for i in range(bs*Nr)]

            prob = self.get_prob()
            if prob == 0:
                raise ValueError("probility shouldn't be zero")
            num_sol += 1
            self.print_path()
            print('diff. prob. = {}'.format(prob))
            # N_D \approx c/p_D
            print('number of paris N_D = {}c'.format(int(1/prob)))





if __name__ == '__main__':
    NR = 4
    key = os.urandom(2*NR + 2)
    cipher = basicSPN(key, nrounds=NR)
    SBOX = cipher.sbox_dct['sbox']
    PBOX_INV = cipher.pbox_inv
    DDT = difference_distribution_table(SBOX)

    print_table(DDT, nonzero=True)
    # for dct in separate_ddt(DDT):
    #     print(to_product_of_sum(dct))
    #     print_table(dct['table'], nonzero=True)
    # exit()

    DDTS = separate_ddt(DDT)
    for dct in DDTS:
        if dct['b'] == 2:
            dct['MIN_POS'] = "F2 = (X3+Y2+X0+Y0)(X2+Y2+X1+X0)(X3+X2+X1+Y0)(Y2+X1'+Y1+Y0)(X3'+Y2'+X1+X0'+Y0')(X3+Y3+Y1+X0)(X3+Y2'+X1'+X0'+Y0')(X3'+Y2'+X1'+X0+Y0')(X3'+Y3'+Y1+X0)(X3+Y2'+X1+Y1'+X0)(X2+Y2'+Y1'+X0'+Y0)(X3+X2'+Y2+X1+X0')(X3'+Y3'+Y2+Y1+Y0')(X2'+Y2+X1'+X0+Y0')(X3'+X2+Y2'+X1'+Y0)(X3'+Y3+Y2+Y1'+Y0')(X3+Y3'+Y2+Y1'+Y0')(X2+Y2+X1'+X0'+Y0')(X3'+Y3+Y2+X0'+Y0)(X2+Y1+X0+Y0)(X2'+Y2'+X1'+X0'+Y0')(X3'+Y3+X2'+X0'+Y0)(Y3+X2'+Y2'+Y1'+Y0')(X2'+Y2+X1+X0'+Y0')(X3+Y3'+X2+Y2'+X0+Y0')(Y3+X2+Y2'+X1+Y1+Y0')(Y3+X2'+Y2'+X1'+Y1'+X0)(X3+Y3+Y2+Y1)(Y2+X1+X0+Y0)(X3'+Y3'+Y2+X1'+Y0)(Y3'+X2'+Y2'+X1+Y1+Y0')(X3+Y3'+Y2'+Y1'+X0'+Y0)(Y3'+X2'+X1'+Y1+X0'+Y0)(Y3'+X2'+Y2'+X1+Y1'+X0)(X3+Y3'+X2+X1+Y1')(X3+Y3+X1+Y1+Y0);"
        elif dct['b'] == 4:
            dct['MIN_POS'] = "F4 = (X3+X2+X0)(X3'+X1+Y0)(X1'+X0'+Y0)(X2'+X0+Y0)(Y2'+Y1+X0'+Y0')(X3'+Y2'+X1'+Y0')(X3+Y3+Y2'+X0')(X3+Y3'+X2'+Y1)(Y3+Y2+Y1'+X0')(Y3+X2+X0)(Y3'+Y2+X1'+Y1')(Y3'+Y2'+X1+Y1')(Y3+Y2+X1'+Y1)(Y3'+X1+X0'+Y0')(X3'+X2'+X1+X0)(X3'+X2+Y2+Y1)(Y3'+X1'+Y1+X0)(Y3+X1+X0)(Y2+Y1+Y0)(Y2'+Y1'+X0+Y0')(Y3+X2+Y1);"
        elif dct['b'] == 6:
            dct['MIN_POS'] = "F6 = (X0'+Y0')(X1'+Y1')(Y2+X1)(X1+X0')(Y3'+Y2'+Y0)(X3'+Y3+X0)(X3+X2'+Y1)(X2+Y2'+Y0)(Y3'+X1'+Y0')(Y1'+Y0')(Y2+X0')(X2+X1)(X3+Y2);"
        elif dct['b'] == 8:
            dct['MIN_POS'] = "F8 = (Y0')(Y2')(X2')(Y3')(X3)(X1)(Y1)(X0);"
        dct['CONSTR'] = from_product_of_sum(dct['MIN_POS'])
        # print(dct)

    MILP_OPTM = MilpOptim(Nr=NR, DDTs=DDTS, pbox_inv=PBOX_INV)
    MILP_OPTM.solve(max_solutions=1)

    """
    best differential characteristic for basicSPN (4 round)

    optimal solution cost 14.0 found
    1 . 1 1 . . . . 1 . 1 1 . . . .
    . . . . . . . . 1 . 1 . . . . .
    . . 1 . . . . . . . . . . . . .
    . . . . 1 . . . . . . . 1 . . .
    diff. prob. = 9/256
    number of paris N_D = 28c
    """
