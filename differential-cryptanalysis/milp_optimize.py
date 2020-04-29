# -*- coding: utf-8 -*-
"""
try to use MILP to Optimize Probability of Differential Characteristics

REFERENCES:

 - Abdelkhalek, A. et al. MILP Modeling for (Large) S-boxes to Optimize
   Probability of Differential Characteristics.
   https://doi.org/10.13154/tosc.v2017.i4.99-129
"""
import os
import re
from fractions import Fraction

import mip
from differential_util import difference_distribution_table, print_table
from spn import basicSPN


class MilpOptim:
    """
    MILP Optimizer for S-boxes
    """
    def __init__(self, Nr, table, pbox_inv, **kwargs):
        input_length = len(table)
        output_length = len(table[0])
        self.input_size = (input_length - 1).bit_length()
        self.output_size = (output_length - 1).bit_length()
        self.table = table
        self.Nr = Nr
        self.pbox_inv = pbox_inv
        model = mip.Model(sense=mip.MINIMIZE, solver_name=mip.CBC)
        model.verbose = kwargs.get('verbose', 0)
        model.preprocess = kwargs.get('preprocess', 1)
        model.threads = kwargs.get('threads', 4)
        ignore_last_round = kwargs.get('ignore_last_round', True)
        x_ = [model.add_var(name=f'x{i}', var_type=mip.BINARY)
              for i in range(16*Nr)]
        A_ = [model.add_var(name=f'S{t}', var_type=mip.BINARY)
              for t in range(4*Nr)]
        sum_At = mip.xsum(A_[t] for t in range(4*Nr if ignore_last_round
                                               else 4*Nr-4))
        model.objective = sum_At
        # avoid trivial solution
        model += sum_At >= 1
        self.model = model
        self.x = x_
        self.A = A_
        # print(self.to_product_of_sum(table))
        self.sbox_constraints = None

    def to_product_of_sum(self, table=None):
        if table is None:
            table = self.table
        sz_i = self.input_size
        sz_o = self.output_size
        POS = "F = "
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

    def from_product_of_sum(self, minimized_pos):
        constraints = []
        for constr in re.findall(r'\([^\)]+\)', minimized_pos):
            pattern = {}
            for v in constr[1:-1].split('+'):
                if v.endswith("'"):
                    pattern[v[:-1]] = 1
                else:
                    pattern[v] = 0
            constraints.append(pattern.copy())
        # print(constraints)
        self.sbox_constraints = constraints

    def add_constr(self, pattern, X, Y):
        sz_i = self.input_size
        sz_o = self.output_size
        self.model += mip.xsum(
            X[i]*(1-pattern.get(f"X{sz_i-1-i}", 1))
            + (1-X[i])*pattern.get(f"X{sz_i-1-i}", 0)
            for i in range(sz_i)
        ) + mip.xsum(
            Y[j]*(1-pattern.get(f"Y{sz_o-1-j}", 1))
            + (1-Y[j])*pattern.get(f"Y{sz_o-1-j}", 0)
            for j in range(sz_o)
        ) >= 1
        # self.model.write('spn.lp')
        # exit()

    def add_sbox(self, xs, A_t, **kwargs):
        ys = kwargs.get('ys', None)
        m = self.model
        sum_xs = mip.xsum(xs[i] for i in range(len(xs)))
        m += sum_xs >= A_t
        for i in range(len(xs)):
            m += A_t >= xs[i]
        if ys is None:
            return
        if self.sbox_constraints is None:
            raise ValueError("sbox_constraints is not set, "
                             "call .from_product_of_sum first")
        for pattern in self.sbox_constraints:
            self.add_constr(pattern, xs, ys)

    def get_prob(self, sol_x_, sol_y_, sol_A_):
        DDT = self.table
        prob = Fraction(1, 1)
        for t in range(len(sol_A_) - 4):
            if prob == 0:
                return 0
            if sol_A_[t] == 1:
                in_diff = int(8*sol_x_[4*t] + 4*sol_x_[4*t+1]
                              + 2*sol_x_[4*t+2] + sol_x_[4*t+3])
                out_diff = int(8*sol_y_[4*t] + 4*sol_y_[4*t+1]
                               + 2*sol_y_[4*t+2] + sol_y_[4*t+3])
                print(in_diff, out_diff, DDT[in_diff][out_diff])
                prob *= DDT[in_diff][out_diff]
                prob /= 16
        return prob

    def print_path(self, sol_x, sol_A, **kwargs):
        Nr = self.Nr
        sol_y = kwargs.get('sol_y', None)
        for r in range(Nr):
            for i in range(16*r, 16*(r+1)):
                if i%4 == 0:
                    print(' ', end='')
                if sol_x[i] == 1:
                    print('1', end=' ')
                else:
                    print('.', end=' ')
            print()
            for i in range(4):
                t = 4*r + i
                if sol_A[t] == 1:
                    print('| S-BOX |', end='')
                else:
                    print(' '*9, end='')
            print()
            if (r == Nr - 1) or (sol_y is None):
                continue
            for i in range(16*r, 16*(r+1)):
                if i%4 == 0:
                    print(' ', end='')
                if sol_y[i] == 1:
                    print('1', end=' ')
                else:
                    print('.', end=' ')
            print()
            print(' P-BOX '.center(36, '='))

    def solve(self, **kwargs):
        max_seconds = kwargs.get('max_seconds', 30)
        max_solutions = kwargs.get('max_solutions', 1)
        sol_x_ = None
        num_sol = 0
        PBOX_INV = self.pbox_inv
        m = self.model
        NR = self.Nr
        x_ = self.x
        A_ = self.A
        self.solutions = []

        for r in range(NR):
            block_in = x_[16*r : 16*(r+1)]
            if r == NR - 1:
                for i in range(4):
                    t = 4*r + i
                    self.add_sbox(xs=block_in[4*i : 4*(i+1)], A_t=A_[t])
                break
            next_block_in = x_[16*(r+1) : 16*(r+2)]
            block_out = [next_block_in[PBOX_INV[i]] for i in range(16)]
            for i in range(4):
                t = 4*r + i
                self.add_sbox(xs=block_in[4*i : 4*(i+1)], A_t=A_[t],
                         ys=block_out[4*i : 4*(i+1)],)

        while num_sol < max_solutions:
            try:
                if sol_x_ is not None:
                    # remove previous solution
                    m += mip.xsum(
                        sol_x_[i]*(1-x_[i]) + (1-sol_x_[i])*x_[i]
                        for i in range(16*NR)
                    ) >= 1
                m.write('spn.lp')
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

            sol_x_ = [int(x_[i].x) for i in range(16*NR)]
            sol_A_ = [int(A_[t].x) for t in range(4*NR)]
            sol_y_ = sol_x_[16:]
            for r in range(NR - 1):
                next_block_in = sol_y_[16*r : 16*(r+1)]
                sol_y_[16*r : 16*(r+1)] = [next_block_in[PBOX_INV[i]]
                                           for i in range(16)]

            prob = self.get_prob(sol_x_, sol_y_, sol_A_)
            if prob == 0:
                continue
            num_sol += 1
            self.solutions.append({
                'Nr': NR,
                'sol_x': sol_x_,
                'sol_y': sol_y_,
                'sol_A': sol_A_,
                'prob': prob,
            })

        for sol in sorted(self.solutions, key=lambda x: x['prob'], reverse=True):
            self.print_path(**sol)
            print('diff. prob. = {}\n'.format(sol['prob']))


NR = 4
key = os.urandom(2*NR + 2)
cipher = basicSPN(key, nrounds=NR)
SBOX = cipher.sbox_dct['sbox']
PBOX_INV = cipher.pbox_inv
DDT = difference_distribution_table(SBOX)
# print_table(DDT)

MILP_OPTM = MilpOptim(Nr=NR, table=DDT, pbox_inv=PBOX_INV)
# print(MILP_OPTM.to_product_of_sum())
MIN_POS = """
Minimized Product of Sums:
F = (X3'+X2'+Y2'+X1'+Y0')(X3+Y3+Y2+X1'+Y1)(Y3+Y2+X1+X0'+Y0)(X3'+Y2+X1+X0+Y0)(X3+Y3'+Y2+Y1+Y0)(X2+Y2'+Y1+X0+Y0)(X2+Y2'+Y1'+X0'+Y0)(X3+Y2+Y1'+X0+Y0)(X3+X2+X1+Y1'+X0)(X3+Y2'+X1'+Y1+X0'+Y0')(X3'+Y3'+Y2'+X1+X0'+Y0')(X3'+Y3'+Y2+Y1+X0+Y0')(Y3+X2+Y2'+X1+Y1+Y0')(Y3'+X2'+X1+Y1+X0'+Y0')(Y2+X1'+Y1+X0'+Y0)(X3'+X2+Y2+X1'+X0'+Y0')(X2+Y2+X1+Y1+X0+Y0')(X3'+Y3'+Y2'+X1'+Y1+Y0)(X2'+Y2+X1+Y1'+X0'+Y0')(X3+Y3+X2'+Y2'+Y1'+Y0')(X3'+Y3+X2+Y2+Y1'+Y0')(Y3'+X2'+Y2+X1'+X0+Y0')(X3+Y3'+Y2+X1'+Y1'+Y0')(Y3'+X2+Y2'+X1'+X0+Y0')(X3'+Y3+X2+X1'+Y1+Y0)(X3+Y3+X2'+X1+Y1+X0)(X3'+Y3'+Y2+X1'+Y1'+Y0)(X3'+Y3+X2'+Y1'+X0'+Y0)(Y3'+X2'+Y2'+X1+Y1'+X0)(X3'+Y3+X2'+Y2'+X1+Y1+X0')(X3+Y2'+X1+X0+Y0')(X3'+Y2'+X1'+X0+Y0')(Y3+X2'+Y2+X1'+Y1+X0)(X3'+Y3+X2'+X1+Y1'+X0+Y0')(X3+Y3+X1+Y1+X0'+Y0)(X3'+Y3+Y2'+X1'+Y1'+Y0)(X3+Y3+X2+X1'+X0'+Y0')(X3+Y3+X2'+X1'+X0+Y0)(X3+Y3'+X2+X1+Y1'+Y0')(X3+Y3'+Y2'+Y1'+X0'+Y0)(X3'+Y3+X2'+Y2+Y1'+X0')(X3'+Y3'+Y2+X1+Y1+Y0')(X3+Y3+X2+X1+Y1+Y0')(X3'+Y3'+X1+Y1+X0+Y0)(Y3'+X2'+X1'+Y1+X0'+Y0);

"""
MILP_OPTM.from_product_of_sum(MIN_POS)
MILP_OPTM.solve()
