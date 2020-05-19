# -*- coding: utf-8 -*-
"""
base class and some other utils for MILP Optimizer
"""
import sys
import subprocess
from itertools import product

from mip import Model, BINARY, xsum, OptimizationStatus

__all__ = ['to_pla', 'from_pla', 'simplify', 'parse_constr', 'MilpOptim',
           'linear_approximation_table']


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


def to_pla(table, output_file=None):
    """
    write table modeling to PLA file
    """
    nrows = len(table)
    ncols = len(table[0])
    size = (nrows-1).bit_length() + (ncols-1).bit_length()
    table = [table[i][j] for i in range(nrows) for j in range(ncols)]
    data = ".i {}\n".format(size)
    data += ".o 1\n"
    data += ".p {}\n".format(nrows * ncols)
    data += ".ilb {} \n".format(
        " ".join("x{}".format(i) for i in range(size)))
    data += ".ob y0\n"
    data += ".type fr\n"
    for idx in range(1 << size):
        data += "{0:0{1}b} {2}\n".format(idx, size, 1 if table[idx] == 0 else 0)
    data += ".e\n"
    if output_file is None:
        sys.stdout.write(data)
        return
    open(output_file, 'w', buffering=1).write(data)

def from_pla(input_file):
    """
    read PLA file
    """
    pla_file = open(input_file, 'r')
    data = [line[:-3] for line in pla_file if not line.startswith('.')]
    return data

def simplify(file_i, file_o):
    """
    run command for Espresso logic minimization
    """
    command = "espresso {} > {}".format(file_i, file_o)
    subprocess.run(command, shell=True, check=True)

def parse_constr(minimized_cons, sbox_size=None):
    """
    build constraints for a given (minimized) POS

    INPUT:

    - ``minimized_cons`` - list of constraints string
    - ``sbox_size`` - tuple (sbox input size, sbox output size)
    """
    if sbox_size is None:
        sz_i = sz_o = len(minimized_cons[0])//2
    else:
        sz_i, sz_o = sbox_size
    sbox_size = len(minimized_cons[0])//2
    constraints = []
    for con in minimized_cons:
        pattern = {}
        for i, ch in enumerate(con[:sz_i]):
            if ch != '-':
                pattern['X{}'.format(i)] = int(ch)
        for j, ch in enumerate(con[-sz_o:]):
            if ch != '-':
                pattern['Y{}'.format(j)] = int(ch)
        constraints.append(pattern.copy())
    return constraints


class MilpOptim(Model):
    """
    MILP Optimizer for S-Box based cipher
    """

    def __repr__(self):
        return "MIP ({}) model with {} constraints".format(
            self.solver_name, self.num_rows)

    def new_var(self, name, nbits):
        """
        generate bit-vecter, represented by array of BINARY
        """
        var = [self.add_var(name=f'{name}_{i}', var_type=BINARY)
               for i in range(nbits)]
        return var

    def add_xor(self, A, B, C):
        r"""
        add xor constraint for expr ``A \xor B = C``
        """
        if not isinstance(A, list):
            A, B, C = [A], [B], [C]
        for a, b, c in zip(A, B, C):
            self += a + b + (1-c) >= 1
            self += a + (1-b) + c >= 1
            self += (1-a) + b + c >= 1
            self += (1-a) + (1-b) + (1-c) >= 1

    def add_impos(self, pattern, X, Y):
        """
        remove impossible case for a given pattern
        """
        sz_i = len(X)
        sz_o = len(Y)
        self += xsum(
            X[i] * (1-pattern.get(f"X{i}", 1))
            + (1-X[i]) * pattern.get(f"X{i}", 0)
            for i in range(sz_i)
        ) + xsum(
            Y[j] * (1-pattern.get(f"Y{j}", 1))
            + (1-Y[j]) * pattern.get(f"Y{j}", 0)
            for j in range(sz_o)
        ) >= 1

    def check(self, **kwds):
        """
        optimize the model and output the status,
        raise ValueError when no solution found
        """
        status = self.optimize(**kwds)
        if status == OptimizationStatus.OPTIMAL:
            sys.stdout.write('optimal solution cost {} found\n'.format(
                self.objective_value))
        elif status == OptimizationStatus.FEASIBLE:
            sys.stdout.write('sol.cost {} found, best possible: {:.4f}\n'.format(
                self.objective_value, self.objective_bound))
        elif status == OptimizationStatus.NO_SOLUTION_FOUND:
            raise ValueError('no feasible solution found, '\
                             'lower bound is: {}'.format(self.objective_bound))
        else:
            raise ValueError("not considered status code")

    def solve(self, num_solutions=1, max_seconds=300, target_vars=None, write_file=None):
        """
        enumerate multiple optimal solutions
        """
        if target_vars is None:
            target_vars = self.vars
        if write_file is not None:
            self.write(write_file)
        solutions = []
        count = 0
        while count < num_solutions:
            if target_vars[0].x is not None:
                # remove previous solution
                self += xsum(
                    int(v.x)*(1-v) + int(1-v.x)*v
                    for v in target_vars
                ) >= 1
            try:
                self.check(max_seconds=max_seconds)
            except ValueError as e:
                sys.stderr.write(str(e)+'\n')
                break
            count += 1
            solution = dict((v.name, v.x) for v in target_vars)
            solutions.append(solution)
        return solutions


if __name__ == '__main__':
    import mip

    table = [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
    ]

    to_pla(table, '/tmp/1.pla')
    simplify('/tmp/1.pla', '/tmp/2.pla')
    DATA_SP = from_pla('/tmp/2.pla')
    print(DATA_SP)
    CONSTR = parse_constr(DATA_SP)

    SLVR = MilpOptim(sense=mip.MINIMIZE, solver_name=mip.CBC)
    SLVR.verbose = 0
    SLVR.preprocess = 1
    SLVR.threads = 4

    x = SLVR.new_var(name='x', nbits=2)
    y = SLVR.new_var(name='y', nbits=2)
    sum_x = mip.xsum(x)
    sum_y = mip.xsum(y)

    SLVR.objective = sum_x

    SLVR += sum_y == 1 # select from column 1(0b01) and 2(0b10)
    for pattern in CONSTR:
        SLVR.add_impos(pattern, x, y)

    sols = SLVR.solve(num_solutions=4, write_file='tst.lp')

    for sol in sols:
        print(sol)

    # solution:
    # {'x_0': 0.0, 'x_1': 0.0, 'y_0': 0.0, 'y_1': 1.0}
    #   ==> table[0b00][0b01] != 0
    # {'x_0': 0.0, 'x_1': 1.0, 'y_0': 1.0, 'y_1': 0.0}
    #   ==> table[0b01][0b10] != 0
