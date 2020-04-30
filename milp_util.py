# -*- coding: utf-8 -*-
"""
base class and some other utils for MILP Optimizer
"""
import subprocess

import mip

__all__ = ['mip', 'to_pla', 'from_pla', 'simplify', 'parse_constr', 'MilpOptim']


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
        print(data)
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

def parse_constr(minimized_cons):
    """
    build constraints for a given (minimized) POS
    """
    sbox_size = len(minimized_cons[0])//2
    constraints = []
    for con in minimized_cons:
        pattern = {}
        for i, ch in enumerate(con[:sbox_size]):
            if ch != '-':
                pattern['X{}'.format(sbox_size-i)] = int(ch)
        for j, ch in enumerate(con[-sbox_size:]):
            if ch != '-':
                pattern['Y{}'.format(sbox_size-j)] = int(ch)
        constraints.append(pattern.copy())
    return constraints


class MilpOptim:
    """
    base class of MILP Optimizer
    """

    def __init__(self, **kwargs):
        model = mip.Model(sense=mip.MINIMIZE, solver_name=mip.CBC)
        model.verbose = kwargs.get('verbose', 0)
        model.preprocess = kwargs.get('preprocess', 1)
        model.threads = kwargs.get('threads', 4)
        self.model = model

    def __repr__(self):
        return "MIP ({}) model with {} constraints".format(
            self.model.solver_name, self.model.num_rows)

    def optimize(self, max_seconds):
        model = self.model
        status = model.optimize(max_seconds=max_seconds)
        if status == mip.OptimizationStatus.OPTIMAL:
            print('optimal solution cost {} found'.format(
                model.objective_value))
        elif status == mip.OptimizationStatus.FEASIBLE:
            print('sol.cost {} found, best possible: {}'.format(
                model.objective_value, model.objective_bound))
        elif status == mip.OptimizationStatus.NO_SOLUTION_FOUND:
            raise ValueError('no feasible solution found, '\
                             'lower bound is: {}'.format(model.objective_bound))
        else:
            raise ValueError("not considered status code")


if __name__ == '__main__':
    LAT = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 2, 2, 0, 0, 2, 6, 2, 2, 0, 0, 2, 2, 0, 0],
        [0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 6, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 2, 6, 2, 2, 2, 2, 2, 2],
        [0, 2, 0, 2, 2, 4, 2, 0, 0, 2, 0, 2, 2, 4, 2, 0],
        [0, 2, 2, 0, 2, 0, 4, 2, 2, 0, 4, 2, 0, 2, 2, 0],
        [0, 2, 2, 4, 2, 0, 0, 2, 0, 2, 2, 4, 2, 0, 0, 2],
        [0, 2, 0, 2, 2, 4, 2, 0, 2, 0, 2, 0, 4, 2, 0, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 6],
        [0, 0, 2, 2, 0, 0, 2, 2, 4, 0, 2, 2, 0, 4, 2, 2],
        [0, 4, 2, 2, 4, 0, 2, 2, 2, 2, 0, 0, 2, 2, 0, 0],
        [0, 4, 0, 4, 4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 4, 2, 2, 0, 2, 0, 2, 0, 2, 4, 0, 2, 0, 2],
        [0, 2, 2, 0, 2, 4, 0, 2, 4, 2, 2, 0, 2, 0, 0, 2],
        [0, 2, 2, 0, 2, 4, 0, 2, 2, 0, 0, 2, 4, 2, 2, 0],
        [0, 2, 4, 2, 2, 0, 2, 0, 0, 2, 4, 2, 2, 0, 2, 0],
    ]
    to_pla(LAT, 'lat.pla')
    simplify('lat.pla', 'lat_sp.pla')
    DATA_SP = from_pla('lat_sp.pla')
    print(DATA_SP)
