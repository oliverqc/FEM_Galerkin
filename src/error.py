import numpy as np


def calc_error(analytical_solution, FEM_solution):
    abs_max_difference = np.max(np.abs(analytical_solution - FEM_solution))
    return abs_max_difference
