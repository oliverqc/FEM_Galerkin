import numpy as np


def build_p_matrix(x):
    return np.array([x])


def build_t_matrix(num_elements):
    t_matrix = np.zeros((2, num_elements), dtype=int)
    for k in range(num_elements):
        t_matrix[0, k] = k
        t_matrix[1, k] = k + 1
    return t_matrix


def build_p_b_matrix(basis_type, p_matrix):
    if basis_type == "linear":
        return p_matrix.copy()
    elif basis_type == "quadratic":
        pass


def build_t_b_matrix(basis_type, t_matrix):
    if basis_type == "linear":
        return t_matrix.copy()
    elif basis_type == "quadratic":
        pass
