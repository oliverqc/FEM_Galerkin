import numpy as np
from gaussian_quadrature import gauss_integrate


def build_local_stiffness_matrix(
    c_func, basis_type, left_node, right_node, gauss_point_number, N_b
):
    local_stiffness_matrix = np.zeros((N_b, N_b))
    for beta in range(1, N_b + 1):
        for alpha in range(1, N_b + 1):

            def integrand(x):
                d_psi_beta = basis_type.derivative(beta, left_node, right_node, x)
                d_psi_alpha = basis_type.derivative(alpha, left_node, right_node, x)
                return c_func(x) * d_psi_beta * d_psi_alpha

            local_stiffness_matrix[beta - 1, alpha - 1] = gauss_integrate(
                integrand, left_node, right_node, gauss_point_number
            )
    return local_stiffness_matrix


def build_stiffness_matrix(N, N_b, K, P_b, T_b, c_func, basis_type, gauss_point_number):
    # Using np.zeros to create a dense matrix (numpy.linalg.solve requires dense matrices)
    stiffness_matrix = np.zeros((N + 2, N + 2))
    for k in range(K):
        left_node = P_b[0, T_b[0, k]]
        right_node = P_b[0, T_b[1, k]]

        local_stiffness_matrix = build_local_stiffness_matrix(
            c_func, basis_type, left_node, right_node, gauss_point_number, N_b
        )

        for beta in range(1, N_b + 1):
            for alpha in range(1, N_b + 1):
                stiffness_matrix[
                    T_b[beta - 1, k], T_b[alpha - 1, k]
                ] += local_stiffness_matrix[beta - 1, alpha - 1]

    return stiffness_matrix


def build_load_matrix(N, N_b, K, P_b, T_b, f_func, basis_type, gauss_point_number):
    load_matrix = np.zeros((N + 2, 1))
    for k in range(K):
        left_node = P_b[0, T_b[0, k]]
        right_node = P_b[0, T_b[1, k]]

        for beta in range(1, N_b + 1):

            def integrand(x):
                psi_beta = basis_type.__call__(beta, left_node, right_node, x)
                return f_func(x) * psi_beta

            load_matrix[T_b[beta - 1, k]] += gauss_integrate(
                integrand, left_node, right_node, gauss_point_number
            )

    return load_matrix
