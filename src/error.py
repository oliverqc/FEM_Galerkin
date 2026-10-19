import numpy as np
from gaussian_quadrature import gauss_integrate

def compute_errors(K, T_b, P_b, basis_type, gauss_point_number, N_b, FEM_solution, u_func, u_deriv_func):
    L2_error_sq = 0.0
    H1_semi_error_sq = 0.0
    L_inf_error = 0.0

    for k in range(K):
        left_node = P_b[0, T_b[0, k]]
        right_node = P_b[0, T_b[1, k]]

        u_h_local = [FEM_solution[T_b[beta - 1, k]] for beta in range(1, N_b + 1)]

        # L2 Norm integrand: (u - u_h)^2
        def l2_integrand(x):
            u_h_val = sum(u_h_local[beta - 1] * basis_type.__call__(beta, left_node, right_node, x) for beta in range(1, N_b + 1))
            return (u_func(x) - u_h_val)**2

        # H1 Semi-norm integrand: (u' - u'_h)^2
        def h1_semi_integrand(x):
            u_h_deriv_val = sum(u_h_local[beta - 1] * basis_type.derivative(beta, left_node, right_node, x) for beta in range(1, N_b + 1))
            return (u_deriv_func(x) - u_h_deriv_val)**2

        L2_error_sq += gauss_integrate(l2_integrand, left_node, right_node, gauss_point_number)
        H1_semi_error_sq += gauss_integrate(h1_semi_integrand, left_node, right_node, gauss_point_number)

        # L_inf norm: max absolute difference over the element. 
        sample_points = np.linspace(left_node, right_node, 20)
        for x in sample_points:
            u_h_val = sum(u_h_local[beta - 1] * basis_type.__call__(beta, left_node, right_node, x) for beta in range(1, N_b + 1))
            diff = abs(u_func(x) - u_h_val)
            if diff > L_inf_error:
                L_inf_error = diff

    L2_error = np.sqrt(L2_error_sq)
    H1_semi_error = np.sqrt(H1_semi_error_sq)
    H1_error = np.sqrt(L2_error_sq + H1_semi_error_sq)

    return L_inf_error, L2_error, H1_semi_error, H1_error
