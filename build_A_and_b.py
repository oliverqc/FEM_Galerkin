import numpy as np


def local_a(k, p_b_matrix, t_b_matrix, basis_type):
    if basis_type == "linear":
        node_left = t_b_matrix[0, k]
        node_right = t_b_matrix[1, k]

        x_left = p_b_matrix[0, node_left]
        x_right = p_b_matrix[0, node_right]
        h = x_right - x_left

        local_stiff_a = np.zeros((2, 2))
        dpsi_dx = np.array([-1 / h, 1.0 / h])
        num_gauss_points = 3
        xi_points, weights = np.polynomial.legendre.leggauss(num_gauss_points)

        for i in range(num_gauss_points):
            x_phys = (h / 2.0) * xi_points[i] + (x_right + x_left) / 2.0
            w_phys = (h / 2.0) * weights[i]

            c_val = np.exp(x_phys)

            for beta in range(2):
                for alpha in range(2):
                    local_stiff_a[beta, alpha] += (
                        w_phys * c_val * dpsi_dx[beta] * dpsi_dx[alpha]
                    )

        return local_stiff_a

    elif basis_type == "quadratic":
        print("Quadratic basis type not implemented yet.")


def build_a(num_basis_nodes, num_elements, p_b_matrix, t_b_matrix, basis_type):
    global_stiff_a = np.zeros((num_basis_nodes, num_basis_nodes))

    for k in range(num_elements):
        local_stiff_a = local_a(k, p_b_matrix, t_b_matrix, basis_type)
        num_local_basis = t_b_matrix.shape[0]  # num local basis functions per element

        for beta in range(num_local_basis):
            for alpha in range(num_local_basis):

                m = t_b_matrix[beta, k]
                n = t_b_matrix[alpha, k]

                global_stiff_a[m, n] += local_stiff_a[beta, alpha]

    return global_stiff_a


def f_func(x):
    return -np.exp(x) * (np.cos(x) - 2 * np.sin(x) - x * np.cos(x) - x * np.sin(x))


def local_b(k, p_b_matrix, t_b_matrix, basis_type):
    if basis_type == "linear":
        node_left = t_b_matrix[0, k]
        node_right = t_b_matrix[1, k]

        x_left = p_b_matrix[0, node_left]
        x_right = p_b_matrix[0, node_right]
        h = x_right - x_left

        local_load_b = np.zeros(2)

        def psi_1(x):
            return (x_right - x) / h

        def psi_2(x):
            return (x - x_left) / h

        num_gauss_points = 5
        xi_points, weights = np.polynomial.legendre.leggauss(num_gauss_points)

        for i in range(num_gauss_points):
            x_phys = (h / 2.0) * xi_points[i] + (x_right + x_left) / 2.0
            w_phys = (h / 2.0) * weights[i]

            f_val = f_func(x_phys)

            psi_vals = np.array([psi_1(x_phys), psi_2(x_phys)])

            for beta in range(2):
                local_load_b[beta] += w_phys * f_val * psi_vals[beta]

        return local_load_b

    elif basis_type == "quadratic":
        print("Quadratic basis type not implemented yet.")


def build_b(num_basis_nodes, num_elements, p_b_matrix, t_b_matrix, basis_type):
    global_b = np.zeros(num_basis_nodes)

    for k in range(num_elements):
        loc_b = local_b(k, p_b_matrix, t_b_matrix, basis_type)
        num_local_basis = t_b_matrix.shape[0]

        # Route local vector into global load vector
        for beta in range(num_local_basis):
            m = t_b_matrix[beta, k]  # Global row index
            global_b[m] += loc_b[beta]

    return global_b
