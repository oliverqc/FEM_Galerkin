import basis
import error
import functions
import numpy as np
from boundary_conditions import dirichlet
from build_stiffness_and_load_matrices import build_load_matrix, build_stiffness_matrix
from information_matrices import linear, quadratic

# Define Parameters
N = 10
K = N + 1  # number of elements
N_b = 3  # number of basis functions per element (2 for linear)
D = 1  # dimension of the problem
gauss_point_number = 4  # number of gauss points for integration

# Define functions to pass
c_func = functions.c_func.hw_1
f_func = functions.f_func.hw_1
basis_type = basis.quadratic  # linear or quadratic

K_values = [4, 8, 16, 32, 64, 128]
print(f"{'h':<8} | {'Maximum absolute error at all nodes'}")
print("-" * 45)

errors = []

for K in K_values:
    N = K - 1  # Calculate number of internal nodes
    # Domain
    P_b = quadratic.P_b_quadratic(N)  # linear or quadratic
    T_b = quadratic.T_b_quadratic(N)

    # Build Stiffness and Load Matrices
    stiffness_matrix = build_stiffness_matrix(
        N, N_b, K, P_b, T_b, c_func, basis_type, gauss_point_number
    )
    load_matrix = build_load_matrix(
        N, N_b, K, P_b, T_b, f_func, basis_type, gauss_point_number
    )

    # Apply Boundary Conditions
    left_bc_value = functions.solution_func.hw_1(0.0)
    right_bc_value = functions.solution_func.hw_1(1.0)
    left_node_index = 0
    right_node_index = N + 1

    stiffness_matrix, load_matrix = dirichlet.apply_dirichlet(
        stiffness_matrix,
        load_matrix,
        left_bc_value,
        right_bc_value,
        left_node_index,
        right_node_index,
    )

    # Solve
    analytical_solution = functions.solution_func.hw_1(P_b[0])
    FEM_solution = np.linalg.solve(stiffness_matrix, load_matrix)

    # Calculate Error
    final_error = error.calc_error(
        analytical_solution.flatten(), FEM_solution.flatten()
    )

    errors.append(final_error)

    h_string = f"1/{K}"
    print(f"{h_string:<8} | {final_error:.4e}")

order_of_convergence = np.log(errors[-2] / errors[-1]) / np.log(
    K_values[-1] / K_values[-2]
)
print(order_of_convergence)
