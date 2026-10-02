import basis
import error
import functions
import numpy as np
from boundary_conditions import dirichlet
from build_stiffness_and_load_matrices import build_load_matrix, build_stiffness_matrix
from information_matrices import linear

# Define Parameters
N = 10  # number of internal nodes in domain
K = N + 1  # number of elements
N_b = 2  # number of basis functions per element (2 for linear)
D = 1  # dimension of the problem
gauss_point_number = 2  # number of gauss points for integration

# Domain
x = np.linspace(0, 1, N + 2)

# Generate Information Matrices
P_b = linear.P_b_linear(N)
T_b = linear.T_b_linear(N)

# Define functions to pass
c_func = functions.c_func.hw_1
f_func = functions.f_func.hw_1
basis_type = basis.linear

K_values = [4, 8, 16, 32, 64, 128]
print(f"{'h':<8} | {'Maximum absolute error at all nodes'}")
print("-" * 45)

for K in K_values:
    N = K - 1  # Calculate number of internal nodes
    # Domain
    x = np.linspace(0, 1, N + 2)
    P_b = linear.P_b_linear(N)
    T_b = linear.T_b_linear(N)

    # Build Stiffness and Load Matrices
    stiffness_matrix = build_stiffness_matrix(
        N, N_b, K, P_b, T_b, c_func, basis_type, gauss_point_number
    )
    load_matrix = build_load_matrix(
        N, N_b, K, P_b, T_b, f_func, basis_type, gauss_point_number
    )

    # Apply Boundary Conditions
    left_boundary = functions.solution_func.hw_1(0.0)
    right_boundary = functions.solution_func.hw_1(1.0)
    stiffness_matrix, load_matrix = dirichlet.apply_dirichlet_linear(
        stiffness_matrix, load_matrix, left_boundary, right_boundary
    )

    # Solve
    analytical_solution = functions.solution_func.hw_1(x)
    FEM_solution = np.linalg.solve(stiffness_matrix, load_matrix)

    # Calculate Error
    final_error = error.calc_error(
        analytical_solution.flatten(), FEM_solution.flatten()
    )

    h_string = f"1/{K}"
    print(f"{h_string:<8} | {final_error:.4e}")
