import basis
import error
import functions
import numpy as np
from boundary_conditions import dirichlet, neumann, robin
from build_stiffness_and_load_matrices import build_load_matrix, build_stiffness_matrix
from information_matrices import linear, quadratic

N_b = 3  # number of basis functions per element (2 for linear)
D = 1  # dimension of the problem
gauss_point_number = 4  # number of gauss points for integration

# Define functions to pass
c_func = functions.c_func.hw_1
f_func = functions.f_func.hw_1
u_func = functions.solution_func.hw_1
u_deriv_func = functions.solution_func.hw_1_deriv
basis_type = basis.quadratic  # linear or quadratic

K_values = [4, 8, 16, 32, 64, 128]

print(f"{'h':<8} | {'L_inf Error':<15} | {'L2 Error':<15} | {'H1 Semi Error':<15} | {'H1 Error':<15}")
print("-" * 80)

for K in K_values:
    N = K - 1  # Calculate number of internal nodes
    
    P_b = quadratic.P_b_quadratic(N)
    T_b = quadratic.T_b_quadratic(N)

    stiffness_matrix = build_stiffness_matrix(
        N, N_b, K, P_b, T_b, c_func, basis_type, gauss_point_number
    )
    load_matrix = build_load_matrix(
        N, N_b, K, P_b, T_b, f_func, basis_type, gauss_point_number
    )

    # Apply Boundary Conditions
    left_bc_value = u_func(0.0)
    right_bc_value = u_func(1.0)
    left_node_index = 0
    right_node_index = N + 1

    # Example 3
    # Left Dirichlet: u(0) = 0
    
    #stiffness_matrix, load_matrix = dirichlet.apply_left(
    #    stiffness_matrix, load_matrix, 0.0, left_node_index
    #)

    # Right Neumann: u'(1) = cos(1) - sin(1)

    #c_b = c_func(1.0)
    #r_b = np.cos(1.0) - np.sin(1.0)
    #load_matrix = neumann.apply_right(
    #    load_matrix, c_b, r_b, right_node_index
    #)

    # Example 4
    # Left Robin

    c_a = c_func(0.0)
    stiffness_matrix, load_matrix = robin.apply_left(
        stiffness_matrix, load_matrix, c_a, 1.0, 1.0, left_node_index
    )

    # Right Dirichlet

    stiffness_matrix, load_matrix = dirichlet.apply_right(
        stiffness_matrix, load_matrix, u_func(1.0), right_node_index
    )

    # Solve
    FEM_solution = np.linalg.solve(stiffness_matrix, load_matrix)

    # Calculate Continuous Errors over the entire domain
    L_inf, L2, H1_semi, H1 = error.compute_errors(
        K, T_b, P_b, basis_type, gauss_point_number, N_b, 
        FEM_solution.flatten(), u_func, u_deriv_func
    )

    h_string = f"1/{K}"
    print(f"{h_string:<8} | {L_inf:.4e}{' '*4} | {L2:.4e}{' '*4} | {H1_semi:.4e}{' '*4} | {H1:.4e}")
