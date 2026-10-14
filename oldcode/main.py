import numpy as np
from analytic_functions import analytic_solution
from boundary_conditions import apply_bcs_A, apply_bcs_b
from build_A_and_b import build_a, build_b
from error import error
from information_matrices import (
    build_p_b_matrix,
    build_p_matrix,
    build_t_b_matrix,
    build_t_matrix,
)


def run_fem_solver(num_elements, basis_type="linear"):
    num_total_nodes = num_elements + 1
    x = np.linspace(0, 1, num_total_nodes)

    # Build information matrices
    p_matrix = build_p_matrix(x)
    t_matrix = build_t_matrix(num_elements)
    p_b_matrix = build_p_b_matrix(basis_type, p_matrix)
    t_b_matrix = build_t_b_matrix(basis_type, t_matrix)

    num_basis_nodes = p_b_matrix.shape[1]

    # Assembly
    global_stiff_a = build_a(
        num_basis_nodes, num_elements, p_b_matrix, t_b_matrix, basis_type
    )
    global_b = build_b(
        num_basis_nodes, num_elements, p_b_matrix, t_b_matrix, basis_type
    )

    global_stiff_a = apply_bcs_A(global_stiff_a)
    global_b = apply_bcs_b(global_b)

    FEM_solution = np.linalg.solve(global_stiff_a, global_b)

    exact_solution = analytic_solution(x)
    max_error = error(exact_solution, FEM_solution)

    return max_error


def main():
    element_counts = [4, 8, 16, 32, 64, 128]

    for elements in element_counts:
        max_error = run_fem_solver(elements, basis_type="linear")
        print(f"Elements: {elements}, Max Error: {max_error:.6e}")


if __name__ == "__main__":
    main()
