class dirichlet:
    @staticmethod
    def apply_dirichlet_linear(
        stiffness_matrix, load_matrix, left_boundary, right_boundary
    ):
        stiffness_matrix[0, :] = 0.0
        stiffness_matrix[0, 0] = 1.0
        load_matrix[0] = left_boundary

        stiffness_matrix[-1, :] = 0.0
        stiffness_matrix[-1, -1] = 1.0
        load_matrix[-1] = right_boundary

        return stiffness_matrix, load_matrix


class robin:
    @staticmethod
    def apply_robin_linear(stiffness_matrix, load_matrix):
        pass  # Implement the Robin boundary condition logic here
