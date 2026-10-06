class dirichlet:
    @staticmethod
    def apply_dirichlet(
        stiffness_matrix,
        load_matrix,
        left_bc_value,
        right_bc_value,
        left_node_index,
        right_node_index,
    ):
        stiffness_matrix[left_node_index, :] = 0.0
        stiffness_matrix[left_node_index, left_node_index] = 1.0
        load_matrix[left_node_index] = left_bc_value

        stiffness_matrix[right_node_index, :] = 0.0
        stiffness_matrix[right_node_index, right_node_index] = 1.0
        load_matrix[right_node_index] = right_bc_value

        return stiffness_matrix, load_matrix


class robin:
    @staticmethod
    def apply_robin_linear(stiffness_matrix, load_matrix):
        pass  # Implement the Robin boundary condition logic here
