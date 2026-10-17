class dirichlet:
    @staticmethod
    def apply_left(stiffness_matrix, load_matrix, left_bc_value, left_node_index):
        stiffness_matrix[left_node_index, :] = 0.0
        stiffness_matrix[left_node_index, left_node_index] = 1.0
        load_matrix[left_node_index] = left_bc_value
        return stiffness_matrix, load_matrix

    @staticmethod
    def apply_right(stiffness_matrix, load_matrix, right_bc_value, right_node_index):
        stiffness_matrix[right_node_index, :] = 0.0
        stiffness_matrix[right_node_index, right_node_index] = 1.0
        load_matrix[right_node_index] = right_bc_value
        return stiffness_matrix, load_matrix

class neumann:
    @staticmethod
    def apply_right(load_matrix, c_val_at_b, r_b, right_node_index):
        # Adds r_b * c(b) to the rightmost load vector entry
        load_matrix[right_node_index] += r_b * c_val_at_b
        return load_matrix

class robin:
    @staticmethod
    def apply_left(stiffness_matrix, load_matrix, c_val_at_a, q_a, p_a, left_node_index):
        stiffness_matrix[left_node_index, left_node_index] -= q_a * c_val_at_a
        load_matrix[left_node_index] -= p_a * c_val_at_a
        return stiffness_matrix, load_matrix
