class linear:
    @staticmethod
    def derivative(beta, left_node, right_node, x):
        h = right_node - left_node
        if beta == 1:
            return -1.0 / h
        elif beta == 2:
            return 1.0 / h
        else:
            raise ValueError("Invalid beta index for linear basis function")

    @staticmethod
    def __call__(beta, left_node, right_node, x):
        h = right_node - left_node
        if beta == 1:
            return (right_node - x) / h
        elif beta == 2:
            return (x - left_node) / h
        else:
            raise ValueError("Invalid beta index for linear basis function")


class quadratic:
    @staticmethod
    def derivative(beta, left_node, right_node, x):
        h = right_node - left_node
        if beta == 1:
            return (1.0 / h) * (4.0 * ((x - left_node) / h) - 3.0)
        elif beta == 2:
            return (1.0 / h) * (4.0 * ((x - left_node) / h) - 1.0)
        elif beta == 3:
            return (1.0 / h) * (-8.0 * ((x - left_node) / h) + 4.0)
        else:
            raise ValueError("Invalid beta index for quadratic basis function")

    @staticmethod
    def __call__(beta, left_node, right_node, x):
        h = right_node - left_node
        if beta == 1:
            return 2.0 * ((x - left_node) / h) ** 2 - 3.0 * ((x - left_node) / h) + 1.0
        elif beta == 2:
            return 2.0 * ((x - left_node) / h) ** 2 - ((x - left_node) / h)
        elif beta == 3:
            return -4.0 * ((x - left_node) / h) ** 2 + 4.0 * ((x - left_node) / h)
        else:
            raise ValueError("Invalid beta index for quadratic basis function")
