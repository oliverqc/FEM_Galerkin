import numpy as np


class linear:
    @staticmethod
    def P(N):
        return np.array([np.linspace(0, 1, N + 2)])

    @staticmethod
    def T(N):
        K = N + 1
        return np.array([[k for k in range(K)], [k + 1 for k in range(K)]])

    @staticmethod
    def P_b_linear(N):
        return linear.P(N)

    @staticmethod
    def T_b_linear(N):
        return linear.T(N)


class quadratic:
    @staticmethod
    def P(N):
        return np.array([np.linspace(0, 1, N + 2)])

    @staticmethod
    def T(N):
        K = N + 1
        return np.array([[k for k in range(K)], [k + 1 for k in range(K)]])

    @staticmethod
    def P_b_quadratic(N):
        points = np.linspace(0, 1, N + 2)
        midpoints = (points[:-1] + points[1:]) / 2.0
        return np.array([np.concatenate((points, midpoints))])

    @staticmethod
    def T_b_quadratic(N):
        K = N + 1
        row1 = [k for k in range(K)]  # left point
        row2 = [k + 1 for k in range(K)]  # right point
        row3 = [(N + 2) + k for k in range(K)]  # midpoint
        return np.array([row1, row2, row3])
