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
        pass

    @staticmethod
    def T(N):
        pass

    @staticmethod
    def P_b_quadratic(N):
        pass

    @staticmethod
    def T_b_quadratic(N):
        pass
