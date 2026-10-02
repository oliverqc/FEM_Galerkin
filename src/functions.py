import numpy as np


class c_func:
    @staticmethod
    def hw_1(x):
        return np.exp(x)


class f_func:
    @staticmethod
    def hw_1(x):
        return -np.exp(x) * (np.cos(x) - 2 * np.sin(x) - x * np.cos(x) - x * np.sin(x))


class solution_func:
    @staticmethod
    def hw_1(x):
        return x * np.cos(x)
