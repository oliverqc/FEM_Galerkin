import numpy as np


def apply_bcs_A(global_stiff_a):
    global_stiff_a[0, :] = 0
    global_stiff_a[0, 0] = 1

    global_stiff_a[-1, :] = 0
    global_stiff_a[-1, -1] = 1

    return global_stiff_a


def apply_bcs_b(global_b):
    global_b[0] = 0.0
    global_b[-1] = np.cos(1.0)

    return global_b
