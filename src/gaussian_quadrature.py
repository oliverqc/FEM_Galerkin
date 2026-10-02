import numpy as np


def generate_gauss_reference_1D(gauss_point_number):
    if gauss_point_number == 2:
        gauss_coefficient_reference_1D = [1.0, 1.0]
        gauss_point_reference_1D = [-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)]
    elif gauss_point_number == 4:
        gauss_coefficient_reference_1D = [
            0.3478548451,
            0.3478548451,
            0.6521451549,
            0.6521451549,
        ]
        gauss_point_reference_1D = [
            0.8611363116,
            -0.8611363116,
            0.3399810436,
            -0.3399810436,
        ]
    elif gauss_point_number == 8:
        gauss_coefficient_reference_1D = [
            0.1012285363,
            0.1012285363,
            0.2223810345,
            0.2223810345,
            0.3137066459,
            0.3137066459,
            0.3626837834,
            0.3626837834,
        ]
        gauss_point_reference_1D = [
            0.9602898565,
            -0.9602898565,
            0.7966664774,
            -0.7966664774,
            0.5255324099,
            -0.5255324099,
            0.1834346425,
            -0.1834346425,
        ]
    else:
        raise ValueError("Unsupported number of Gauss points. Choose 2, 4, or 8.")
    return gauss_coefficient_reference_1D, gauss_point_reference_1D


def generate_gauss_local_1D(
    gauss_coefficient_reference_1D, gauss_point_reference_1D, lower_bound, upper_bound
):
    gauss_coefficient_local_1D = []
    gauss_point_local_1D = []

    for i in range(len(gauss_coefficient_reference_1D)):
        local_coeff = (
            (upper_bound - lower_bound) * gauss_coefficient_reference_1D[i] / 2.0
        )
        local_point = (upper_bound - lower_bound) * gauss_point_reference_1D[
            i
        ] / 2.0 + (upper_bound + lower_bound) / 2.0
        gauss_coefficient_local_1D.append(local_coeff)
        gauss_point_local_1D.append(local_point)
    return gauss_coefficient_local_1D, gauss_point_local_1D


def gauss_integrate(integrand_func, left_node, right_node, gauss_point_number):
    ref_coeffs, ref_points = generate_gauss_reference_1D(gauss_point_number)
    local_coeffs, local_points = generate_gauss_local_1D(
        ref_coeffs, ref_points, left_node, right_node
    )

    integral = 0.0
    for i in range(gauss_point_number):
        integral += local_coeffs[i] * integrand_func(local_points[i])
    return integral
