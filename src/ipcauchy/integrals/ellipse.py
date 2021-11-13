import numpy as np
import math
import quadpy

from src.ipcauchy.helpers.integrand import (
    H_ellipse,
    gamma_ellipse,
    gamma_der_ellipse
)

def integrand_ellipse_adaptive(t, a, b, r1, r2, gamma, gamma_der):
    """
    Vectorized integrand function returns integrand value along
    elliptical path
    """
    input_shape = t.shape
    sol = []
    length = input_shape[0]
    for i in range(0, length):
        sub_sol = []
        input_sub_shape = t[i].shape
        for k in t[i]:
            sub_sol.append(
                H_ellipse(a, k, r1, r2, gamma)
                * (1 / (gamma(k, r1, r2) ** (b + 1)))
                * gamma_der(k, r1, r2)
            )
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)

    return sol_array


def evaluate_integral_ellipse(a, b, R1, R2,
                              alpha=0, beta=2*math.pi, error=1):
    """
    Evaluates the integral (3) on an elliptic path parameterized by the
    radii R1 (real axis) and R2. Uses the quadpy library to compute the
    integral with an adaptive quadrature approach.
    """
    val = quadpy.line_segment.integrate_adaptive(
        lambda t: integrand_ellipse_adaptive(t, a, b,
                                             R1, R2, gamma_ellipse,
                                             gamma_der_ellipse),
        [alpha, beta],
        error)

    return ((1/(2*math.pi*1j)) * val[0], val[1])
