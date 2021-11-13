import math
import numpy as np
import quadpy

from src.ipcauchy.helpers.integrand import (
    H,
    gamma,
    gamma_der
)


def integrand_adaptive(t, a, b, R, gamma, gamma_der):
    """
    Vectorized integrand function returns integrand value along
    circular path
    """
    input_shape = t.shape
    sol = []
    length = input_shape[0]
    for i in range(0, length):
        sub_sol = []
        input_sub_shape = t[i].shape
        for k in t[i]:
            sub_sol.append(
                H(a, k, R, gamma)
                * (1 / (gamma(k, R) ** (b+1)))
                * gamma_der(k, R)
            )
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)

    return sol_array


def evaluate_integral_circle(a, b, R, alpha=0, beta=2*math.pi, error=1):
    """
    Evaluates the integral (3) on a cricular path parameterized by the
    radius R. Uses the quadpy library to compute the integral with an adaptive
    quadrature approach.
    """
    val = quadpy.line_segment.integrate_adaptive(
        lambda t: integrand_adaptive(t, a, b, R, gamma, gamma_der),
        [alpha, beta],
        error)

    return (1/(2*math.pi*1j) * val[0], val[1])
