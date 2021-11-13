import datetime
import math
import numpy as np
import quadpy
from dijkstar import Graph, find_path

from src.ipcauchy.helpers.shortest_path import (
    generate_graph,
    fill_edges,
    convert,
    gamma_radial,
    gamma_der_radial,
    gamma_circle_rad,
    gamma_der_circle_rad,
)

from src.ipcauchy.helpers.integrand import (
    H,
    H_radial
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


def integrand_adaptive_line(t, a, b, x1, y1, x2, y2, gamma, gamma_der):
    """Along a line between two points With MP """
    input_shape = t.shape
    sol = []
    length = input_shape[0]
    for i in range(0, length):
        sub_sol = []
        input_sub_shape = t[i].shape
        for k in t[i]:
            sub_sol.append(
                H_radial(a, k, x1, x2, y1, y2, gamma)
                * (1 / (gamma(k, x1, x2, y1, y2) ** (b + 1)))
                * gamma_der(k, x1, x2, y1, y2)
            )
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)

    return sol_array


def angle_adjustment_integral(a, b, alpha, beta, R, gamma, gamma_der, error):
    """Adaptive integral from theta1 to theta2 along a constant radiu  """
    val = quadpy.line_segment.integrate_adaptive(
        lambda t: integrand_adaptive(t, a, b, R, gamma, gamma_der),
        [alpha, beta],
        error)

    return ((1/(2*math.pi * 1j)) * val[0], val[1])


def line_integral(a, b, x1, y1, x2, y2, gamma, gamma_der, error):
    """
    Adaptive integral from r1 to r2 along the same radial line
    """
    val = quadpy.line_segment.integrate_adaptive(
        lambda t: integrand_adaptive_line(t, a, b,
                                          x1, y1, x2,
                                          y2, gamma, gamma_der),
        [0, 1],
        error)
    return ((1/(2*math.pi * 1j)) * val[0], val[1])


def shortest_path_arc_integral(a, b, start_radius, radius_delta, N, error=1):
    """
    Solves the shortest path given a,b and outputs the integral, moves along
    arcs not straight lines. The max radius should be larger than the min
    radius such that there are not zero paths. THe precision of the
    max_radius should be larger than the delta since this sometimes
    throws an error.
    """
    start_radius = round(start_radius, 15)
    epsilon = round((2*math.pi)/720, 15)
    path_start = (start_radius, epsilon)
    angles = [epsilon]
    current = epsilon
    delta = round((2*math.pi - 2*epsilon)/(N), 15)

    for i in range(0, N):
        current += delta
        angles.append(current)

    radii = []
    r = start_radius
    while r < 0.9999:
        radii.append(r)
        r += radius_delta
    path_end = (start_radius, angles[-1])

    nodes = generate_graph(radii, angles)
    graph_ = Graph()
    graph = fill_edges(graph_, radii, angles, nodes, a, b)

    cost_func = lambda u, v, e, prev_e: e['cost']
    path_nodes = list(find_path(graph, nodes[path_start], nodes[path_end],
                                cost_func=cost_func).nodes)
    edge_labels = {y: x for x, y in nodes.items()}

    currentDT = datetime.datetime.now()
    integral = 0
    for n in range(1, len(path_nodes)):
        r1 = edge_labels[path_nodes[n-1]][0]
        r2 = edge_labels[path_nodes[n]][0]

        theta1 = edge_labels[path_nodes[n-1]][1]
        theta2 = edge_labels[path_nodes[n]][1]

        # if there is no radius adjustment, move along an arc
        if r1 == r2:
            integral += (
                angle_adjustment_integral(a, b, theta1,
                                          theta2, r1, gamma_circle_rad,
                                          gamma_der_circle_rad, error)[0]
            )
        else:
            z1 = convert(r1, theta1)
            z2 = convert(r2, theta2)
            x1 = z1.real
            y1 = z1.imag
            x2 = z2.real
            y2 = z2.imag
            integral += line_integral(a, b, x1,
                                      y1, x2, y2,
                                      gamma_radial,
                                      gamma_der_radial,
                                      error)[0]

    # add final integral between start and end node
    integral += angle_adjustment_integral(a, b, 0,
                                          angles[0], radii[0],
                                          gamma_circle_rad,
                                          gamma_der_circle_rad, error)[0]
    integral += angle_adjustment_integral(a, b, angles[-1],
                                          2*math.pi, radii[0],
                                          gamma_circle_rad,
                                          gamma_der_circle_rad, error)[0]

    currentDT2 = datetime.datetime.now()
    timed = currentDT2 - currentDT
    t = timed.seconds + (timed.microseconds*1e-6)

    return integral, t
