import argparse
import datetime
import json
import math

from src.ipcauchy.integrals.circle import evaluate_integral_circle
from src.ipcauchy.integrals.ellipse import evaluate_integral_ellipse
from src.ipcauchy.integrals.shortest_path import shortest_path_arc_integral
from src.ipcauchy.helpers.radius import optimal_radius


def get_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, required=True,
                        choices=['circle',
                                 'ellipse',
                                 'shortest_path'],
                        help='integration method')
    parser.add_argument('--file', type=str, required=True,
                        help='path to instance json')

    parser.add_argument('--N', type=int, required=False,
                        default=36,
                        help='number of angular nodes')
    parser.add_argument('--r', type=float, required=False,
                        default=0.001,
                        help='radial distance')

    return parser


def load_instance(path):
    f = open(path, 'r')
    instance = json.load(f)

    return instance['a'], instance['b'], instance['name']


if __name__ == '__main__':
    parser = get_argparse()
    cmd, other = parser.parse_known_args()

    a, b, name = load_instance(cmd.file)
    print(f'Instance Loaded: {name}')

    print(f'Integration Method Selected: {cmd.method}')
    if cmd.method == 'circle':
        start = datetime.datetime.now()
        r = optimal_radius(0, a, b)
        output = evaluate_integral_circle(a, b, r)
        end = datetime.datetime.now()

        integration_value = output[0]
        solutions = round(integration_value.real, 0)
        time = end - start
        computation_time = round(
            time.seconds
            + (time.microseconds * 1e-6), 3
        )

        print(f'Integration Value: {integration_value}')
        print(f'Solution Count: {solutions}')
        print(f'Computation Time: {computation_time} sec')

    elif cmd.method == 'ellipse':
        start = datetime.datetime.now()
        R1 = optimal_radius(0, a, b)
        R2 = optimal_radius(math.pi/2, a, b)
        output = evaluate_integral_ellipse(a, b, R1, R2)
        end = datetime.datetime.now()

        integration_value = output[0]
        solutions = round(integration_value.real, 0)
        time = end - start
        computation_time = round(
            time.seconds
            + (time.microseconds * 1e-6), 3
        )

        print(f'Integration Value: {integration_value}')
        print(f'Solution Count: {solutions}')
        print(f'Computation Time: {computation_time} sec')

    elif cmd.method == 'shortest_path':
        N, r_ = cmd.N, cmd.r
        start = datetime.datetime.now()
        r = optimal_radius(0, a, b)
        output = shortest_path_arc_integral(a, b, r, r_, N)
        end = datetime.datetime.now()

        integration_value = output[0]
        integration_time = round(output[1], 3)
        solutions = round(integration_value.real, 0)
        time = end - start
        computation_time = round(
            time.seconds
            + (time.microseconds * 1e-6), 3
        )

        print(f'Integration Value: {integration_value}')
        print(f'Solution Count: {solutions}')
        print(f'Integration Time: {integration_time} sec')
        print(f'Computation Time: {computation_time} sec')
