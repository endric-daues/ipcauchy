import cmath
import math


# FUNCTIONS TO DETERMINE THE OPTIMAL RADIUS
def summer(t, r, a):
    '''Integrand summation function'''
    temp = 0
    n = len(a)
    for i in range(0, n):
        temp += (
            a[i] * cmath.exp(1j*t) *
            ((r*cmath.exp(1j*t))**(a[i]-1))/(1-(r*cmath.exp(1j*t))**a[i])
        )

    return temp


def g(r, t, a, b):
    '''Objective function for the optimal radius'''
    return r - abs(((r ** 2) / b) * summer(t, r, a))


# FUNCTIONS TO EVALUATE THE INTEGRAND
def gamma(t, R):
    '''General gamma function using cmath and exponential'''
    return R * cmath.exp(t * 1j)


def gamma_der(t, R):
    '''derivative of above'''
    return 1j * R * cmath.exp(t * 1j)


def gamma_ellipse(t, a, b):
    '''Plain version taking two radii'''
    return a * math.cos(t) + 1j * b * math.sin(t)


def gamma_der_ellipse(t, a, b):
    '''Derivative of above'''
    return -a * math.sin(t) + 1j * b * math.cos(t)


def H(a, t, R, gamma):
    '''
    Returns generating function value. Takes in the a vector,
    the current position on the parametrized variable t and the
    gamma function.
    '''
    temp = 1
    n = len(a)
    for k in range(0, n):
        temp = temp * 1 / (1 - ((gamma(t, R))) ** a[k])

    return temp


def H_ellipse(a, t, r1, r2, gamma):
    '''
    The ellipse method requires two radii, also takes the gamma
    function and the a vector
    '''
    temp = 1
    n = len(a)
    for k in range(0, n):
        temp = temp * 1 / (1 - ((gamma_ellipse(t, r1, r2)) ** a[k]))

    return temp


def H_radial(a, R, x0, x1, y0, y1, gamma):
    '''
    Used for integrals along a ray from the inside to outside of the
    circle, perpendicular to the direction of integration
    '''
    temp = 1
    n = len(a)
    for k in range(0, n):
        temp = temp * 1 / (1 - ((gamma(R, x0, x1, y0, y1))) ** a[k])

    return temp


# FUNCTIONS TO RETRIEVE PISINGER INSTANCES
def read_pisinger_file(f):
    '''Retrieves Pisinger instances from the Pisinger File'''
    a = []
    b = []
    suba = []

    start = False
    for line in f:
        line = line.rstrip().split(' ')
        if len(line) > 0 and line[0] == 'c':
            b.append(int(line[1]))
        if line[0] == '-----':
            a.append(suba)
            suba = []
            start = False
        elif line[0].startswith('time'):
            start = True
        elif start:
            suba.append(int(line[0].split(',')[2]))

    return a, b
