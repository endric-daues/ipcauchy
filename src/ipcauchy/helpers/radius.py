import scipy.optimize as opt
from src.ipcauchy.helpers.integrand import g

def optimal_radius(angle, a, b, lower_bound=0.5, upper_bound=0.999):
    """
    Takes in the knapsack instance and returns the optimal radius using
    Brent's method to find the zero of the derivative at a specific angle
    """
    return opt.brentq(lambda r: g(r, angle, a, b), lower_bound, upper_bound)
