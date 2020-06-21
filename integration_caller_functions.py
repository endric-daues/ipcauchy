import math
import quadpy
import scipy.optimize as opt
from integrand_functions import integrand_adaptive
from integrand_functions import integrand_ellipse_adaptive
from helpers import gamma,gamma_der, gamma_ellipse,gamma_der_ellipse
from shortest_path_helpers import fill_edges,generate_graph,convert
from shortest_path_helpers import angle_adjustment_integral, line_integral
from shortest_path_helpers import gamma_circle_rad,gamma_der_circle_rad
from shortest_path_helpers import gamma_radial,gamma_der_radial
from helpers import g

from dijkstar import Graph, find_path
import datetime

def optimal_radius(angle,a,b,lower_bound=0.5,upper_bound=0.999):
    '''Takes in the knapsack instance and returns the optimal radius using
    Brent´s method to find the zero of the derivative at a specific angle'''
    return opt.brentq(lambda r: g(r,angle,a,b), lower_bound, upper_bound)

def count_solutions_circle(a,b,R,alpha=0,beta=2*math.pi,error=1):
    '''Evaluates the integral (3) on a cricular path parameterized by the 
    radius R. Uses the quadpy library to compute the integral with an adaptive
    quadrature approach. Returns the integral value as well as an error estimate.'''

    val = quadpy.line_segment.integrate_adaptive(lambda t: integrand_adaptive(t,a,b,R,gamma,gamma_der),[alpha,beta],error)
    return ((1/(2*math.pi*1j)) * val[0],val[1])


def count_solutions_ellipse(a,b,R1,R2,alpha=0,beta=2*math.pi,error=1):
    '''Evaluates the integral (3) on an elliptic path parameterized by the 
    radii R1 (real axis) and R2. Uses the quadpy library to compute the integral
    with an adaptive quadrature approach. Returns the integral value as well as 
    an error estimate.'''
    
    val = quadpy.line_segment.integrate_adaptive(lambda t: integrand_ellipse_adaptive(t,a,b,R1,R2,gamma_ellipse,gamma_der_ellipse),[alpha,beta],error)
    return ((1/(2*math.pi*1j)) * val[0],val[1])


def count_solutions_shortest_path(a,b,start_radius,radius_delta,N,error=1):
    '''Solves the shortest path given a,b and outputs the integral, moves along
    arcs. The minimum radius is determined automatically given the optimal radius 
    procedure.'''
    
    start_radius = round(start_radius,15)
    epsilon = round((2*math.pi)/(720),15)
    
    path_start = (start_radius,epsilon)
    
    angles = [epsilon]
    current = epsilon
    delta = round((2*math.pi-2*epsilon)/(N),15)
    for i in range(0,N):
        current+= delta
        angles.append(current)
    
    radii = []
    r = start_radius
    while r<0.9999:
        radii.append(r)
        r+= radius_delta
        
    
    path_end = (start_radius,angles[-1])
    
    nodes = generate_graph(radii,angles)
    graph_ = Graph()
    graph = fill_edges(graph_,radii,angles,nodes,a,b)
    
    cost_func = lambda u,v,e,prev_e: e['cost']
    path_nodes = list(find_path(graph,nodes[path_start],nodes[path_end],cost_func=cost_func).nodes)
    edge_labels = {y:x for x,y in nodes.items()}
    
    
    
    integral = 0
    
    currentDT = datetime.datetime.now()
    for n in range(1,len(path_nodes)):
        r1=edge_labels[path_nodes[n-1]][0]
        r2=edge_labels[path_nodes[n]][0]
        
        theta1 = edge_labels[path_nodes[n-1]][1]
        theta2 = edge_labels[path_nodes[n]][1]
        
        #if there is no radius adjustment, move along an arc
        if r1 == r2:
            
            integral += angle_adjustment_integral(a,b,theta1,theta2,r1,gamma_circle_rad,gamma_der_circle_rad,error)[0]
            
        else:
            z1 = convert(r1,theta1)
            z2 = convert(r2,theta2)    
            x1=z1.real
            y1=z1.imag
            x2=z2.real
            y2=z2.imag
            integral += line_integral(a,b,x1,y1,x2,y2,gamma_radial,gamma_der_radial,error)[0]
    
    #add final integral between start and end node
    integral += angle_adjustment_integral(a,b,0,angles[0],radii[0],gamma_circle_rad,gamma_der_circle_rad,error)[0]
    integral += angle_adjustment_integral(a,b,angles[-1],2*math.pi,radii[0],gamma_circle_rad,gamma_der_circle_rad,error)[0]
    
    currentDT2 = datetime.datetime.now()
    timed = currentDT2 - currentDT
    t = timed.seconds + (timed.microseconds * 1e-6)
    
    
    return integral,t



