import math
import gmpy2
from gmpy2 import mpc
import quadpy
from integrand_functions import integrand_adaptive, integrand_adaptive_line

def H_mpc(a,z):
    '''Returns generating function value. Takes in the a vector, the current position 
    on the parametrized variable t and the gamma function. Uses MPC and is
    not vectorized.'''
    temp=1
    n=len(a)
    for k in range(0,n):
        temp = gmpy2.div(temp, mpc((1-((z))**a[k])))
    return temp

def gamma_circle_rad(t,R):
    '''Circle gamma function using sin and cosine'''
    return R * math.cos(t) + R*1j*math.sin(t)
    
def gamma_der_circle_rad(t,R):
    '''Derivative of above'''
    return -R * math.sin(t) + 1j * R * math.cos(t)

def gamma_radial(R,x0,x1,y0,y1):
    '''Linear Parametrization for R between 0 and 1'''
    return complex(x0+R*(x1-x0),y0+R*(y1-y0))

def gamma_der_radial(R,x0,x1,y0,y1):
    return complex((x1-x0),(y1-y0))

def f_(z,a,b):
    '''Returns the function value at z using MPC'''
    return gmpy2.div(mpc(H_mpc(a,z)),mpc((z**(b+1))))

def convert(R,angle):
    '''Convert radians to complex plane coordinates'''
    return complex(R*math.cos(angle),R*math.sin(angle))


def dist(p1,p2):
    '''Returns the euclidean distance between the two imaginary numbers in the plane'''
    return math.sqrt((p1.real-p2.real)**2 + (p1.imag - p2.imag)**2)

def dist2(p1,p2):
    '''Returns the arc distance between the two imaginary numbers in the plane'''
    r = math.sqrt(p1.real**2 +p1.imag**2)
    theta1 = math.atan(p1.imag/p1.real)
    theta2 = math.atan(p2.imag/p2.real)
    return r*abs(theta2-theta1)
    


def edge_value(p1,p2,a,b):
    '''Determines the cost of the edge between points p1 and p2 given a,b'''    
    return abs(0.5*(f_(p1,a,b) + f_(p2,a,b))*dist(p2,p1))


def edge_value2(p1,p3,p2,a,b):
    '''Determines the cost of the edge between points p1 and p2 given a,b using the three point rule'''    
    return (abs(0.5*(f_(p1,a,b) + f_(p3,a,b))*dist2(p3,p1)) + abs(0.5*(f_(p3,a,b) + f_(p2,a,b))*dist2(p3,p2))) 


def generate_graph(radii,angles):
    '''Generates the graph on which the shortest path algorithm can be applied'''
    nodes = {}
    node_counter=1
    
    for R in radii:
        for angle in angles:
            nodes[R,angle] = node_counter  
            node_counter+=1
    return nodes


def fill_edges(graph,radii,angles,nodes,a,b):
    '''Fills the edges of the given graph with weights using the three point rule and adds them to the Dijkstar library'''

    for r_index in range(0,len(radii)-1):
        for angle_index in range(0,len(angles)-1):
            
            R = radii[r_index]
            R2 = radii[r_index+1]
            
            angle = angles[angle_index]
            angle2 = angles[angle_index+1]
           
            p1 = (R, angle)
            p2 = (R2, angle)
            p3 = (R, angle2)
            
            z1 = convert(R,angle)
            z2 = convert(R2,angle)
            z3 = convert(R,angle2)
            z4 = convert(R,0.5*(angle+angle2))
            
            v1 = edge_value(z1,z2,a,b)
            v2 = edge_value2(z1,z4,z3,a,b)
            
            graph.add_edge(nodes[p1],nodes[p2],{'cost':v1})
            graph.add_edge(nodes[p2],nodes[p1],{'cost':v1})
            graph.add_edge(nodes[p1],nodes[p3],{'cost':v2})
            
    angle = angles[-1]
    for r_index in range(0,len(radii)-1):
        R = radii[r_index]
        R2 = radii[r_index+1]
        
        p1 = (R, angle)
        p2 = (R2, angle)
            
        z1 = convert(R,angle)
        z2 = convert(R2,angle)
            
        v1 = edge_value(z1,z2,a,b)
        
        graph.add_edge(nodes[p1],nodes[p2],{'cost':v1})
        graph.add_edge(nodes[p2],nodes[p1],{'cost':v1})
        
    R = radii[-1]
    for angle_index in range(0,len(angles)-1):
        angle = angles[angle_index]
        angle2 = angles[angle_index+1]
        
        p1 = (R, angle)
        p2 = (R, angle2)
        
        z1 = convert(R,angle)
        z2 = convert(R,angle2)
        z3 = convert(R,0.5*(angle+angle2))
        
        v1 = edge_value2(z1,z3,z2,a,b)
        graph.add_edge(nodes[p1],nodes[p2],{'cost':v1})
        
        
    return graph

def angle_adjustment_integral(a,b,alpha,beta,R,gamma,gamma_der,error):
    '''Adaptive integral from theta1 to theta2 along a constant radius'''
    val=quadpy.line_segment.integrate_adaptive(lambda t: integrand_adaptive(t,a,b,R,gamma,gamma_der),[alpha,beta],error)
    return ((1/(2*math.pi*1j)) * val[0],val[1])


def line_integral(a,b,x1,y1,x2,y2,gamma,gamma_der,error):
    '''Adaptive integral from r1 to r2 along the same radial line'''
    val=quadpy.line_segment.integrate_adaptive(lambda t: integrand_adaptive_line(t,a,b,x1,y1,x2,y2,gamma,gamma_der),[0,1],error)
    return ((1/(2*math.pi*1j)) * val[0],val[1])