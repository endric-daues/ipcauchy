import numpy as np
from helpers import H,H_ellipse, H_radial

def integrand_adaptive(t,a,b,R,gamma,gamma_der):
    '''Vectorized integrand function returns integrand value along circular path'''
    input_shape=t.shape
    sol=[]
    length = input_shape[0]
    for i in range(0,length):
        sub_sol=[]
        input_sub_shape=t[i].shape
        #sub_length = len(t[i])
        for k in t[i]:
            sub_sol.append(H(a,k,R,gamma) * (1/ (gamma(k,R)**(b+1))) * gamma_der(k,R))
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)
    return sol_array


def integrand_ellipse_adaptive(t,a,b,r1,r2,gamma,gamma_der):
    '''Vectorized integrand function returns integrand value along elliptical path'''
    input_shape=t.shape
    sol=[]
    length = input_shape[0]
    for i in range(0,length):
        sub_sol=[]
        input_sub_shape=t[i].shape
        #sub_length = len(t[i])
        for k in t[i]:
            sub_sol.append((H_ellipse(a,k,r1,r2,gamma) * (1/ (gamma(k,r1,r2)**(b+1))) * gamma_der(k,r1,r2)))
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)
    return sol_array

def integrand_adaptive_line(t,a,b,x1,y1,x2,y2,gamma,gamma_der):
    '''Along a line between two points With MPC'''
    input_shape=t.shape
    sol=[]
    length = input_shape[0]
    for i in range(0,length):
        sub_sol=[]
        input_sub_shape=t[i].shape
        #sub_length = len(t[i])
        for k in t[i]:
            sub_sol.append(H_radial(a,k,x1,x2,y1,y2,gamma) * (1/ (gamma(k,x1,x2,y1,y2)**(b+1))) * gamma_der(k,x1,x2,y1,y2))
        sub_array = np.array(sub_sol)
        sub_array.reshape(input_sub_shape)
        sol.append(sub_array)
    sol_array = np.array(sol)
    sol_array.reshape(input_shape)
    return sol_array
