from integration_caller_functions import optimal_radius
from integration_caller_functions import count_solutions_circle
from integration_caller_functions import count_solutions_ellipse
from integration_caller_functions import count_solutions_shortest_path
from helpers import read_pisinger_file
import math

# Initiate a knapsack constraint
a = [2,2,3,4,5]
b = 10

print('Coefficients: {}'.format(a))
print('b value: {}'.format(b))


#determine the optimal bypass point for this instance, r
r = optimal_radius(0,a,b)

#count solutions using a circular parameterization with radius r
output_circle = count_solutions_circle(a,b,r)
solutions_circle = round(output_circle[0].real,0)
print('Circle Solution Count: ',solutions_circle)

#count solutions using a elliptic parameterization with radii R1 and R2
R1 = r
R2 = optimal_radius(math.pi/2,a,b)
output_ellipse = count_solutions_ellipse(a,b,R1,R2)
solutions_ellipse = round(output_ellipse[0].real,0)
print('Ellipse Solution Count: ',solutions_ellipse)


#count solutions using a shortest parameterization with N angular nodes and a 
#radial distance of r_ and a bypass point of r
N = 36
r_ = 0.001

output_sp = count_solutions_shortest_path(a,b,r,r_,N)
solutions_sp = round(output_sp[0].real,0)
print('SP Solution Count: ', solutions_sp)



#We can also retrieve the Pisinger instances used in the numerical study
f = open('./smallcoeff_pisinger/knapPI_1_50_1000.csv', "r")
a_,b_= read_pisinger_file(f)

#Select the first Pisinger instance, p1
a = a_[0]
b = b_[0]
print('Pisinger Instance ')
#print('a: ',a)
#print('b: ',b)

#Run the shortest path parameterization for this instance
r = optimal_radius(0,a,b)
output_sp = count_solutions_shortest_path(a,b,r,r_,N)
solutions_sp = round(output_sp[0].real,0)
print('SP Solution Count to p1: ', solutions_sp)
print('Integration Time: {}s'.format(round(output_sp[1]),3))
