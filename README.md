![Screenshot](ipcauchy.jpeg)


The code contained here follows the theory and numerical study described in (link to paper). Users can download the directory and use the **example.py** file to count solutions to Knapsack instances with a single constraint using a complex path integral that follows a circular, elliptic or shortest path. The code is written in Python 3.

First, the quadpy numerical integration library must be installed. All other libraries should be included with Anaconda.

```
pip install quadpy
```

The following code is included in the **example.py** file, and can be used to test our implementation on any knapsack constraint. Here, *a* represents the coefficient vector and *b* represents the right-hand side constraint value.

```
from integration_caller_functions import optimal_radius
from integration_caller_functions import count_solutions_circle
from integration_caller_functions import count_solutions_ellipse
from integration_caller_functions import count_solutions_shortest_path
import math

# Initiate a knapsack constraint
a = [2,2,3,4,5]
b = 10

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
```


Finally, we have included a function to retrieve the coefficients and RHS values from the Pisinger instances (http://hjemmesider.diku.dk/~pisinger/codes.html) used in the paper. These can also be used in the above code. Again, this code is set up for use in the **example.py** file.

```
f = open('./smallcoeff_pisinger/knapPI_1_50_1000.csv', "r")
a_,b_= read_pisinger_file(f)

#Select the first Pisinger instance, p1
a = a_[0]
b = b_[0]
#print('a: ',a)
#print('b: ',b)

#Run the shortest path parameterization for this instance
r = optimal_radius(0,a,b)
output_sp = count_solutions_shortest_path(a,b,r,r_,N)
solutions_sp = round(output_sp[0].real,0)
print('SP Solution Count to p1: ', solutions_sp)
```