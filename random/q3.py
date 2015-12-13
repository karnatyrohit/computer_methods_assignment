from scipy import *
from matplotlib.pyplot import *
from q1 import *
from numpy import *



def func1(x,y):
	alpha = alpha1(x,y)
	u = x*cos(alpha) + y*sin(alpha) - 0.5
	v = -x*sin(alpha) + y*cos(alpha) - 0.5
	return u*u + v*v

def alpha1(x,y):
	alpha = pi * sin(10*(sqrt(x*x + y*y ) - 0.5))
	return alpha

def contour1():	
	value = func1(x,y)
	return (x, y, value)

R1 = 1 + 1/sqrt(2) 
#x = linspace(-R1,R1,1000)
#y = linspace(-R1,R1,1000)
x = []
y = []
for i in range(1000):
	x.append(2*R1*rand0() - R1)
	y.append(2*R1*rand0() - R1)
x = sort(x)
y = sort(y)

X , Y = meshgrid(x,y)	
Z = func1(X , Y)

contourf(X,Y,Z,[1.0,0.0])
savefig("contour.png")
#show()

def get_random_in_range(xran):
	(mn, mx) = xran
	return (mx - mn) * rand0() + mn

# montecarlo
def montecarlo_area(xran, yran, f, num_steps):
	count = 0
	for i in range(num_steps):
		(x, y) = (get_random_in_range(xran), get_random_in_range(yran))
		val = f(x, y)
		if(abs(val) < 1):
			count += 1
	return float(count) / num_steps * (xran[1] - xran[0])*(yran[1] - yran[0])

for i in [1,10,100,1000,10000]:
	print "Monte carlo area %d: " % (1000*i), montecarlo_area((-R1, R1), (-R1, R1), func1, 1000*i)

