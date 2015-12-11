from scipy import *
from matplotlib.pyplot import *
from numpy.polynomial.polynomial import Polynomial as P
from numpy import *
import math
import sys

def dotP_type1(a, b):
	value=0
	for i in range(a.degree()+1):
		for j in range(b.degree()+1):
			value += a.coef[i] * b.coef[j] / (i + j + 1)
	return value


## _get_one_orthogonal(dotP_type1, 3, [P([0,0,1]), P([0,1]), P([0,0,0,1]))
def _get_one_orthogonal(dotP, deg, orthoSet):
	if(deg == 0):
		return P([1])
	else:
		k = P(append(zeros(deg), 1))
		for j in orthoSet:
			k = k - dotP(k, j) * j
		k = k/math.sqrt(dotP(k, k))
		return k 

def get_orthogonal(dotP, deg):
	p = []
	for i in range(deg+1):
		p.append(_get_one_orthogonal(dotP, i, p))
	return p


# test function
a=P([0,1,1])
b=P([0,0,1])
print "Test dotP_type1: ", dotP_type1(a,b)
c = (a*b).integ()
print "Test mult check for type1: ", c, c(1) - c(0)

#plotting functions

d=get_orthogonal(dotP_type1, 5)
x=linspace(0,1,1000)

fig = figure("bases from 1 to 5", figsize=(18,10))
for i in d:
	print i
	plot(x,i(x),label = str(i))
legend()
savefig('orthoPlot_type1.png', )
#show()
for i in d:
	for j in d:
		 sys.stdout.write(" %f " % dotP_type1(i, j))
	print ""

for i in d:
	print dotP_type1(i,i)

