from scipy import *
from matplotlib.pyplot import *
from numpy.polynomial.polynomial import Polynomial as P
from numpy import *
import math
import sys

def dotP_type0(a, b, t):
	value=0
	for i in range(a.degree()+1):
		for j in range(b.degree()+1):
			value += a.coef[i] * b.coef[j] / (i + j + 1 + t)
	return value

def dotP_type1(a, b):
	return dotP_type0(a, b, 0)

def dotP_type2(a, b):
	return dotP_type0(a, b, 1)


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

print "Test dotP_type2: ", dotP_type2(a,b)
c = (a*b*P([0, 1])).integ()
print "Test mult check for type2: ", c, c(1) - c(0)

#plotting functions
def plot_orthonormal_for_5deg(dotP, fname, figname):
	d=get_orthogonal(dotP, 5)
	x=linspace(0,1,1000)

	fig = figure(figname, figsize=(18,10))
	for i in d:
		print i
		plot(x,i(x),label = str(i))
		legend()
	savefig(fname)
#show()
	for i in d:
		for j in d:
			sys.stdout.write(" %f " % dotP(i, j))
		print ""
plot_orthonormal_for_5deg(dotP_type1,'orthoPlot_type1.png', "bases from 1 to 5 type 1")
plot_orthonormal_for_5deg(dotP_type2,'orthoPlot_type2.png', "bases from 1 to 5 type 2")

