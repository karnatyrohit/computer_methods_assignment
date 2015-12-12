from scipy import *
from scipy import integrate
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
		k = P([1])
		return k / math.sqrt(dotP(k, k))
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

def get_basis(f, d):
	a = [ (integrate.quad(lambda x : f(x)*i(x), 0, 1))[0]  for i in d]
	return a

def get_function(a, d):
	k = P([0])
	for i in range(len(d)):
		k += a[i]*d[i]	
	return k

def sin_taylor(x):
	return 1 - 0.5*power(pi*(x-0.5), 2) + 1/4/3/2*power(pi*(x-0.5), 4)

d1 = get_orthogonal(dotP_type1, 5)
a1 = get_basis(lambda x : sin(pi*x), d1)
e1 = get_function(a1, d1)

d2 = get_orthogonal(dotP_type2, 5)
a2 = get_basis(lambda x : sin(pi*x), d2)
e2 = get_function(a2, d2)

x = linspace(0, 1, 100)
figure("test", figsize=(18,10))
plot(x, e1(x),'ro', label="approx type 1")
plot(x, e2(x),'go', label="approx type 2")
plot(x, sin(pi*x), label="sin")
plot(x, sin_taylor(x),'bo', label="taylor")
legend()
savefig("approx.png")

figure ("error", figsize=(18,10))
loglog(abs(x-0.5), abs(sin(pi*x)-e1(x)), label="orthogonal fit")
loglog(abs(x-0.5), abs(sin(pi*x)-sin_taylor(x)), label="taylor fit")
legend()
savefig("error_log_log.png")
#show()
