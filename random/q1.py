from scipy import *
import math
from matplotlib.pyplot import *

A = 16807
MAX = 2147483647
MASK = 12345987
idum = 21314
iset = 0
gauss2 = 0

def rand0():
	global idum
	idum ^= MASK
	value = A*idum
	value = value % MAX
	idum = value
	idum ^= MASK
	return float(value) / MAX

def gauss1():
	rsq=0.0
	global gauss2 
	global iset 
	if(iset==0):
		while(rsq>=1.0 or rsq ==0):
			v1 = 2*rand0() - 1
			v2 = 2*rand0() - 1
			rsq = v1*v1 + v2*v2
		fac = math.sqrt(-2*log(rsq)/rsq)
		gauss2 = v2*fac
		iset = 1
		return v1*fac
	else:
		iset = 0
		return gauss2
if "__main__" == __name__:

	x = []
	for i in range(20000):
		x.append(gauss1())

#hist(x, bins=50)
	n,bins,patches=hist(x, 50, normed=1,facecolor='green',alpha=0.75,label='Histogram')
	y=mlab.normpdf(bins,0,1)
	plot(bins,y,'r--',linewidth=1,label='Standard Normal Distribution')
	print n ,"A" , bins, "B", patches 
	show()	
