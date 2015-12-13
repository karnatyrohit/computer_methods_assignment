from scipy import *
import math
from matplotlib.pyplot import *
from scipy.special import *

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

def poisson1():
	value = - log(rand0())
	return value

	
x = []

for j in range(100):
	t = 0
	for i in range(30):
		t += (poisson1())
	x.append(t)

#hist(x, bins=50)
figure("erlang distibution",figsize=(18,10))
n,bins,patches=hist(x, 30, normed=1,facecolor='red',alpha=0.5,label='Histogram')
shape = 30
lam = 1
y = bins**(shape - 1) * exp(-bins * lam) * lam**(shape)  / gamma(shape)
plot(bins, y,'b--',label = "erlang distribution" )
legend()
savefig("erlang_10000_30.png")

show()
