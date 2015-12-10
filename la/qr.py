from matplotlib.pylab import *
from scipy import *
from numpy import *
import math

x=linspace(0,pi,9)
x=x[0:-1]
A=c_[ones_like(x),cos(x),cos(2*x),cos(3*x),cos(4*x), \
cos(5*x),cos(6*x),cos(7*x)]
B=array([(15, -3, -5, -7),(-3.0, 25, -4, -6), \
(-5.0, -4, 10, -5),(-7.0, -6, -5, 20)])
c=(1+1j)/sqrt(2.0)
d=(1-1j)/sqrt(2.0)
C=array([(1, 1, 1, 1),(1, 1, -1j, d),(1, 1j, -1, -1j), \
(1, c, 1j, 1)])

def householder(a): 
    n = len(a)
    for k in range(n-2):
        u = a[k+1:n,k]
        alpha = math.sqrt(dot(u,u))
        if u[0] < 0.0: alpha = -alpha
        u[0] = u[0] + alpha
        h = dot(u,u)/2.0
        p = dot(a[k+1:n,k+1:n],u)/h
        k1 = dot(u,p)/(2.0*h)
        q = p - k1*u
	print u
	print "as"
	print q
	print " "
	print a , n , k
        a[k+1:n,k+1:n] = a[k+1:n,k+1:n] - outer(q,u) - outer(u,q)
	print "wer"
        a[k,k+1] = alpha
    return np.diagonal(a),np.diagonal(a,1)

def computeP(a): 
    n = len(a)
    p = np.identity(n)*1.0
    for k in range(n-2):
        u = a[k+1:n,k]
        h = np.dot(u,u)/2.0
        v = np.dot(p[1:n,k+1:n],u)/h           
        p[1:n,k+1:n] = p[1:n,k+1:n] - np.outer(v,u)
    return p
      
k=householder(B)
P=computeP(B)
print P
print B*P
print P*B
print P*B*P
