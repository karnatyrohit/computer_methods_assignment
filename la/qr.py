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
    print "A", a
    n = len(a)
    for k in range(n-2):
        u = (a[k+1:n,k]).copy()
        alpha = math.sqrt(dot(u,u))
        if u[0] < 0.0: alpha = -alpha
        u[0] = u[0] + alpha
        h = dot(u,u)/2.0
        p = dot(a[k+1:n,k+1:n],u)/h
        k1 = dot(u,p)/(2.0*h)
        q = p - dot(k1,u)
	print u
	print "as"
	print q
	print " "
	print a , n , k ,k1
        a[k+1:n,k+1:n] = a[k+1:n,k+1:n] - outer(q,u) - outer(u,q)
	print "wer"
        a[k,k+1] = -alpha
    return np.diagonal(a),np.diagonal(a,1)

def computeP(a): 
    print "B", a
    n = len(a)
    p = np.identity(n)
    for k in range(n-2):
    	p = np.identity(n)
        u = (a[k+1:n,k]).copy()
        alpha = math.sqrt(dot(u,u))
        if u[0] < 0.0: alpha = -alpha
        u[0] = u[0] + alpha
        h = np.dot(u,u)/2.0
        v = u.transpose()           
        p[k+1:n,k+1:n] = p[k+1:n,k+1:n] - np.outer(u,v)/h

	#print p
	a= dot(dot(p,a),p)
	print around(a,decimals=6)
	#print a
    return a
      
k=householder(A.copy())
print "====="
print k
B=computeP(A.copy())
def qrStep(b):
	n = len(b)
	q=identity(n)
	for k in range(n-1):
		p = np.identity(n)  
		u = (b[k:n,k]).copy() # getU(b, k)
		alpha = math.sqrt(dot(u,u)) 
		if u[0] < 0.0: alpha = -alpha
		u[0] = u[0] + alpha
		h = np.dot(u,u)/2.0
		v = u.transpose()           
		p[k:n,k:n] = p[k:n,k:n] - np.outer(u,v)/h # getP(u, k)
		q=dot(q,p.transpose())
		b = dot(p,b)
	outcome=(2*trace(b,offset=1))/trace(b)
	return (q,b,outcome)

def qrDecomposition(b):
	outcome = 1
	q = None
	r = None
	i=0
	while(abs(outcome) > 1e-10):
		(q,r,outcome) = qrStep(b)
		b=dot(r,q)
		i=i+1
		#print i
	print b
	print around(b,decimals=6)
	return (q,r)

print "A", around(A,decimals=6)
print "B", around(B,decimals=6)
(Q,R) = qrDecomposition(B)
g=linalg.eig(B)
print g[0]
#print around(Q)
#print around(R)
