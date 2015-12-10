from scipy import *
from matplotlib.pyplot import *
from q1 import *

def f(x):
	y = zeros(x.shape)
	for k in range(1,6):
		y += sin((2*k + 1) * x) / (2*k + 1) 
	return y

x = linspace(-pi,pi,101)
#y = [f(x1) for x1 in x]
y = f(x)
figure(3)
plot(x, y)
max_exact=zeros(21)
max_est=zeros(21)
for n in range(2,21):
	#n=5 # order of interpolation
#xarr=linspace(0,1,30)
#yarr=sin(xarr+xarr*xarr);
#t=linspace(0,pi,111)
	delta = 0
	xx=linspace(-pi- delta,pi + delta,999)
#xx=array([-pi])
	z=polint(x,y,xx,n)
	yy=z[0];dyy=z[1]
	y0=f(xx)
	figure(0)
	plot(xx,yy,'ro',xx,y0,'k')
	title("Interpolation by %dth order polynomial" % n)
	figure(1)
	semilogy(xx,abs(yy-y0),'ro',xx,abs(dyy),'k')
	title("Error in interpolation")
	legend(["Actual error","Error Est"])
	max_exact[n]=amax(abs(yy-y0))
	max_est[n]=amax(abs(dyy))
	figure("err")
	semilogy(n,max_exact[n],'ro',n,max_est[n],'bo')
	legend(["Actual error","Error Est"])
show()
