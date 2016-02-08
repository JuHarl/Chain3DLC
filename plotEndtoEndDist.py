from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
from pylab import *
from uncertainties import ufloat

def fit(f,x,y,guess=None):
	params,covariance=curve_fit(f,x,y,p0=guess)
	print(covariance)
	errors=np.sqrt(np.diag(covariance))
	a=ufloat(params[0],errors[0])
	b=ufloat(params[1],errors[1])
	return a,b

def dist(r,C,K):
	return C*r**0.274*np.exp(-(K*r)**2.451)

N,M,V,rho,l,K,Beta,x0,tMax,Schrittweite=np.loadtxt("PlotKonfig.txt")

plt.clf()
x1=np.loadtxt("EndToEnd.txt",unpack=True)
R=np.mean(x1)
RN=np.sqrt(np.mean(x1))


b=np.linspace(0,RN*3,50)
data=hist(np.sqrt(x1),bins=b,label="EndToEnd",normed=True)
#	print(data)

x=np.linspace(0,RN*3,1000)
t=2.451
theta=0.274
g=0.279
gamma=1.162
nu=0.592
K=np.sqrt(math.gamma((theta+5)/t)/math.gamma((theta+3)/t))
C=t*math.gamma((theta+5)/t)**((theta+3)/2)/math.gamma((theta+3)/t)**((theta+5)/2)
plt.plot(x,x**2*1/RN**3*C*(x/RN)**theta*np.exp(-(K*x/RN)**t),label="SAW dist")
plt.plot(x,4*np.pi*x**2*(3/(2*np.pi*(N-1)))**(1.5)*np.exp(-3*x**2/(2*(N-1))),label="Distribution ideale Kette")
xlabel(r"$|R|$")
ylabel("Haeufigkeit")
legend(loc="best")

plt.title(r"$N="+str(int(N))+r"\, \rho=" +str(rho)+" <R^2>=" +str(R)+"\, V="+str(int(V))+"^3$")

plt.savefig("EndToEnd.pdf") 
