import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat

def fit(f,x,y,guess=None):
	params,covariance=curve_fit(f,x,y,p0=guess)
	print(covariance)
	errors=np.sqrt(np.diag(covariance))
	a=ufloat(params[0],errors[0])
	#b=ufloat(params[1],errors[1])
	return a


def gerade(x,m,b):
	return m*x

def wurzel(x,m,b):
	return m*np.sqrt(x)

plt.clf()

N,r2,sigmar2=np.loadtxt("EndToEnd.txt",unpack=True)

plt.errorbar(N,r2,marker="+",linestyle="none",label="Messwerte")
t=np.linspace(N[0],N[N.size-1],1000)
#plt.plot(t,1*t,label=r"$N^1$")
#plt.plot(t,1*t**(6./5.),label=r"$N^{6/5}$")
#plt.plot(t,.01*t**(2*6./5.),label=r"$N^{2*6/5}$")
#plt.plot(t,.05*t**(2),label=r"$N^{2}$")

plt.legend(loc="best")

#plt.xscale("log") 
#plt.yscale("log") 

plt.xlabel(r"$t$")
plt.ylabel(r"$<R^2>$")

plt.savefig("PlotEndToEnd.pdf")
