from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pylab import *
from uncertainties import ufloat

def fit(f,x,y,guess=None):
	params,covariance=curve_fit(f,x,y,p0=guess)
	print(covariance)
	errors=np.sqrt(np.diag(covariance))
	a=ufloat(params[0],errors[0])
	#b=ufloat(params[1],errors[1])
	return a


def gerade(x,m):
	return m*x/x

def wurzel(x,m):
	return m*np.sqrt(x)

def hocheinviertel(x,m):
	return m*x**.25

plt.clf()

N,M,V,rho,l,K,Beta,x0,tMax,Schrittweite=np.loadtxt("PlotKonfig.txt")

deltaT,g1,sigmar2=np.loadtxt("g1.txt",unpack=True)
plt.errorbar(deltaT,g1,marker="+",linestyle="none",label=r"$g_1(t)$")


deltaT,g2,sigmar2=np.loadtxt("g2.txt",unpack=True)
plt.errorbar(deltaT,g2,marker="+",linestyle="none",label=r"$g_2(t)$")

deltaT,g3,sigmar2=np.loadtxt("g3.txt",unpack=True)
plt.errorbar(deltaT,g3,marker="+",linestyle="none",label=r"$g_3(t)$")

deltaT,g4,sigmar2=np.loadtxt("g4.txt",unpack=True)
plt.errorbar(deltaT,g4/g1,marker="+",linestyle="none",label=r"$g_4/g_1(t)$")

deltaT,g5,sigmar2=np.loadtxt("g5.txt",unpack=True)
plt.errorbar(deltaT,g5,marker="+",linestyle="none",label=r"$g_5(t)$")

deltaT,g6,sigmar2=np.loadtxt("g6.txt",unpack=True)
plt.errorbar(deltaT,g6,marker="+",linestyle="none",label=r"$g_6(t)$")

t=np.linspace(deltaT[0],deltaT[deltaT.size-1],1000)

#m1=fit(wurzel,deltaT[0:10],r2[0:10])
#m2=fit(gerade,deltaT[35:],r2[35:])

m1=fit(wurzel,deltaT[5:20],g2[5:20])

plt.plot(deltaT[5:20],wurzel(deltaT[5:20],m1.n),"g-")#, label=r"$\Delta t^{1/2}$")
plt.plot(deltaT[1:26],wurzel(deltaT[1:26],m1.n),"g--")

m2=fit(hocheinviertel,deltaT[20:33],g2[20:33])

plt.plot(deltaT[20:33],hocheinviertel(deltaT[20:33],m2.n),"b-")#, label=r"$\Delta t^{1/4}$")


m3=fit(wurzel,deltaT[39:43],g2[39:43])

plt.plot(deltaT[39:43],wurzel(deltaT[39:43],m3.n),"g-")#, label=r"$\Delta t^{1/2}$")
plt.plot(deltaT[20:50],wurzel(deltaT[20:50],m3.n),"g--")

m4=fit(gerade,deltaT[47:64],g2[47:64])

plt.plot(deltaT[47:64],gerade(deltaT[47:64],m4.n),"r-")#, label=r"$\Delta t^{1}$")
plt.plot(deltaT[40:70],gerade(deltaT[40:70],m4.n),"r--")

TauD=(m4/m3)**2



print(TauD)

plt.plot([TauD.n,TauD.n],[0.1,1e4],"--",label=r"$\tau_D="+"{:.1f}".format(TauD.n) +"\pm"+"{:.1f}".format(TauD.s)+"$")

TauR=(m2/m3)**4
plt.plot([TauR.n,TauR.n],[0.1,1e4],"--",label=r"$\tau_R="+"{:.1f}".format(TauR.n) +"\pm"+"{:.1f}".format(TauR.s)+"$")


#plt.plot(t,1.2*t**.5)#, label=r"$\Delta T^{ 1/2}$")
#plt.plot(t,0.02*(t)**1)#, label=r"$\Delta T^{ 1}$")
#plt.plot(t,.8*(t)**.25)#, label=r"$\Delta T^{1/4}$")

#print(m1.n)
#print(m2.n)

#tR=(m1/m2)**2

#print("t_R=",tR.n,"+-",tR.s)

plt.ylim(1e-1,1e4)


plt.legend(loc="best")

plt.xscale("log") 
plt.yscale("log") 

plt.xlabel(r"$\Delta t$")
plt.ylabel(r"$g_1(t)$")

plt.title(r"$\rho="+str(rho)+" \, N="+str(int(N))+"\, V="+str(int(V))+"^3$")

plt.savefig("PlotGs.pdf")
