from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pylab import *
from uncertainties import ufloat
from uncertainties import unumpy as unp

def fit(f,x,y,guess=None):
	params,covariance=curve_fit(f,x,y,p0=guess)
	#print(covariance)
	errors=np.sqrt(np.diag(covariance))
	a=ufloat(params[0],errors[0])
	#b=ufloat(params[1],errors[1])
	return a

def leastsqerror(f,m,x,y):
	Delta=0
	for i in range(x.size):
		Delta=Delta+(f(x[i],m)-y[i])**2
	return 1.*Delta/x.size**3
		


def gerade(x,b):
	return x+b

def wurzel(x,b):
	return .5*x+b

def hocheinviertel(x,b):
	return .25*x+b


N,M,V,rho,l,K,Beta,x0,tMax,Schrittweite=np.loadtxt("PlotKonfig.txt")


plt.clf()

deltaT,r2,sigmar2=np.loadtxt("g1.txt",unpack=True)

#plt.plot(deltaT,r2,"bx",label="Messwerte")
plt.errorbar(deltaT,r2,yerr=sigmar2,fmt="k+",linestyle="none")


t=np.linspace(deltaT[0],deltaT[deltaT.size-1],1000)

#m1=fit(wurzel,deltaT[0:10],r2[0:10])
#m2=fit(gerade,deltaT[35:],r2[35:])

minStart=0
maxStart=deltaT.size-5
#minEnd=30
maxEnd=deltaT.size


minFehler=1000000000

for i in range(minStart,maxStart):
	for j in range(i+5,maxEnd):
		m2=fit(hocheinviertel,log(deltaT[i:j]),log(r2[i:j]))
		Delta=leastsqerror(hocheinviertel,m2.n,log(deltaT[i:j]),log(r2[i:j]))
		if(Delta<minFehler):
			minFehler=Delta
			Start=i
			End=j

print("m2:")
print(Start)
print(End)

ViertelStart=Start
ViertelEnd=End

m2=fit(hocheinviertel,log(deltaT[Start:End]),log(r2[Start:End]))


plt.plot(deltaT[Start-10:End+10],np.exp(m2.n)*(deltaT[Start-10:End+10])**.25,"b--", linewidth=0.6)



minStart=ViertelEnd
maxStart=deltaT.size-5
maxEnd=deltaT.size


minFehler=100000000

for i in range(minStart,maxStart):
	for j in range(i+5,maxEnd):
		m4=fit(gerade,log(deltaT[i:j]),log(r2[i:j]))
		Delta=leastsqerror(gerade,m4.n,log(deltaT[i:j]),log(r2[i:j]))
		if(Delta<minFehler):
			minFehler=Delta
			Start=i
			End=j
			#print(minFehler)

print("m4:")
print(Start)
print(End)

GeradeStart=Start
GeradeEnd=End

m4=fit(gerade,log(deltaT[Start:End]),log(r2[Start:End]))	


plt.plot(deltaT[Start-10:maxEnd],np.exp(m4.n)*deltaT[Start-10:maxEnd],"r--", linewidth=0.6)


minStart=0
maxStart=ViertelStart-5
maxEnd=ViertelStart


minFehler=100000000

for i in range(minStart,maxStart):
	for j in range(i+5,maxEnd):
		m1=fit(wurzel,log(deltaT[i:j]),log(r2[i:j]))
		Delta=leastsqerror(wurzel,m1.n,log(deltaT[i:j]),log(r2[i:j]))
		if(Delta<minFehler):
			minFehler=Delta
			Start=i
			End=j
			#print(i," ", j, " Fehler= ", minFehler)
print("m1:")
print(Start)
print(End)

m1Start=Start
m1End=End

m1=fit(wurzel,log(deltaT[Start:End]),log(r2[Start:End]))	
		




plt.plot(deltaT[1:End+10],np.exp(m1.n)*(deltaT[1:End+10])**.5,"g--", linewidth=0.6)



minStart=ViertelEnd
maxStart=GeradeStart-5
maxEnd=GeradeStart


minFehler=100000000

for i in range(minStart,maxStart):
	for j in range(i+5,maxEnd):
		m3=fit(wurzel,log(deltaT[i:j]),log(r2[i:j]))
		Delta=leastsqerror(wurzel,m3.n,log(deltaT[i:j]),log(r2[i:j]))
		if(Delta<minFehler):
			minFehler=Delta
			Start=i
			End=j

print("m3:")
print(Start)
print(End)

m3=fit(wurzel,log(deltaT[Start:End]),log(r2[Start:End]))	


plt.plot(deltaT[Start-10:End+10],np.exp(m3.n)*(deltaT[Start-10:End+10])**.5,"g--", linewidth=0.6)

plt.plot(deltaT[m1Start:m1End],np.exp(m1.n)*(deltaT[m1Start:m1End])**.5,"g-")
plt.plot(deltaT[ViertelStart:ViertelEnd],np.exp(m2.n)*(deltaT[ViertelStart:ViertelEnd])**.25,"b-", label=r"$\Delta t^{1/4}$")
plt.plot(deltaT[Start:End],np.exp(m3.n)*(deltaT[Start:End])**.5,"g-", label=r"$\Delta t^{1/2}$")
plt.plot(deltaT[GeradeStart:GeradeEnd],np.exp(m4.n)*deltaT[GeradeStart:GeradeEnd],"r-", label=r"$\Delta t^{1}$")

TauD=(unp.exp(m3-m4))**2

print(TauD)

plt.plot([TauD.n,TauD.n],[0.1,1e5],"m--",label=r"$\tau_D="+"{:.1f}".format(TauD.n) +"\pm"+"{:.1f}".format(TauD.s)+"$")

TauR=(unp.exp(m2-m3))**4
plt.plot([TauR.n,TauR.n],[0.1,1e5],"y--",label=r"$\tau_R="+"{:.1f}".format(TauR.n) +"\pm"+"{:.1f}".format(TauR.s)+"$")


#plt.plot(t,.8*t**.5, label=r"$\Delta T^{ 1/2}$")
#plt.plot(t,0.01*(t)**1, label=r"$\Delta T^{ 1}$")
#plt.plot(t,3.5*(t)**.25, label=r"$\Delta T^{1/4}$")



#print(m1.n)
#print(m2.n)

#tR=(m1/m2)**2

#print("t_R=",tR.n,"+-",tR.s)

plt.ylim(1e-1,1e5)


plt.legend(loc="upper left")

plt.xscale("log") 
plt.yscale("log") 

plt.xlabel(r"$\Delta t$")
plt.ylabel(r"$g_1(t)$")


plt.title(r"$\rho="+str(rho)+r" \, N="+str(int(N))+"\, V="+str(int(V))+"^3$")

plt.savefig("Plot.pdf")
