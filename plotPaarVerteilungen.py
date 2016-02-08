from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pylab import *
from uncertainties import ufloat

plt.clf()

N,M,V,rho,l,K,Beta,x0,tMax,Schrittweite=np.loadtxt("PlotKonfig.txt")

r,g=np.loadtxt("paarVerteilungBeads.txt", unpack=True)

plt.plot(r,g,label="g(r) Beads")

r,g=np.loadtxt("paarVerteilungBonds.txt", unpack=True)

plt.plot(r,g,label="g(r) Bonds")

plt.legend(loc="best")

plt.xlabel("r")
plt.ylabel("g(r)")

plt.title(r"$\rho="+str(rho)+" \, N="+str(int(N))+"\, V="+str(int(V))+"^3$")

plt.yscale("log")
plt.grid()

plt.savefig("paarVerteilung.pdf")
