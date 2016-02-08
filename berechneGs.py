import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat

def g1(Mx,My,Mz,M):
	


Sx,Sy,Sz,Mx,My,Mz,Fx,Fy,Fz,Lx,Ly,Lz=np.loadtxt("MidBeads.txt", unpack=True)

N,M,V,rho,l,K,Beta,x0,tMax,Schrittweite=np.loadtxt("PlotKonfig.txt")

 
