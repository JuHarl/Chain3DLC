from scipy.optimize import curve_fit
from scipy.stats import sem
from scipy import stats
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from pylab import *
from uncertainties import ufloat
from mpl_toolkits.mplot3d import Axes3D

plt.clf()

t,i1,i2=np.loadtxt("heat01.txt", unpack=True)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.scatter(t, i1, i2)
ax.set_xlabel('t')
ax.set_ylabel('i1')
ax.set_zlabel('i2')
plt.show()

#plt.savefig("heat01.pdf")
